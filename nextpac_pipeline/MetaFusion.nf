#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------------------
// Parameter Definitions
// ----------------------------
params.seqs         = '/nfs/turbo/lsa-tyleho/mycology/hkz/PebbleScout/SRR9034100/SRR9034100_{1,2}.fastq'
params.outdir       = '/scratch/tyleho_root/tyleho0/MMC/Results'
params.catdb        = '/nfs/turbo/lsa-tyleho/mycology/hkz/CAT_prepare_20210107/2021-01-07_CAT_database'
params.taxdb        = '/nfs/turbo/lsa-tyleho/mycology/hkz/CAT_prepare_20210107/2021-01-07_taxonomy'
params.kaiju_names  = '/scratch/tyleho_root/tyleho0/MMC/kaijudb/names.dmp'
params.kaiju_nodes  = '/scratch/tyleho_root/tyleho0/MMC/kaijudb/nodes.dmp'
params.kaiju_db     = '/scratch/tyleho_root/tyleho0/MMC/kaijudb/kaiju_db_nr_euk.fmi'
params.threads      = 32   // Global thread parameter for processes

// ----------------------------
// Process: Trimmomatic (Quality trimming)
// ----------------------------
process trimmomatic {
    input:
    tuple val(sra_id), path(reads)
    
    output:
    tuple path("${sra_id}_paired_1.fastq"),
          path("${sra_id}_single_1.fastq"),
          path("${sra_id}_paired_2.fastq"),
          path("${sra_id}_single_2.fastq"), emit: trimmed_reads
    
    script:
    """
    trimmomatic PE -phred33 \\
      ${reads[0]} ${reads[1]} \\
      ${sra_id}_paired_1.fastq ${sra_id}_single_1.fastq \\
      ${sra_id}_paired_2.fastq ${sra_id}_single_2.fastq \\
      ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

// ----------------------------
// Process: Assembly (Using metaSPAdes)
// ----------------------------
process assembly {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::spades'
    cpus params.threads

    input:
    // Pass the sample ID along with the trimmed reads tuple
    tuple val(sra_id), path(cleanedReads)
    
    output:
    // Output directory renamed to include the sample ID
    path "spades_out_${sra_id}" into spades_out_ch
    
    script:
    """
    spades.py --meta -k 21,25,31,55,65,75,85,105 -o spades_out_${sra_id} \\
      -1 ${cleanedReads[0]} -2 ${cleanedReads[2]} \\
      -s ${cleanedReads[1]} -s ${cleanedReads[3]} \\
      --threads ${task.cpus}
    # Modify scaffold headers for downstream compatibility
    cat spades_out_${sra_id}/scaffolds.fasta | sed 's/[[:blank:]].*//' > spades_out_${sra_id}/scaffolds.modified.fasta
    """
}

// ----------------------------
// Process: MaxBin (Binning using MaxBin2)
// ----------------------------
process maxbin {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::maxbin2'
    cpus 8

    input:
    // Receive assembly output and trimmed reads (ensuring matching sample IDs)
    tuple val(sra_id), path(spades_dir)
    tuple val(sra_id2), path(cleanedReads)
    
    output:
    path "maxbin_out_${sra_id}" into maxbin_out_ch
    
    script:
    """
    mkdir maxbin_out_${sra_id}
    run_MaxBin.pl -contig ${spades_dir}/scaffolds.modified.fasta \\
      -reads ${cleanedReads[0]} -reads2 ${cleanedReads[2]} \\
      -thread ${task.cpus} -out maxbin_out_${sra_id}/bin
    mkdir maxbin_out_${sra_id}/maxbin_bins
    cp maxbin_out_${sra_id}/*.fasta maxbin_out_${sra_id}/maxbin_bins
    """
}

// ----------------------------
// Process: Minimap (Mapping reads to assembly)
// ----------------------------
process minimap {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::minimap2 bioconda::samtools'
    cpus 12

    input:
    tuple val(sra_id), path(spades_dir)
    tuple val(sra_id2), path(cleanedReads)
    
    output:
    path "minimap_out_${sra_id}" into minimap_out_ch
    
    script:
    """
    mkdir minimap_out_${sra_id}
    minimap2 -ax sr -t ${task.cpus} \\
      ${spades_dir}/scaffolds.modified.fasta \\
      ${cleanedReads[0]} ${cleanedReads[2]} | \\
      samtools sort -@${task.cpus} -O BAM -o minimap_out_${sra_id}/${sra_id}.bam
    """
}

// ----------------------------
// Process: Depth (Contig depth estimation)
// ----------------------------
process depth {
    publishDir params.outdir, mode: 'copy'
    // Ensure jgi_summarize_bam_contig_depths is available (consider adding a conda package if needed)
    cpus 4

    input:
    tuple val(sra_id), path(minimap_dir)
    
    output:
    path "depth_${sra_id}.txt" into depth_ch
    
    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth_${sra_id}.txt \\
      ${minimap_dir}/${sra_id}.bam
    """
}

// ----------------------------
// Process: MetaBAT (Binning using MetaBAT2)
// ----------------------------
process metabat {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::metabat2'
    cpus 4

    input:
    tuple val(sra_id), path(spades_dir)
    path depth_file from depth_ch.filter { it.name.contains(sra_id) }
    
    output:
    path "metabat_out_${sra_id}" into metabat_out_ch
    
    script:
    """
    metabat2 -i ${spades_dir}/scaffolds.modified.fasta \\
      -a ${depth_file} -o metabat_out_${sra_id}/bin
    mkdir metabat_out_${sra_id}/metabat_bins
    mv metabat_out_${sra_id}/*.fa metabat_out_${sra_id}/metabat_bins
    """
}

// ----------------------------
// Process: Concoct (Binning using CONCOCT via metaWRAP)
// ----------------------------
process concoct {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::metawrap'
    cpus 4

    input:
    tuple val(sra_id), path(spades_dir)
    tuple val(sra_id2), path(cleanedReads)
    
    output:
    path "concoct_out_${sra_id}" into concoct_out_ch
    
    script:
    """
    metaWRAP binning --concoct -a ${spades_dir}/scaffolds.modified.fasta \\
      -o concoct_out_${sra_id} \\
      ${cleanedReads[0]} ${cleanedReads[2]}
    rm concoct_out_${sra_id}/concoct_bins/unbinned.fa
    """
}

// ----------------------------
// Process: Table (Generate contig-to-bin mapping)
// ----------------------------
process table {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::das_tool'
    cpus 2

    input:
    // Combine outputs from MaxBin, MetaBAT, and CONCOCT (assumes matching sample IDs)
    tuple val(sra_id), path(maxbin_dir) from maxbin_out_ch
    tuple val(sra_id2), path(metabat_dir) from metabat_out_ch
    tuple val(sra_id3), path(concoct_dir) from concoct_out_ch
    
    output:
    tuple path("maxbin_contig2bin_${sra_id}.tsv"),
          path("metabat_contig2bin_${sra_id}.tsv"),
          path("concoct_contig2bin_${sra_id}.tsv") into table_out_ch
    
    script:
    """
    Fasta_to_Contig2Bin.sh -e fasta -i ${maxbin_dir}/maxbin_bins > maxbin_contig2bin_${sra_id}.tsv
    Fasta_to_Contig2Bin.sh -e fa -i ${metabat_dir}/metabat_bins > metabat_contig2bin_${sra_id}.tsv
    Fasta_to_Contig2Bin.sh -e fa -i ${concoct_dir}/concoct_bins > concoct_contig2bin_${sra_id}.tsv
    """
}

// ----------------------------
// Process: DAS_Tool (Integrate binning results)
// ----------------------------
process dastool {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::das_tool'
    cpus 12

    input:
    tuple val(sra_id), path(table_files) from table_out_ch
    tuple val(sra_id2), path(spades_dir)
    // Use sample ID for naming
    val(sra_id_str) from Channel.value(sra_id)
    
    output:
    path "dastool_out_${sra_id}" into dastool_out_ch
    
    script:
    """
    mkdir dastool_out_${sra_id}
    DAS_Tool --write_bin_evals --write_bins --score_threshold=0.1 \\
      -t ${task.cpus} -i ${table_files} --labels maxbin,metabat,concoct \\
      -c ${spades_dir}/scaffolds.modified.fasta \\
      -o dastool_out_${sra_id}/${sra_id}
    """
}

// ----------------------------
// Process: BAT (Taxonomic annotation for bins using CAT)
// ----------------------------
process BAT {
    publishDir params.outdir, mode: 'copy'
    cpus 4

    input:
    tuple val(sra_id), path(dastool_dir) from dastool_out_ch
    val(sra_id_str) from Channel.value(sra_id)
    
    output:
    path "BAT_out_${sra_id}" into BAT_out_ch
    
    script:
    """
    mkdir BAT_out_${sra_id}
    CAT bins -b ${dastool_dir}/${sra_id}_DASTool_bins \\
      -d ${params.catdb} -t ${params.taxdb} --bin_suffix .fa \\
      -o BAT_out_${sra_id}/${sra_id}
    """
}

// ----------------------------
// Process: CAT (Taxonomic annotation for contigs)
// ----------------------------
process CAT {
    publishDir params.outdir, mode: 'copy'
    cpus 4

    input:
    tuple val(sra_id), path(spades_dir)
    val(sra_id_str) from Channel.value(sra_id)
    
    output:
    path "CAT_out_${sra_id}" into CAT_out_ch
    
    script:
    """
    mkdir CAT_out_${sra_id}
    CAT contigs -c ${spades_dir}/scaffolds.modified.fasta \\
      -d ${params.catdb} -t ${params.taxdb} -o CAT_out_${sra_id}/${sra_id}
    """
}

// ----------------------------
// Process: Adding Names for BAT results
// ----------------------------
process adding_names_BAT {
    publishDir params.outdir, mode: 'copy'
    cpus 2

    input:
    tuple val(sra_id), path(BAT_dir) from BAT_out_ch
    val(sra_id_str) from Channel.value(sra_id)
    
    output:
    path "bins_taxonomy_${sra_id}.txt" into BAT_names_ch
    
    script:
    """
    CAT add_names -i ${BAT_dir}/${sra_id}.bin2classification.txt \\
      -o bins_taxonomy_${sra_id}.txt -t ${params.taxdb} --only_official
    """
}

// ----------------------------
// Process: Adding Names for CAT results
// ----------------------------
process adding_names_CAT {
    publishDir params.outdir, mode: 'copy'
    cpus 2

    input:
    tuple val(sra_id), path(CAT_dir) from CAT_out_ch
    val(sra_id_str) from Channel.value(sra_id)
    
    output:
    path "contigs_taxonomy_${sra_id}.txt" into CAT_names_ch
    
    script:
    """
    CAT add_names -i ${CAT_dir}/${sra_id}.contig2classification.txt \\
      -o contigs_taxonomy_${sra_id}.txt -t ${params.taxdb} --only_official
    """
}

// ----------------------------
// Process: Summarise CAT results
// ----------------------------
process summarize_CAT {
    publishDir params.outdir, mode: 'copy'
    cpus 2

    input:
    tuple val(sra_id), path(spades_dir)
    path contigs_taxonomy from CAT_names_ch.filter { it.name.contains(sra_id) }
    
    output:
    path "CAT_summary_${sra_id}.txt" into CAT_summary_ch
    
    script:
    """
    CAT summarise -c ${spades_dir}/scaffolds.modified.fasta \\
      -i ${contigs_taxonomy} -o CAT_summary_${sra_id}.txt
    """
}

// ----------------------------
// Process: Assembly Statistics
// ----------------------------
process assembly_stat {
    publishDir params.outdir, mode: 'copy'
    cpus 2

    input:
    tuple val(sra_id), path(spades_dir)
    
    output:
    path "assembly.stats_${sra_id}.txt"
    
    script:
    """
    gaas_fasta_statistics.pl -f ${spades_dir}/scaffolds.modified.fasta > assembly.stats_${sra_id}.txt
    """
}

// ----------------------------
// Process: Bins Abundance Estimation
// ----------------------------
process bins_abundance {
    publishDir params.outdir, mode: 'copy'
    cpus 4

    input:
    tuple val(sra_id), path(dastool_dir) from dastool_out_ch
    tuple val(sra_id2), path(spades_dir)
    tuple val(sra_id3), path(cleanedReads)
    
    output:
    path "bins_abundance_${sra_id}" into bins_abundance_ch
    
    script:
    """
    mkdir bins_abundance_${sra_id}
    metaWRAP quant_bins -b ${dastool_dir}/${sra_id}_DASTool_bins \\
      -a ${spades_dir}/scaffolds.modified.fasta \\
      -o bins_abundance_${sra_id} \\
      ${cleanedReads[0]} ${cleanedReads[2]}
    """
}

// ----------------------------
// Process: Metaxa2 (Ribosomal marker analysis)
// ----------------------------
process metaxa {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::metaxa2'
    cpus 4

    input:
    tuple val(sra_id), path(cleanedReads)
    val(sra_id_str) from Channel.value(sra_id)
    
    output:
    path "metaxa_out_${sra_id}" into metaxa_out_ch
    
    script:
    """
    mkdir metaxa_out_${sra_id}
    metaxa2 -1 ${cleanedReads[0]} -2 ${cleanedReads[2]} -E 1e-5 -t all \\
      -o metaxa_out_${sra_id}/${sra_id}
    """
}

// ----------------------------
// Process: ITSx (Fungal ITS extraction)
// ----------------------------
process ITSx {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::itsx'
    cpus 12

    input:
    tuple val(sra_id), path(spades_dir)
    val(sra_id_str) from Channel.value(sra_id)
    
    output:
    path "ITSx_out_${sra_id}" into ITSx_out_ch
    
    script:
    """
    mkdir ITSx_out_${sra_id}
    ITSx -i ${spades_dir}/scaffolds.modified.fasta \\
      -o ITSx_out_${sra_id}/${sra_id} --cpu ${task.cpus}
    """
}

// ----------------------------
// Process: Kaiju (Taxonomic classification of reads)
// ----------------------------
process kaiju {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::kaiju'
    cpus 4

    input:
    tuple val(sra_id), path(cleanedReads)
    val(sra_id_str) from Channel.value(sra_id)
    
    output:
    path "kaiju_out_${sra_id}" into kaiju_out_ch
    
    script:
    """
    mkdir kaiju_out_${sra_id}
    kaiju -E 1e-5 -t ${params.kaiju_nodes} -f ${params.kaiju_db} \\
      -i ${cleanedReads[0]} -j ${cleanedReads[2]} \\
      -o kaiju_out_${sra_id}/${sra_id}_kaiju.out
    """
}

// ----------------------------
// Process: Kaiju Summary (Generate classification table)
// ----------------------------
process kaiju_summary {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::kaiju'
    cpus 2

    input:
    tuple val(sra_id), path(kaiju_dir) from kaiju_out_ch
    val(sra_id_str) from Channel.value(sra_id)
    
    output:
    path "kaiju_summary_${sra_id}.tsv"
    
    script:
    """
    kaiju2table -t ${params.kaiju_nodes} -n ${params.kaiju_names} -r genus \\
      -o kaiju_summary_${sra_id}.tsv \\
      ${kaiju_dir}/${sra_id}_kaiju.out
    """
}

// ----------------------------
// Process: Reassemble (Refinement of bins)
// ----------------------------
process reassemble {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::metawrap'
    cpus 12

    input:
    tuple val(sra_id), path(dastool_dir) from dastool_out_ch
    tuple val(sra_id2), path(cleanedReads)
    
    output:
    path "reassembly_${sra_id}" into reassembly_ch
    
    script:
    """
    metaWRAP reassemble_bins -o reassembly_${sra_id} \\
      -b ${dastool_dir}/${sra_id}_DASTool_bins \\
      -1 ${cleanedReads[0]} -2 ${cleanedReads[2]} -t ${task.cpus} \\
      --skip-checkm --parallel
    mkdir -p reassembly_${sra_id}/reassembled_bins/strict
    mkdir -p reassembly_${sra_id}/reassembled_bins/permissive
    mv reassembly_${sra_id}/reassembled_bins/*.strict.fa reassembly_${sra_id}/reassembled_bins/strict
    mv reassembly_${sra_id}/reassembled_bins/*.permissive.fa reassembly_${sra_id}/reassembled_bins/permissive
    """
}

// ----------------------------
// Workflow Definition
// ----------------------------
workflow {
    // Read input paired files; 'flat: true' ensures proper tuple creation
    reads_ch = Channel.fromFilePairs(params.seqs, flat: true)
    
    // Quality trimming
    trimmed_ch = reads_ch | trimmomatic
    
    // Assembly from trimmed reads
    assembly_ch = trimmed_ch | assembly
    
    // Binning processes
    maxbin_ch   = assembly_ch.combine(trimmed_ch) | maxbin
    minimap_ch  = assembly_ch.combine(trimmed_ch) | minimap
    depth_ch    = minimap_ch.map { spades_dir -> tuple(spades_dir.baseName, spades_dir) } | depth
    metabat_ch  = assembly_ch.combine(depth_ch) | metabat
    concoct_ch  = assembly_ch.combine(trimmed_ch) | concoct
    
    // Generate contig-to-bin mapping table
    table_ch    = maxbin_ch.combine(metabat_ch).combine(concoct_ch) | table
    
    // Integrate binning results with DAS_Tool
    dastool_ch  = table_ch.combine(assembly_ch) | dastool
    
    // Reassembly of bins for refinement
    reassemble_ch = dastool_ch.combine(trimmed_ch) | reassemble
    
    // Taxonomic annotation and naming
    BAT_ch          = dastool_ch | BAT
    CAT_ch          = assembly_ch | CAT
    adding_names_BAT_ch = BAT_ch | adding_names_BAT
    adding_names_CAT_ch = CAT_ch | adding_names_CAT
    summarize_CAT_ch    = assembly_ch.combine(adding_names_CAT_ch) | summarize_CAT
    
    // Assembly statistics
    assembly_stat(assembly_ch)
    
    // Abundance estimation
    bins_abundance_ch = dastool_ch.combine(assembly_ch).combine(trimmed_ch) | bins_abundance
    
    // Additional taxonomic profiling and marker extraction
    metaxa(trimmed_ch)
    ITSx(assembly_ch)
    kaiju_ch   = trimmed_ch | kaiju
    kaiju_summary(kaiju_ch)
}

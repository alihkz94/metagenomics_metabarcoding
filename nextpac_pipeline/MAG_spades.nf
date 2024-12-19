#!/usr/bin/env nextflow

params.seqs = '/nfs/turbo/lsa-tyleho/mycology/hkz/PebbleScout/SRR9034100/SRR9034100_{1,2}.fastq'
params.outdir = '/scratch/tyleho_root/tyleho0/MMC/Results'
params.catdb = '/nfs/turbo/lsa-tyleho/mycology/hkz/CAT_prepare_20210107/2021-01-07_CAT_database'
params.taxdb = '/nfs/turbo/lsa-tyleho/mycology/hkz/CAT_prepare_20210107/2021-01-07_taxonomy'
params.kaiju_names = '/scratch/tyleho_root/tyleho0/MMC/kaijudb/names.dmp'
params.kaiju_nodes = '/scratch/tyleho_root/tyleho0/MMC/kaijudb/nodes.dmp'
params.kaiju_db = '/scratch/tyleho_root/tyleho0/MMC/kaijudb/kaiju_db_nr_euk.fmi'

process trimmomatic {
	publishDir = params.outdir
	conda 'bioconda::trimmomatic'

	input:
	tuple val (SRA_id), path(Reads)
	
	output:
	tuple path(fq_1_paired), path(fq_1_unpaired), path(fq_2_paired), path(fq_2_unpaired)

	script:
	fq_1_paired = SRA_id  + '_paired_1.fastq'
	fq_1_unpaired = SRA_id + '_single_1.fastq'
	fq_2_paired = SRA_id + '_paired_2.fastq'
	fq_2_unpaired = SRA_id + '_single_2.fastq'

	"""
	trimmomatic \
	PE -phred33 \
	${Reads[0]} \
	${Reads[1]} \
	$fq_1_paired \
	$fq_1_unpaired \
	$fq_2_paired \
	$fq_2_unpaired \
	ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
	"""
}

process assembly {
	publishDir = params.outdir
	conda 'bioconda::spades'
	
	input:
	path(cleanedReads)
	
	output:
	path(spades_out) 
	
	"""
	spades.py \
	--meta \
	-k 21,25,31,55,65,75,85,105 \
	-o spades_out \
	-1 ${cleanedReads[0]} \
	-2 ${cleanedReads[2]} \
	-s ${cleanedReads[1]} \
	-s ${cleanedReads[3]} \
	--threads 32
	cat spades_out/scaffolds.fasta | sed 's/[[:blank:]].*//' > spades_out/scaffolds.modified.fasta
	"""
}

process maxbin {
	publishDir = params.outdir
	conda 'bioconda::maxbin2'

	input:
	path (spades_out)
	path (cleanedReads)

	output:
	path (maxbin_out)

	"""
	mkdir maxbin_out
	run_MaxBin.pl \
	-contig $spades_out/scaffolds.modified.fasta \
	-reads ${cleanedReads[0]} \
	-reads2 ${cleanedReads[2]} \
	-thread 32 \
	-out maxbin_out/bin
	mkdir maxbin_out/maxbin_bins
	cp maxbin_out/*.fasta maxbin_out/maxbin_bins
	"""
}

process minimap {
	publishDir = params.outdir
	conda 'bioconda::minimap2 bioconda::samtools'

	input:
	path (spades_out)
	path (cleanedReads)
	val (id)

	output:
	path (minimap_out)

	"""
	mkdir minimap_out
	minimap2 -ax sr \
	-t 12 \
	$spades_out/scaffolds.modified.fasta \
	${cleanedReads[0]} \
	${cleanedReads[2]} \
	| samtools sort -@12 \
	-O BAM \
	-o minimap_out/$id".bam"
	"""
}

process depth {
	publishDir = params.outdir
	
	input:
	path (minimap_out)
	val (id)

	output:
	path 'depth.txt'

	"""
	jgi_summarize_bam_contig_depths \
	--outputDepth depth.txt \
	$minimap_out/$id".bam"
	"""
}

process metabat {
	publishDir = params.outdir

	input:
	path (spades_out)
	path (depth_out)

	output:
	path (metabat_out)

	"""
	metabat2 \
	-i $spades_out/scaffolds.modified.fasta \
	-a $depth_out \
	-o metabat_out/bin
	mkdir metabat_out/metabat_bins
	mv metabat_out/*.fa metabat_out/metabat_bins
	"""
}

process concoct {
	publishDir = params.outdir

	input:
	path (spades_out)
	path (cleanedReads)

	output:
	path (concoct_out)

	"""
	metaWRAP binning --concoct \
	-a $spades_out/scaffolds.modified.fasta \
	-o concoct_out \
	${cleanedReads[0]} \
	${cleanedReads[2]}
	rm concoct_out/concoct_bins/unbinned.fa 
	"""
}

process table {
        publishDir = params.outdir
        conda 'bioconda::das_tool'

        input:
        path (maxbin_out)
        path (metabat_out)
        path (concoct_out)

        output:
        path 'maxbin_contig2bin.tsv'
        path 'metabat_contig2bin.tsv'
        path 'concoct_contig2bin.tsv'

        """
        Fasta_to_Contig2Bin.sh -e fasta -i $maxbin_out/maxbin_bins > maxbin_contig2bin.tsv
        Fasta_to_Contig2Bin.sh -e fa -i $metabat_out/metabat_bins > metabat_contig2bin.tsv
        Fasta_to_Contig2Bin.sh -e fa -i $concoct_out/concoct_bins > concoct_contig2bin.tsv
        """
}

process dastool {
        publishDir = params.outdir
        conda 'bioconda::das_tool'

        input:
        path (maxbin)
        path (metabat)
        path (concoct)
        path (spades_out)
        val (id)

        output:
        path (dastool_out)

        """
        mkdir dastool_out
        DAS_Tool --write_bin_evals --write_bins --score_threshold=0.1 -t 12 -i $maxbin,$metabat,$concoct --labels maxbin,metabat,concoct -c $spades_out/scaffolds.modified.fasta -o dastool_out/$id
        """
}

process BAT {
	publishDir = params.outdir

	input:
	path (dastool_out)
	val(id)

	output:
	path (BAT_out)

	script:
	database = params.catdb
	taxonomy = params.taxdb
	
	"""
	mkdir BAT_out
	CAT bins \
	-b $dastool_out/$id"_DASTool_bins" \
	-d $database \
	-t $taxonomy \
	--bin_suffix .fa \
	-o BAT_out/$id
	"""
}

process CAT {
	publishDir = params.outdir

	input:
	path (spades_out)
	val (id)	

	output:
	path (CAT_out)

	script:
	database = params.catdb
        taxonomy = params.taxdb

	"""
	mkdir CAT_out
	CAT contigs \
	-c $spades_out/scaffolds.modified.fasta \
	-d $database \
	-t $taxonomy \
	-o CAT_out/$id
	"""
}

process adding_names_BAT {
	publishDir = params.outdir

	input:
	path (BAT_out)
	val (id)

	output:
	path 'bins_taxonomy.txt'

	script:
	taxonomy = params.taxdb
	
	"""
	CAT add_names \
	-i $BAT_out/$id".bin2classification.txt" \
	-o bins_taxonomy.txt \
	-t $taxonomy \
	--only_official
	"""
}

process adding_names_CAT {
	publishDir = params.outdir

	input:
	path (CAT_out)
	val (id)
	
	output:
	path 'contigs_taxonomy.txt'

	script:
	taxonomy = params.taxdb

	"""
	CAT add_names \
	-i $CAT_out/$id".contig2classification.txt" \
	-o contigs_taxonomy.txt \
	-t $taxonomy \
	--only_official
	"""
}

process summarize_CAT {
	publishDir = params.outdir

	input:
	path (spades_out)
	path 'contigs_taxonomy.txt'
	val (id)

	output:
	path 'CAT_summary.txt'
	
	"""
	CAT summarise \
	-c $spades_out/scaffolds.modified.fasta \
	-i contigs_taxonomy.txt \
	-o CAT_summary.txt
	"""
}

process assembly_stat {
	publishDir = params.outdir

	input:
	path (spades_out)

	output:
	path 'assembly.stats.txt'

	"""
	gaas_fasta_statistics.pl \
	-f $spades_out/scaffolds.modified.fasta > assembly.stats.txt
	"""
}

process bins_abundance {
	publishDir = params.outdir
	
	input:
	path (dastool_out)
	val (id)
	path (spades_out)
	path (cleanedReads)

	output:
	path (bins_abundance)
	
	"""
	mkdir bins_abundance
	metaWRAP quant_bins \
	-b $dastool_out/$id"_DASTool_bins" \
	-a $spades_out/scaffolds.modified.fasta \
	-o bins_abundance/$id \
	${cleanedReads[0]} \
	${cleanedReads[2]} \		
	"""
}

process metaxa {
	publishDir = params.outdir
	
	input:
	path (cleanedReads)
	val (id)
	
	output:
	path (metaxa_out)

	"""
	mkdir metaxa_out
	metaxa2 \
	-1 ${cleanedReads[0]} \
	-2 ${cleanedReads[2]} \
	-E 1e-5 \
	-t all \
	-o metaxa_out/$id
	"""
}

process ITSx {
	publishDir = params.outdir
	conda 'bioconda::itsx'

	input:
	path (spades_out)
	val (id)

	output:
	path (ITSx_out)

	"""
	mkdir ITSx_out
	ITSx \
	-i $spades_out/scaffolds.modified.fasta \
	-o ITSx_out/$id \
	--cpu 12
	"""
}

process kaiju {
	publishDir = params.outdir
	conda 'bioconda::kaiju'

	input:
	path (cleanedReads)
	val (id)

	output:
	path (kaiju_out)	
	
	script:
	kaiju_db = params.kaiju_db
	kaiju_nodes = params.kaiju_nodes

	"""
	mkdir kaiju_out
	kaiju \
	-E 1e-5 \
	-t $kaiju_nodes \
	-f $kaiju_db \
	-i ${cleanedReads[0]} \
	-j ${cleanedReads[2]} \
	-o kaiju_out/$id"_kaiju.out"
	"""
}

process kaiju_summary {
	publishDir = params.outdir
	conda 'bioconda::kaiju'

	input:
	val (id)
	path (kaiju_out)
	
	output:
	path (kaiju_out)

	script:
	kaiju_nodes = params.kaiju_nodes
	kaiju_names = params.kaiju_names

	"""
	kaiju2table \
	-t $kaiju_nodes \
	-n $kaiju_names \
	-r genus \
	-o kaiju_out/$id"_summary.tsv" \
	$kaiju_out/$id"_kaiju.out"
	"""
}

process reassemble {
        publishDir = params.outdir

        input:
        path (dastool_out)
	val (id)
        path (cleanedReads)

        output:
        path (reassembly)
	
        """
        metaWRAP reassemble_bins \
        -o reassembly \
        -b $dastool_out/$id"_DASTool_bins" \
        -1 ${cleanedReads[0]} \
        -2 ${cleanedReads[2]} \
        -t 12 \
        --skip-checkm \
        --parallel
        mkdir reassembly/reassembled_bins/strict
        mkdir reassembly/reassembled_bins/permissive
        mv reassembly/reassembled_bins/*.strict.fa reassembly/reassembled_bins/strict
        mv reassembly/reassembled_bins/*.permissive.fa reassembly/reassembled_bins/permissive
        """
}

workflow {
		
	reads_ch = Channel.fromFilePairs(params.seqs)
	
	trimmomatic(reads_ch)

	assembly(trimmomatic.out)
	
	maxbin(assembly.out, trimmomatic.out)

	minimap(assembly.out, trimmomatic.out, reads_ch.map {id, reads -> id})

	depth(minimap.out, reads_ch.map{id, reads -> id})

	metabat(assembly.out, depth.out)

	concoct(assembly.out, trimmomatic.out)

	table(maxbin.out, metabat.out, concoct.out)

    dastool(table.out, assembly.out, reads_ch.map {id, reads -> id})

	reassemble(dastool.out, reads_ch.map {id, reads -> id}, trimmomatic.out)

	BAT(dastool.out, reads_ch.map{id, reads -> id})

	CAT(assembly.out, reads_ch.map{id, reads -> id})
	
	adding_names_BAT(BAT.out, reads_ch.map {id, reads -> id})
	
	adding_names_CAT(CAT.out, reads_ch.map {id, reads -> id})

	summarize_CAT(assembly.out, adding_names_CAT.out, reads_ch.map{id, reads -> id})

	assembly_stat(assembly.out)

	bins_abundance(dastool.out, reads_ch.map {id, reads -> id}, assembly.out, trimmomatic.out)

	metaxa(trimmomatic.out, reads_ch.map {id, reads -> id})

	ITSx(assembly.out, reads_ch.map {id, reads -> id})

	kaiju(trimmomatic.out, reads_ch.map {id, reads -> id})

	kaiju_summary(reads_ch.map {id, reads -> id}, kaiju.out)
}
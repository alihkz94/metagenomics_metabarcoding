#!/usr/bin/env nextflow

params.seq = '/nfs/turbo/lsa-tyleho/mycology/Yongjie/Entorrhiza/New_Reads/PacBio_Data/FSXU0295R.hifi_reads.fastq'
params.outdir = '/scratch/tyleho_root/tyleho0/hkz/Results'
params.catdb = '/nfs/turbo/lsa-tyleho/mycology/MMC/20231120_CAT_nr/db'
params.taxdb = '/nfs/turbo/lsa-tyleho/mycology/MMC/20231120_CAT_nr/tax'
params.kaiju_names = '/scratch/tyleho_root/tyleho0/hkz/kaijudb/names.dmp'
params.kaiju_nodes = '/scratch/tyleho_root/tyleho0/hkz/kaijudb/nodes.dmp'
params.kaiju_db = '/scratch/tyleho_root/tyleho0/hkz/kaijudb/kaiju_db_nr_euk.fmi'
params.emapper_db = '/nfs/turbo/lsa-tyleho/mycology/MMC/data'
params.emapper_db_name = "Bacteria"

process Flye {
	publishDir = params.outdir

	input:
	path (Read)
	
	output:
	path (Flye_out)

	"""
	mkdir Flye_out
	flye \
	--pacbio-hifi $Read \
	--out-dir Flye_out \
	--meta \
	--threads 24	
	"""
}

process minimap {
        publishDir = params.outdir
        conda 'bioconda::minimap2 bioconda::samtools'

        input:
        path (Flye_out)
        path (Read)
	val (id)

        output:
        path (minimap_out)

        """
        mkdir minimap_out
        minimap2 -ax map-pb \
       	-t 24 \
       	$Flye_out/assembly.fasta \
       	$Read \
       	| samtools sort -@24 \
       	-O BAM \
       	-o minimap_out/$id".bam"
        """
}

process SemiBin {
	publishDir = params.outdir

	input:
	path (Flye_out)
	path (minimap_out)
	val (id)
	
	output:
	path (SemiBin_out)

	"""
	mkdir SemiBin_out
	SemiBin2 \
	single_easy_bin \
	--sequencing-type=long_read \
	-i $Flye_out/assembly.fasta \
	-b $minimap_out/$id".bam" \
	-o SemiBin_out
	gunzip SemiBin_out/output_bins/*.gz
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
	path (Flye_out)
	path (depth_out)

	output:
	path (metabat_out)

	"""
	metabat2 \
	-i $Flye_out/assembly.fasta \
	-a $depth_out \
	-o metabat_out/bin
	mkdir metabat_out/metabat_bins
	mv metabat_out/*.fa metabat_out/metabat_bins
	"""
}

process table {
        publishDir = params.outdir
        conda 'bioconda::das_tool'

        input:
        path (SemiBin_out)
        path (metabat_out)

        output:
        path 'SemiBin_contig2bin.tsv'
        path 'metabat_contig2bin.tsv'

        """
        Fasta_to_Contig2Bin.sh -e fa -i $SemiBin_out/output_bins > SemiBin_contig2bin.tsv
        Fasta_to_Contig2Bin.sh -e fa -i $metabat_out/metabat_bins > metabat_contig2bin.tsv
        """
}

process dastool {
        publishDir = params.outdir
        conda 'bioconda::das_tool'

        input:
        path (SemiBin)
        path (metabat)
        path (Flye_out)
        val (id)

        output:
        path (dastool_out)

        """
        mkdir dastool_out
        DAS_Tool --write_bin_evals --write_bins --score_threshold=0.1 -t 12 -i $SemiBin,$metabat --labels SemiBin,metabat -c $Flye_out/assembly.fasta -o dastool_out/$id
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
	path (Flye_out)
	val (id)	

	output:
	path (CAT_out)

	script:
	database = params.catdb
        taxonomy = params.taxdb

	"""
	mkdir CAT_out
	CAT contigs \
	-c $Flye_out/assembly.fasta \
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
	path (Flye_out)
	path 'contigs_taxonomy.txt'
	val (id)

	output:
	path 'CAT_summary.txt'
	
	"""
	CAT summarise \
	-c $Flye_out/assembly.fasta \
	-i contigs_taxonomy.txt \
	-o CAT_summary.txt
	"""
}

process assembly_stat {
	publishDir = params.outdir

	input:
	path (Flye_out)

	output:
	path 'assembly.stats.txt'

	"""
	gaas_fasta_statistics.pl \
	-f $Flye_out/assembly.fasta > assembly.stats.txt
	"""
}

process ITSx {
	publishDir = params.outdir
	conda 'bioconda::itsx'

	input:
	path (Flye_out)
	val (id)

	output:
	path (ITSx_out)

	"""
	mkdir ITSx_out
	ITSx \
	-i $Flye_out/assembly.fasta \
	-o ITSx_out/$id \
	--cpu 12
	"""
}

workflow {
		
	reads_ch = Channel.fromPath(params.seq, checkIfExists: true)

	Flye(reads_ch)

	minimap(Flye.out, reads_ch, reads_ch.map { fn -> fn.simpleName })

	SemiBin(Flye.out, minimap.out, reads_ch.map { fn -> fn.simpleName })

	depth(minimap.out, reads_ch.map { fn -> fn.simpleName })

	metabat(Flye.out, depth.out)

	table(SemiBin.out, metabat.out)

    dastool(table.out, Flye.out, reads_ch.map { fn -> fn.simpleName })

	BAT(dastool.out, reads_ch.map { fn -> fn.simpleName })

	CAT(Flye.out, reads_ch.map { fn -> fn.simpleName })
	
	adding_names_BAT(BAT.out, reads_ch.map { fn -> fn.simpleName })
	
	adding_names_CAT(CAT.out, reads_ch.map { fn -> fn.simpleName })

	summarize_CAT(Flye.out, adding_names_CAT.out, reads_ch.map { fn -> fn.simpleName })

	assembly_stat(Flye.out)

	ITSx(Flye.out, reads_ch.map { fn -> fn.simpleName })

}
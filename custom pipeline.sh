#!/bin/bash
set -e # Exit on error

# Initial Cleanup: Remove existing directories or files from past analyses
DIRECTORIES=("cut_primers_out" "quality_filtering_out" "fasta_files" "uchime_denovo_out" \
"uchime_ref_out" "ITSx_out" "derep_out" "clustering_out" "blast_out")

echo "Cleaning up old directories if they exist..."
for dir in "${DIRECTORIES[@]}"; do
    if [ -d "$dir" ]; then
        rm -r "$dir"
        echo "Removed old directory: $dir"
    fi
done

# 1. Check and Decompress Input Files (if compressed)
echo "Decompressing files if necessary..."
if ls *.gz 1> /dev/null 2>&1; then
    for file in *.gz; do
        gunzip $file
    done
else
    echo "No compressed files found. Moving to the next step."
fi


INPUT_FILES=(*.fastq)

# 2. Trim Primers with Cutadapt
echo "Trimming Primers..."
mkdir -p cut_primers_out
for file in "${INPUT_FILES[@]}"; do
    READ_IN=$(seqkit stats $file | awk 'NR==2 {print $2}')
    cutadapt -a TGTACACACCGCCCGTCG -g "TCCTSCGCTTATTGATATGC;min_overlap=20" -e 0.1 --discard-untrimmed -j 8 \
    --action=none -o cut_primers_out/$file $file
    READ_OUT=$(seqkit stats cut_primers_out/$file | awk 'NR==2 {print $2}')
done

# 3. Quality Filtering with Vsearch
echo "Quality Filtering..."
mkdir -p quality_filtering_out
for file in cut_primers_out/*.fastq; do
    READ_IN=$(seqkit stats $file | awk 'NR==2 {print $2}')
    vsearch --fastq_filter $file --fastq_maxee 1 --fastq_maxns 0 --fastq_minlen 50 --threads 8 --fastq_qmax 93 \
    --fastq_qmin 0 --fastqout quality_filtering_out/${file##*/}
    READ_OUT=$(seqkit stats quality_filtering_out/${file##*/} | awk 'NR==2 {print $2}')
done

# 4. Convert FASTQ to FASTA
echo "Converting FASTQ to FASTA..."
mkdir -p fasta_files
for file in quality_filtering_out/*.fastq; do
    if [ -s "$file" ]; then
        seqkit fq2fa $file > fasta_files/$(basename $file .fastq).fasta
    else
        echo "Warning: $file is empty. Skipping..."
    fi
done

# 5. Chimera Filtering with UCHIME Denovo
echo "Chimera Filtering with UCHIME Denovo..."
mkdir -p uchime_denovo_out
for file in fasta_files/*.fasta; do
    if [ -s "$file" ]; then  # Check if file is not empty
        vsearch --uchime_denovo $file --mindiv 0.5 --dn 1.6 --threads 8 --uchimeout uchime_denovo_out/${file##*/}.uchime --nonchimeras uchime_denovo_out/${file##*/}.fasta
    else
        echo "Warning: $file is empty. Skipping..."
    fi
done


# 6. Chimera Filtering with UCHIME Ref
echo "Chimera Filtering with UCHIME Ref..."
mkdir -p uchime_ref_out
for file in uchime_denovo_out/*.fasta; do
    if [ -s "$file" ]; then  # Check if file is not empty
        vsearch --uchime_ref $file --mindiv 0.5 --dn 1.6 --threads 8 --db ref.fasta --nonchimeras uchime_ref_out/${file##*/}
    else
        echo "Warning: $file is empty. Skipping..."
    fi
done


# 7. ITS Extraction with ITSx
echo "ITS Extraction..."
mkdir -p ITSx_out
for file in uchime_ref_out/*.fasta; do
    ITSx -i $file --cpu 128 -nhmmer T -E 1e-3 --complement T --only_full T -domains 2 -score 0 --partial 50 \
    --not_found F --graphical F --preserve T --summary F -o ITSx_out/${file##*/}
done

# 8. Dereplication with Vsearch
echo "Dereplication..."
mkdir -p derep_out
for file in ITSx_out/*.full.fasta; do
    vsearch --derep_fulllength $file --sizein --sizeout --threads 8 -o derep_out/${file##*/}
done

# 9. Clustering with Vsearch
echo "Clustering..."
mkdir -p clustering_out
for file in derep_out/*.fasta; do
    vsearch --cluster_size $file --fasta_width 0 --id 0.97 --threads 8 --strands both --remove_singeltons true \
    --sizein --qmask dust --maxaccepts 1 --uc clustering_out/OTUs.uc --centroids clustering_out/OTUs.fasta
done

# 10. Parallel Blast Analysis
echo "Running Blast Analysis..."
mkdir -p blast_out
split -l $(($(wc -l < clustering_out/OTUs.fasta)/10)) clustering_out/OTUs.fasta chunk_
for file in chunk_*; do
    blastn -query $file \
    -db silva_138 \
    -word_size 7 \
    -num_threads 128 \
    -outfmt "6 delim=+ qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" \
    -evalue 0.001 \
    -strand both \
    -max_target_seqs 1 \
    -max_hsps 2 \
    -dust no \
    -soft_masking true \
    -penalty -1 \
    -reward 1 \
    -gapopen 1 \
    -gapextend 2 \
    -out blast_out/${file}_blast_results.txt &
done
wait
cat blast_out/*_blast_results.txt > blast_out/final_blast_results.txt

echo "Pipeline Completed!"

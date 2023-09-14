#!/bin/bash
#SBATCH --job-name="parallel_blast"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --partition amd
#SBATCH --time=90:00:00

DIR="/gpfs/space/home/alihakim/blast"
FOLDERS=("2017" "2019_2020" "2020_2021" "2021_2022")

for folder in "${FOLDERS[@]}"; do
    cd $DIR/$folder

    TOTAL_SEQS=$(grep -c '^>' OTUs.fasta)

    if (( TOTAL_SEQS <= 500 )); then
        NUM_CHUNKS=5
    else
        NUM_CHUNKS=10
    fi

    SEQ_PER_CHUNK=$((TOTAL_SEQS / NUM_CHUNKS))

    # Corrected AWK command for splitting the FASTA file
    awk -v prefix="OTUs_chunk_" -v chunk_size="$SEQ_PER_CHUNK" '/^>/{n++;if(n%chunk_size==1) {m++; close(f); f=prefix m ".fasta"}} {print > f}' OTUs.fasta

    for file in OTUs_chunk_*.fasta; do
        (
            blastn -query $file \
            -db /gpfs/space/home/alihakim/analysis/databse/UNITE \
            -word_size 7 \
            -num_threads 12 \
            -outfmt "6  delim=; qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" \
            -evalue 0.001 \
            -strand both \
            -max_target_seqs 10 \
            -max_hsps 2 \
            -dust no \
            -soft_masking true \
            -penalty -1 \
            -reward 1 \
            -gapopen 1 \
            -gapextend 2 \
            -out ${file%.fasta}_blast_results.txt

            rm $file
        ) &
    done
    wait

    cat OTUs_chunk_*_blast_results.txt > blast_top10hit.txt
    awk -F '\t' '!seen[$1]++' _blast_top10hit.txt > blast_1st_tophit.txt
    rm OTUs_chunk_*_blast_results.txt
done

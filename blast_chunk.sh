#!/bin/bash
#SBATCH --job-name="parallel_blast"
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --partition amd
#SBATCH --time=70:00:00  # Adjust based on your estimation of the total runtime for all folders

# Directory containing the folders
DIR="/gpfs/space/home/alihakim/blast"

# List of folders to process
FOLDERS=("2017" "2019_2020" "2020_2021" "2021_2022")

# Loop through each folder
for folder in "${FOLDERS[@]}"; do
    cd $DIR/$folder

    # Delete any remaining chunks from previous analyses
    rm -f OTUs_chunk_*.fasta

    # Delete any BLAST results from previous runs on chunks
    rm -f OTUs_chunk_*_blast_results.txt


# Loop through each folder
for folder in "${FOLDERS[@]}"; do
    cd $DIR/$folder

    # Count the total number of sequences in the FASTA file
    TOTAL_SEQS=$(grep -c '^>' OTUs.fasta)

    # Determine the number of chunks based on the total sequences (between 5 and 10)
    if (( TOTAL_SEQS <= 500 )); then
        NUM_CHUNKS=5
    else
        NUM_CHUNKS=10
    fi

    SEQ_PER_CHUNK=$((TOTAL_SEQS / NUM_CHUNKS))

    # Split the fasta file into roughly equal chunks based on SEQ_PER_CHUNK
    awk -v prefix="OTUs_chunk_" -v chunk_size="$SEQ_PER_CHUNK" '/^>/{n++;if(n%chunk_size==1 && m<NUM_CHUNKS-1) {m++; close(f); f=prefix m ".fasta"}} {print > f}' OTUs.fasta

    # Run BLAST analysis for each chunk simultaneously
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

            # Delete the chunk after analysis
            rm $file
        ) &
    done
    wait # Wait for all parallel jobs to complete

    # After all jobs are complete, merge the results
    cat OTUs_chunk_*_blast_results.txt > _blast_top10hit.txt

    # Extract the first top hit for each OTU
    # Extract the first top hit for each OTU
    awk -F '\t' '!seen[$1]++' _blast_top10hit.txt > blast_1st_tophit.txt


    # Delete the results chunks
    rm OTUs_chunk_*_blast_results.txt
done

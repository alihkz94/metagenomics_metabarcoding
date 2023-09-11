#!/bin/bash
#SBATCH --job-name="sequential_blast"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --partition amd
#SBATCH --time=70:00:00  # Adjust based on your estimation of the total runtime for all folders

# Directory containing the folders
DIR="/gpfs/space/home/alihakim/blast"

# Define the number of chunks
NUM_CHUNKS=10

# List of folders to process
FOLDERS=("2014" "2015" "2016" "2017" "2019_2020" "2020_2021" "2021_2022")

# Loop through each folder
for folder in "${FOLDERS[@]}"; do
    cd $DIR/$folder

    # Split the fasta file into equal chunks
    awk -v prefix="OTUs_chunk_" '/^>/{n++;if(n%'$NUM_CHUNKS'==1) close(f); f=prefix n ".fasta"} {print > f}' OTUs.fasta

    # Loop through each chunk and run BLAST
    for file in OTUs_chunk_*.fasta; do
        blastn -query $file \
        -db silva_138 \
        -word_size 7 \
        -num_threads 128 \
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
    done

    # After all jobs are complete, merge the results
    cat OTUs_chunk_*_blast_results.txt > _blast_top10hit.txt

    # Extract the first top hit for each OTU
    awk -F 'delim=;' '!seen[$1]++' _blast_top10hit.txt > blast_1st_tophit.txt

    # Delete the chunks after analysis
    rm OTUs_chunk_*.fasta
done

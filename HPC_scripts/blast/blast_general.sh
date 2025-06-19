#!/bin/bash
#SBATCH --job-name="ITS"
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --mem=6G
#SBATCH --partition=amd
#SBATCH --time=90:00:00
#SBATCH --array=0-9

DIR="/gpfs/helios/home/alihakim/leho_euk/ITS"
HEADER="qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident"

# Get the list of fasta files
FASTA_FILES=($(ls ${DIR}/*.fasta))
NUM_FILES=${#FASTA_FILES[@]}

# Calculate which file and chunk we're processing
FILE_IDX=$((SLURM_ARRAY_TASK_ID / 10))
CHUNK_IDX=$((SLURM_ARRAY_TASK_ID % 10))

if [ $FILE_IDX -ge $NUM_FILES ]; then
    echo "No more files to process"
    exit 0
fi

FASTA_FILE=${FASTA_FILES[$FILE_IDX]}
BASENAME=$(basename "$FASTA_FILE" .fasta)

# Create directory for output
mkdir -p "$DIR/$BASENAME"

# Split the file into chunks if not already done
if [ ! -f "$DIR/$BASENAME/${BASENAME}_chunk_${CHUNK_IDX}.fasta" ]; then
    TOTAL_SEQS=$(grep -c '^>' "$FASTA_FILE")
    # Use a ceiling approach if you want strictly 10 chunks:
    # SEQ_PER_CHUNK=$(( (TOTAL_SEQS + 9) / 10 ))
    SEQ_PER_CHUNK=$((TOTAL_SEQS / 10))
    awk -v prefix="$DIR/$BASENAME/${BASENAME}_chunk_" -v chunk_size="$SEQ_PER_CHUNK" \
        '/^>/{n++; if(n%chunk_size==1){m++; close(f); f=prefix m-1 ".fasta"}} {print > f}' "$FASTA_FILE"
fi

# Process the assigned chunk
CHUNK_FILE="$DIR/$BASENAME/${BASENAME}_chunk_${CHUNK_IDX}.fasta"
OUTPUT_FILE="${CHUNK_FILE%.fasta}_blast_results.txt"

# Run BLAST
blastn -query "$CHUNK_FILE" \
    -db /gpfs/helios/home/alihakim/leho_euk/database/ITS/EUK \
    -word_size 7 \
    -task blastn \
    -num_threads 32 \
    -outfmt "6 delim=+ $HEADER" \
    -evalue 0.001 \
    -strand both \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -out "$OUTPUT_FILE"

# Add header to the output file
sed -i '1i'"$HEADER" "$OUTPUT_FILE"

# Clean up the chunk file
rm "$CHUNK_FILE"

# Create a completion flag
touch "$DIR/$BASENAME/.chunk_${CHUNK_IDX}_complete"

# Check if all chunks are complete and combine results if they are
NUM_COMPLETE=$(ls "$DIR/$BASENAME/.chunk_"*"_complete" 2>/dev/null | wc -l)
if [ "$NUM_COMPLETE" -eq 10 ]; then
    # Concatenate all chunk results
    cat "$DIR/$BASENAME/${BASENAME}_chunk_"*_blast_results.txt > "$DIR/$BASENAME/combined_blast_top10hit.txt"

    # Convert the *first line* from space-separated to plus-separated
    sed -i '1s/ /+/g' "$DIR/$BASENAME/combined_blast_top10hit.txt"

    # Make a "best hit only" file by selecting unique qseqid in first column
    awk -F '+' '!seen[$1]++' "$DIR/$BASENAME/combined_blast_top10hit.txt" > "$DIR/$BASENAME/${BASENAME}.txt"

    # Convert the *first line* of best hits file similarly
    sed -i '1s/ /+/g' "$DIR/$BASENAME/${BASENAME}.txt"

    # Clean up
    rm "$DIR/$BASENAME/${BASENAME}_chunk_"*_blast_results.txt
    rm "$DIR/$BASENAME/.chunk_"*"_complete"
fi

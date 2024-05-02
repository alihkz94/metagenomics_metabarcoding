#!/bin/bash

# Specify the directory containing the fasta files
DIR="/home/ali/Pictures/blast"
CPUS=5 # Number of CPUs to use for BLAST
BLAST_DB="/home/ali/Documents/simulated_data/database/unite" # Path to the BLAST database

# Header for BLAST output files
HEADER="qseqid+stitle+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch"

# Identify the large file
LARGE_FILE=$(find "$DIR" -name "*.fasta" -size +316M -print -quit)

# Process each FASTA file in the specified directory, excluding the large file
for fasta_file in "$DIR"/*.fasta; do
    if [ "$fasta_file" == "$LARGE_FILE" ]; then
        continue # Skip the large file for now
    fi

    # Extract the basename of the FASTA file
    BASENAME=$(basename "$fasta_file" .fasta)
    # Create a directory specific for the FASTA file's results
    mkdir -p "$DIR/$BASENAME"

    # Count the total number of sequences in the FASTA file
    TOTAL_SEQS=$(grep -c '^>' "$fasta_file")
    # Determine the number of chunks based on the sequence count
    if (( TOTAL_SEQS <= 500 )); then
        NUM_CHUNKS=5
    else
        NUM_CHUNKS=10
    fi

    # Calculate the number of sequences per chunk, avoiding division by zero
    SEQ_PER_CHUNK=$(( (TOTAL_SEQS + NUM_CHUNKS - 1) / NUM_CHUNKS ))
    if (( SEQ_PER_CHUNK == 0 )); then
        echo "Warning: No sequences to process in $fasta_file"
        continue # Skip this file if there are no sequences to divide into chunks
    fi

    # Split the FASTA file into smaller chunks
    awk -v prefix="$DIR/$BASENAME/${BASENAME}_chunk_" -v chunk_size="$SEQ_PER_CHUNK" \
        '/^>/{n++; if(n % chunk_size == 1) { m++; close(f); f=prefix m ".fasta" }} { print > f }' "$fasta_file"

    # Process each chunk with BLAST, only if the chunk file is not empty
    for file in "$DIR/$BASENAME/${BASENAME}_chunk_"*.fasta; do
        if [ -s "$file" ]; then # Check if file is not empty
            (
                # Run BLASTn and handle outputs
                blastn -query "$file" \
                       -db "$BLAST_DB" \
                       -word_size 7 \
                       -task blastn \
                       -num_threads "$CPUS" \
                       -outfmt "6 delim=+ $HEADER" \
                       -evalue 0.001 \
                       -strand both \
                       -max_target_seqs 10 \
                       -max_hsps 1 \
                       -out "${file%.fasta}_blast_results.txt"

                # Add a header to the output file
                sed -i '1i'"$HEADER" "${file%.fasta}_blast_results.txt"
                # Clean up by removing the chunk file after processing
                rm "$file"
            ) &
        else
            echo "Skipping empty file: $file"
        fi
    done

    # Wait for all background processes to finish
    wait

    # Combine results from all chunks, check if files exist before attempting to concatenate
    if compgen -G "$DIR/$BASENAME/${BASENAME}_chunk_*_blast_results.txt" > /dev/null; then
        cat "$DIR/$BASENAME/${BASENAME}_chunk_"*_blast_results.txt > "$DIR/$BASENAME/combined_blast_top10hit.txt"
        # Deduplicate and finalize the results file
        awk -F '+' '!seen[$1]++' "$DIR/$BASENAME/combined_blast_top10hit.txt" > "$DIR/$BASENAME/${BASENAME}.txt"
        # Remove intermediate result files
        rm "$DIR/$BASENAME/${BASENAME}_chunk_"*_blast_results.txt
    else
        echo "No BLAST results to combine for $BASENAME"
    fi
done

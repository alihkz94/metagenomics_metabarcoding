#!/bin/bash

# Directory containing original FASTA files
DIR="/home/ali/Documents/simulated_data/analysis/final results/rechime_test/uchime/chimeras"

# New directory to store reformatted FASTA files
NEW_DIR="$DIR/new_fasta"

# Create the new directory if it does not exist
mkdir -p "$NEW_DIR"

# Process each FASTA file in the directory
for fasta_file in "$DIR"/*.fasta; do
    # Extract the basename of the fasta file
    BASENAME=$(basename "$fasta_file")
    NEW_FILE_PATH="$NEW_DIR/$BASENAME"

    # Use seqkit to reformat the fasta file and save it to the new directory
    seqkit seq --line-width 0 -w 0 "$fasta_file" > "$NEW_FILE_PATH"
done

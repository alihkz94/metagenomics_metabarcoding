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

    # Initialize a variable to hold the full sequence
    full_sequence=""

    # Read through the file line by line
    while IFS= read -r line; do
        # Check if the line contains a sequence description
        if [[ "$line" == ">"* ]]; then
            # If there is a stored sequence, write it to the new file
            if [[ -n "$full_sequence" ]]; then
                echo "$full_sequence" >> "$NEW_FILE_PATH"
                full_sequence=""  # Reset sequence
            fi
            # Write the sequence description to the new file
            echo "$line" >> "$NEW_FILE_PATH"
        else
            # Remove all newline characters and append to the sequence
            full_sequence+="$line"
        fi
    done < "$fasta_file"

    # Write the last sequence to the new file if not empty
    if [[ -n "$full_sequence" ]]; then
        echo "$full_sequence" >> "$NEW_FILE_PATH"
    fi
done

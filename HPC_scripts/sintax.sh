#!/bin/bash

DIR="/home/ali/Pictures/sintax_uchime_best"  # Path to your FASTA files
DB="/home/ali/Documents/simulated_data/sintax_database/UNITE.udb"  # Path to your UDB database

# Process each FASTA file in the directory
for fasta_file in "$DIR"/*.fasta; do
    # Extract the basename of the fasta file
    BASENAME=$(basename "$fasta_file" .fasta)

    # Prepare output directory and files
    mkdir -p "$DIR/$BASENAME"
    OUTPUT="$DIR/$BASENAME/${BASENAME}_sintax_results.txt"
    ALNOUT="$DIR/$BASENAME/${BASENAME}_alignment.aln"

    # Run Vsearch with Sintax
    vsearch --sintax "$fasta_file" \
            --db "$DB" \
            --tabbedout "$OUTPUT" \
            --sintax_cutoff 0.8 \
            --strand both \
            --wordlength 7 \
            --threads 8

    # Check if output exists and echo status
    if [[ -f "$OUTPUT" ]]; then
        echo "Sintax classification completed for $BASENAME"
    else
        echo "Sintax classification failed for $BASENAME"
    fi
done

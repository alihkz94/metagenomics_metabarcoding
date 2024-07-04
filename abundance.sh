#!/bin/bash

# Output file
output_file="abundances.txt"

# Clear the output file
echo "" > $output_file

# Directory containing the FASTA files
dir="/home/ali/Documents/simulated_data/analysis/final_results/original/DADA2"

# Loop over all FASTA files in the directory
for fasta in $dir/*.fasta; do
    # Calculate the total abundance for the current file
    total_abundance=$(grep -o 'size=[0-9]*' "$fasta" | awk -F= '{sum += $2} END {print sum}')
    
    # Get the base name of the file
    base_name=$(basename $fasta)
    
    # Print the filename and total abundance to the output file
    echo "$base_name: $total_abundance" >> $output_file
done

# Print the contents of the output file
cat $output_file

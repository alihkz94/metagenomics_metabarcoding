#!/bin/bash
#SBATCH --job-name="2018_2019"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition amd
#SBATCH --time=180:00:00


# Create the main output directory if it doesn't exist
mkdir -p itsx_out

for fasta in *.fasta; do
    # Create a subdirectory for the current file
    subdir="itsx_out/${fasta%.*}"
    mkdir -p "$subdir"

    # Run ITSx with the specified settings and save the output in the subdirectory
    ITSx \
    --cpu 128 \
    -i "$fasta" \
    -o "$subdir/${fasta%.*}" \
    -nhmmer T \
    -E 1e-2 \
    --complement T \
    --only_full T \
    -domains 2 \
    -score 0 \
    --partial 50 \
    --save_regions ITS2 \
    --not_found F \
    --graphical F \
    --preserve T \
    --summary F 
done
#!/bin/bash
#SBATCH --job-name="BlasCh_nonchimeric"
#SBATCH --cpus-per-task=64
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition amd
#SBATCH --time=90:00:00

# HPC-Optimized BlasCh Execution
# With 64 CPUs and 8 parallel BLAST jobs, each job uses 8 CPUs
# This allows processing 8 files simultaneously for maximum speed

python Blasch_nonchimeric.py --input_nonchimeric_dir ./nonchimeric \
--reference_db ../../database/bold.fasta \
--threads 64 \
--max_parallel_blast 8 \
--self_fasta_dir ./original \
--output_dir ./nonchimeric_results
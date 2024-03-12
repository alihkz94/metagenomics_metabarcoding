import os
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict

# Function to check if a directory exists
def check_dir(directory):
    if not os.path.isdir(directory):
        print(f"Error: Directory {directory} does not exist.")
        sys.exit(1)

# Function to validate FASTA format
def validate_fasta(file_path):
    with open(file_path, 'r') as file:
        if not any(line.startswith('>') for line in file):
            print(f"Error: File {file_path} does not appear to be in FASTA format.")
            sys.exit(1)

# Parse command-line options
parser = argparse.ArgumentParser(description='Chimera Rescue Script')
parser.add_argument('-c', '--chimeric_dir', required=True, help='Directory of chimeric reads')
parser.add_argument('-n', '--non_chimeric_dir', required=True, help='Directory of non-chimeric reads')
args = parser.parse_args()

# Check if directories exist
check_dir(args.chimeric_dir)
check_dir(args.non_chimeric_dir)

# Automatically set the new_non_chimeric_dir
new_non_chimeric_dir = os.path.join(args.non_chimeric_dir, "Rechime_non_chimeric")

# Create a subfolder for rescued sequences
seq_rescued_dir = os.path.join(new_non_chimeric_dir, "seq_rescued")

# Path for the pooled rescued sequences file
rechimed_pool_file = os.path.join(seq_rescued_dir, "rechimed.pool.fasta")

# Remove the existing new_non_chimeric_dir if it exists and create a new one, including the seq_rescued subdirectory
if os.path.isdir(new_non_chimeric_dir):
    os.system(f"rm -rf {new_non_chimeric_dir}")
os.makedirs(seq_rescued_dir)

# Report file path
report_file = os.path.join(new_non_chimeric_dir, "report.txt")

# Minimum sequence occurrence
min_occurrence = 2

# Initialize report file
with open(report_file, 'w') as report:
    report.write("Filename,Before,After,Rescued\n")
total_rescued = 0

# Main analysis logic
for chimeric_filename in os.listdir(args.chimeric_dir):
    if chimeric_filename.endswith(".chimeras.fasta"):
        chimera_file = os.path.join(args.chimeric_dir, chimeric_filename)

        # Skip validation for empty chimeric files
        if os.stat(chimera_file).st_size == 0:
            print(f"Notice: Skipping empty file {chimera_file}.")
            continue

        # Validate FASTA format
        validate_fasta(chimera_file)

        basename = chimeric_filename.replace(".chimeras.fasta", "")
        non_chimeric_file = os.path.join(args.non_chimeric_dir, f"{basename}.fasta")
        new_non_chimeric_file = os.path.join(new_non_chimeric_dir, f"{basename}.fasta")

        # Copy the original non-chimeric file to the new directory
        os.system(f"cp {non_chimeric_file} {new_non_chimeric_file}")

        # Count sequences in the original non-chimeric file
        count_before = sum(1 for _ in SeqIO.parse(non_chimeric_file, "fasta"))

        # Initialize counter for rescued sequences
        rescued = 0

        # Process chimeric sequences
        seq_header_map = defaultdict(str)
        for record in SeqIO.parse(chimera_file, "fasta"):
            seq_header_map[str(record.seq)] += f">{record.id}\n"

        with open(new_non_chimeric_file, 'a') as new_file, open(rechimed_pool_file, 'a') as pool_file:
            for sequence, headers in seq_header_map.items():
                count = headers.count('>')
                if count >= min_occurrence:
                    header = headers.split('\n', 1)[0]
                    new_file.write(f"{header}\n{sequence}\n")
                    pool_file.write(f"{header}\n{sequence}\n")
                    rescued += 1

        # Count sequences in the enhanced non-chimeric file and report
        count_after = sum(1 for _ in SeqIO.parse(new_non_chimeric_file, "fasta"))
        with open(report_file, 'a') as report:
            report.write(f"{basename},{count_before},{count_after},{rescued}\n")
        total_rescued += rescued

# Final report
with open(report_file, 'a') as report:
    report.write(f"Total Rescued Sequences: {total_rescued}\n")
print("Processing complete. Chimeric sequences processed into non-chimeric files.")

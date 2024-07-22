#!/bin/bash

# Define input and output files
OTU_TABLE="OTU_table.txt"
FASTA_FILE="OTUs.fasta"
SINGLETONS_FILE="singletons.txt"
NON_SINGLETON_HEADERS="non_singleton_headers.txt"
CLEANED_FASTA_FILE="OTU_without_singeltons.fasta"
CLEANED_OTU_TABLE="OTU_without_singeltons_table.txt"
ALL_HEADERS_FILE="all_headers.txt"

# Step 1: Identify Singleton OTUs in the OTU Table
awk 'BEGIN {FS=OFS="\t"} 
NR > 1 {
  count=0
  for (i=2; i<=NF; i++) {
    if ($i == 1) count++
    else if ($i > 0) {count=0; break}
  }
  if (count == 1) print $1
}' $OTU_TABLE > $SINGLETONS_FILE

# Step 2: Filter Singleton OTUs from the OTU Table
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1]; next} !($1 in a)' $SINGLETONS_FILE $OTU_TABLE > $CLEANED_OTU_TABLE

# Step 3: Use embedded Python script to clean the FASTA file
python3 - <<EOF
import sys

def remove_singletons_from_fasta(fasta_file, singletons_file, output_file):
    # Load singleton headers
    with open(singletons_file, 'r') as f:
        singletons = set(line.strip() for line in f)

    # Process FASTA file
    with open(fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
        write_sequence = True
        for line in infile:
            if line.startswith('>'):
                header = line[1:].strip()
                write_sequence = header not in singletons
            if write_sequence:
                outfile.write(line)

if __name__ == "__main__":
    fasta_file = "$FASTA_FILE"
    singletons_file = "$SINGLETONS_FILE"
    output_file = "$CLEANED_FASTA_FILE"
    remove_singletons_from_fasta(fasta_file, singletons_file, output_file)
EOF

echo "Singleton sequences removed successfully."
echo "Cleaned OTU table: $CLEANED_OTU_TABLE"
echo "Cleaned FASTA file: $CLEANED_FASTA_FILE"

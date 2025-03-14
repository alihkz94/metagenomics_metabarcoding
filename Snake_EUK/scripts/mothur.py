#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import csv
import os
from utils import get_ascii_fasta_stream

def clean_taxonomy(taxonomy_str):
    # Split taxonomy string into parts, trim whitespace, and remove empty fields.
    parts = [part.strip() for part in taxonomy_str.split(';') if part.strip()]
    # Remove trailing dot placeholders.
    while parts and parts[-1] == '.':
        parts.pop()
    return ';'.join(parts)

def generate_mothur_files(input_fasta, output_fasta, output_tax):
    stream = get_ascii_fasta_stream(input_fasta)
    with open(output_fasta, 'w') as fasta_file, open(output_tax, 'w', newline='') as tax_file:
        tax_writer = csv.writer(tax_file, delimiter='\t')
        for record in SeqIO.parse(stream, 'fasta'):
            # Clean header by removing quotes and enforcing ASCII.
            header = record.description.encode('ascii', errors='replace').decode('ascii').replace('"', '')
            parts = header.split(';')
            feature_id = parts[0].strip()
            taxonomy_raw = ';'.join(parts[1:])
            taxonomy = clean_taxonomy(taxonomy_raw)
            # Write taxonomy info to the TAX file.
            tax_writer.writerow([feature_id, taxonomy])
            # Write FASTA file with the feature ID as header.
            seq_str = str(record.seq).encode('ascii', errors='replace').decode('ascii')
            fasta_file.write(f">{feature_id}\n{seq_str}\n")
    stream.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Mothur formatted files with cleaned taxonomy fields.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file for Mothur")
    parser.add_argument("--output_tax", required=True, help="Output TAX file for Mothur")
    args = parser.parse_args()
    generate_mothur_files(args.input, args.output_fasta, args.output_tax)

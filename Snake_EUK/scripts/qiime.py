#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import csv
import os
from utils import get_ascii_fasta_stream

def generate_qiime2_files(input_fasta, output_fasta, output_tsv):
    stream = get_ascii_fasta_stream(input_fasta)
    with open(output_fasta, 'w') as fasta_file, open(output_tsv, 'w', newline='') as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        tsv_writer.writerow(["Feature ID", "Taxon"])
        for record in SeqIO.parse(stream, 'fasta'):
            header = record.description.encode('ascii', errors='replace').decode('ascii').replace('"', '')
            parts = header.split(';')
            feature_id = parts[0].strip()
            taxonomy = ';'.join(part.strip() for part in parts[1:] if part.strip())
            tsv_writer.writerow([feature_id, taxonomy])
            seq_str = str(record.seq).encode('ascii', errors='replace').decode('ascii')
            fasta_file.write(f">{feature_id}\n{seq_str}\n")
    stream.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate QIIME2 formatted files.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file for QIIME2")
    parser.add_argument("--output_tsv", required=True, help="Output TSV file for QIIME2")
    args = parser.parse_args()
    generate_qiime2_files(args.input, args.output_fasta, args.output_tsv)

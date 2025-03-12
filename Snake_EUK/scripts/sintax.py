#!/usr/bin/env python3
import argparse
import os
from utils import get_ascii_fasta_stream

def convert_blast_to_sintax(input_fasta, output_fasta):
    with get_ascii_fasta_stream(input_fasta) as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Clean header and split taxonomy.
                line_clean = line[1:].strip().replace('"', '')
                parts = line_clean.split(';')
                seq_id = parts[0].strip()
                taxonomy = parts[1:]
                sintax_taxonomy = []
                rank_dict = {
                    'k__': 'd:',
                    'p__': 'p:',
                    'c__': 'c:',
                    'o__': 'o:',
                    'f__': 'f:',
                    'g__': 'g:',
                    's__': 's:'
                }
                for taxon in taxonomy:
                    for key, value in rank_dict.items():
                        if taxon.strip().startswith(key):
                            taxon_name = taxon.split('__')[1].strip()
                            if taxon_name.lower() != 'unclassified':
                                sintax_taxonomy.append(f"{value}{taxon_name}")
                            break
                sintax_header = f">{seq_id};tax={','.join(sintax_taxonomy)};\n"
                outfile.write(sintax_header)
            else:
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FASTA file to SINTAX format.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file for SINTAX")
    args = parser.parse_args()
    convert_blast_to_sintax(args.input, args.output)

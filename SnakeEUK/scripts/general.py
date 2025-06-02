#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import os
import subprocess
import tempfile
from utils import get_ascii_fasta_stream

def transform_header(header):
    # Remove quotes and enforce ASCII.
    header = header.encode('ascii', errors='replace').decode('ascii').replace('"', '')
    parts = header.split(';')
    transformed_parts = [parts[0].strip()]
    taxonomic_ranks = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    for rank, part in zip(taxonomic_ranks, parts[1:]):
        part = part.strip()
        if part == '.' or part.lower() == 'unused':
            transformed_parts.append(rank + 'unclassified')
        else:
            transformed_parts.append(rank + part)
    # Fill in missing ranks if necessary.
    for i in range(len(parts) - 1, len(taxonomic_ranks)):
        transformed_parts.append(taxonomic_ranks[i] + 'unclassified')
    return ';'.join(transformed_parts)

def process_fasta(input_file, output_file):
    stream = get_ascii_fasta_stream(input_file)
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(stream, 'fasta'):
            new_header = transform_header(record.description)
            # Ensure the sequence is strictly ASCII.
            seq_str = str(record.seq).encode('ascii', errors='replace').decode('ascii')
            out_f.write(f">{new_header}\n{seq_str}\n")
    stream.close()
    # Create a unique temporary file for seqkit processing.
    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
        tmp_filepath = tmp_file.name
    subprocess.run(["seqkit", "seq", "-u", output_file, "-o", tmp_filepath], check=True)
    subprocess.run(["seqkit", "seq", "-w", "0", tmp_filepath, "-o", output_file], check=True)
    os.remove(tmp_filepath)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTA file for general taxonomy assignment.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    args = parser.parse_args()
    process_fasta(args.input, args.output)

#!/usr/bin/env python3
import argparse
from utils import get_ascii_fasta_stream

def convert_to_dada2_format(input_file, output_file):
    stream = get_ascii_fasta_stream(input_file)
    with open(output_file, 'w') as outfile:
        for line in stream:
            if line.startswith('>'):
                # Remove the '>' and any quotes.
                header_line = line[1:].strip().replace('"', '')
                # Split header fields by semicolon.
                fields = header_line.split(';')
                # Discard the first field (accession number).
                taxonomy_fields = fields[1:]
                valid_taxa = []
                # Iterate over taxonomy fields and stop at the first placeholder.
                for field in taxonomy_fields:
                    field = field.strip()
                    if field == '.' or field == '':
                        break
                    valid_taxa.append(field)
                # Construct the new header.
                new_header = ">" + ";".join(valid_taxa) + ";"
                outfile.write(new_header + "\n")
            else:
                outfile.write(line)
    stream.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert FASTA file to DADA2 format by removing accession number and truncating at the first placeholder."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file for DADA2")
    args = parser.parse_args()
    convert_to_dada2_format(args.input, args.output)

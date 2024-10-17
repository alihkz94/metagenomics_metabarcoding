import os

def concatenate_blast_outputs(input_dir, output_file):
    """
    Concatenate multiple BLAST output files into a single file.
    
    Parameters:
        input_dir (str): Directory containing the input BLAST output files.
        output_file (str): Path to the output file.
    """
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.txt')]
    header_written = False
    
    with open(output_file, 'w') as outfile:
        for file in files:
            with open(file, 'r') as infile:
                lines = infile.readlines()
                if not header_written:
                    outfile.write(lines[0])  # Write header from the first file
                    header_written = True
                outfile.writelines(lines[1:])  # Write the rest of the lines, excluding the header

if __name__ == "__main__":
    input_directory = "."  # Replace with the path to your input directory
    output_filename = "./ERR6454470.txt"  # Replace with your desired output file path
    
    concatenate_blast_outputs(input_directory, output_filename)
    print(f"Concatenated BLAST results have been saved to {output_filename}")

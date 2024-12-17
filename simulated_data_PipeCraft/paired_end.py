import os
import subprocess

def run_command(command):
    """Run a shell command and handle errors."""
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e.cmd}\nExit code: {e.returncode}")
        exit(1)

def process_fastq(file, output_dir, read_length):
    """
    Generate a paired-end dataset from single-end reads.
    The forward reads remain as-is (5'-3'), and the reverse reads
    are produced by reverse complementing the forward reads.
    Both are trimmed to the same length to ensure equal read counts.
    Sequence IDs are kept identical to ensure compatibility with seqkit pair.
    """
    # Extract base name without the _R1.fastq suffix
    base_name = os.path.splitext(os.path.basename(file))[0].replace("_R1", "")

    # Define output file names
    paired_forward = os.path.join(output_dir, f"{base_name}_R1.fastq")
    paired_reverse = os.path.join(output_dir, f"{base_name}_R2.fastq")

    # Create a reverse complement file
    reverse_complement = os.path.join(output_dir, f"{base_name}_reverse_complement.fastq")
    print(f"Generating reverse complement for {file}...")
    run_command(f"seqkit seq -t DNA --reverse --complement {file} -o {reverse_complement}")

    # Trim both forward and reverse reads to target length
    forward_trimmed = os.path.join(output_dir, f"{base_name}_forward_trimmed.fastq")
    reverse_trimmed = os.path.join(output_dir, f"{base_name}_reverse_trimmed.fastq")

    print(f"Trimming forward reads for {file} to {read_length} bp...")
    run_command(f"seqkit subseq -r 1:{read_length} {file} -o {forward_trimmed}")

    print(f"Trimming reverse complement reads for {file} to {read_length} bp...")
    run_command(f"seqkit subseq -r 1:{read_length} {reverse_complement} -o {reverse_trimmed}")

    # Rename sequences for paired-end compatibility with identical sequence IDs
    print(f"Renaming sequences in {file} for paired-end compatibility...")
    temp_forward = os.path.join(output_dir, f"{base_name}_temp_R1.fastq")
    temp_reverse = os.path.join(output_dir, f"{base_name}_temp_R2.fastq")

    run_command(f"seqtk rename {forward_trimmed} {base_name} > {temp_forward}")
    run_command(f"seqtk rename {reverse_trimmed} {base_name} > {temp_reverse}")

    # Replace sequence IDs with consistent identifiers
    with open(temp_forward, 'r') as f_in, open(paired_forward, 'w') as f_out:
        for i, line in enumerate(f_in):
            if i % 4 == 0:
                f_out.write(f"@{base_name}:{i//4+1}\n")
            else:
                f_out.write(line)

    with open(temp_reverse, 'r') as f_in, open(paired_reverse, 'w') as f_out:
        for i, line in enumerate(f_in):
            if i % 4 == 0:
                f_out.write(f"@{base_name}:{i//4+1}\n")
            else:
                f_out.write(line)

    # Cleanup intermediate files
    os.remove(reverse_complement)
    os.remove(forward_trimmed)
    os.remove(reverse_trimmed)
    os.remove(temp_forward)
    os.remove(temp_reverse)

    print(f"Processed {file}: Final paired-end files saved as {paired_forward} and {paired_reverse}.")

def main():
    input_dir = "./single_end"
    output_dir = "paired_end_outputs"
    read_length = 150  # Target read length

    os.makedirs(output_dir, exist_ok=True)
    
    # Expecting single-end reads named like "sample_X_R1.fastq" as input
    fastq_files = [f for f in os.listdir(input_dir) if f.endswith(".fastq")]

    if not fastq_files:
        print("No FASTQ files found in the input directory.")
        return

    for fastq_file in fastq_files:
        full_path = os.path.join(input_dir, fastq_file)
        if os.path.exists(full_path):
            process_fastq(full_path, output_dir, read_length)
        else:
            print(f"File not found: {fastq_file}")

    print("All files processed successfully. Final paired-end files are ready.")

if __name__ == "__main__":
    main()
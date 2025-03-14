import os
import glob
from collections import defaultdict
import concurrent.futures
from multiprocessing import cpu_count

def read_fasta_file(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        current_id = None
        current_sequence = []
        
        for line in file:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_sequence)
                current_id = line[1:].strip().split()[0]
                current_sequence = []
            else:
                current_sequence.append(line.strip())
        
        if current_id:
            sequences[current_id] = ''.join(current_sequence)
    
    return file_path, sequences

def process_file(file_path):
    return read_fasta_file(file_path)

def merge_sequences(sequences_dict):
    merged_sequences = defaultdict(list)
    
    for file_sequences in sequences_dict.values():
        for seq_id, sequence in file_sequences.items():
            merged_sequences[seq_id].append(sequence)
    
    return {seq_id: ''.join(seqs) for seq_id, seqs in merged_sequences.items()}

def write_merged_fasta(output_file, merged_sequences):
    with open(output_file, 'w') as file:
        for seq_id, sequence in merged_sequences.items():
            file.write(f">{seq_id}\n{sequence}\n")

def main():
    # Specify the root directory containing the folders with FASTA files
    root_dir = '/path/to/root/directory'
    
    # Specify the output folder
    output_folder = os.path.join(root_dir, 'merged_output')
    os.makedirs(output_folder, exist_ok=True)
    
    # Get all FASTA files
    fasta_files = glob.glob(os.path.join(root_dir, '**/*.fasta'), recursive=True)
    
    # Group files by filename
    file_groups = defaultdict(list)
    for file_path in fasta_files:
        file_name = os.path.basename(file_path)
        file_groups[file_name].append(file_path)
    
    # Process each group of files
    for file_name, file_paths in file_groups.items():
        print(f"Merging files for {file_name}")
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count()) as executor:
            futures = [executor.submit(process_file, file_path) for file_path in file_paths]
            
            sequences_dict = {}
            for future in concurrent.futures.as_completed(futures):
                try:
                    file_path, sequences = future.result()
                    sequences_dict[file_path] = sequences
                except Exception as e:
                    print(f"Error processing {file_path}: {str(e)}")
        
        merged_sequences = merge_sequences(sequences_dict)
        
        output_file = os.path.join(output_folder, f"merged_{file_name}")
        try:
            write_merged_fasta(output_file, merged_sequences)
            print(f"Merged file written: {output_file}")
        except Exception as e:
            print(f"Error writing merged file: {str(e)}")

if __name__ == "__main__":
    main()

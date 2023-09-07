#!/bin/bash

# Define the number of chunks
NUM_CHUNKS=10

# Split the fasta file into equal chunks
awk -v prefix="OTUs_chunk_" '/^>/{n++;if(n%'$NUM_CHUNKS'==1) close(f); f=prefix n ".fasta"} {print > f}' OTUs.fasta

# Loop through each chunk and submit a separate job for BLAST
for file in OTUs_chunk_*.fasta; do
    cat <<EOT > run_blast_$file.sh
#!/bin/bash
#SBATCH --job-name="blast_$file"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --partition amd
#SBATCH --time=10:00:00

blastn -query $file \
-db silva_138 \
-word_size 7 \
-num_threads 128 \
-outfmt "6  delim=; qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" \
-evalue 0.001 \
-strand both \
-max_target_seqs 1 \
-max_hsps 2 \
-dust no \
-soft_masking true \
-penalty -1 \
-reward 1 \
-gapopen 1 \
-gapextend 2 \
-out ${file%.fasta}_blast_results.txt
EOT

    # Grant execute permissions and submit the job
    chmod +x run_blast_$file.sh
    sbatch run_blast_$file.sh
done

# After all jobs are complete, merge the results
cat OTUs_chunk_*_blast_results.txt > final_blast_results.txt
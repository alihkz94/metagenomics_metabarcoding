#!/bin/bash
#SBATCH --job-name="blast"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --partition amd
#SBATCH --time=30:00:00

## for top 10 hit
blastn -query /gpfs/space/home/alihakim/analysis/2014_cluster/OTUs.fasta \
-db /gpfs/space/home/alihakim/analysis/databse/UNITE \
-word_size 7 \
-num_threads 128 \
-outfmt "6  delim=; qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" \
-evalue 0.001 \
-strand both \
-max_target_seqs 10 \
-max_hsps 2 \
-dust no \
-soft_masking true \
-penalty -1 \
-reward 1 \
-gapopen 1 \
-gapextend 2 \
-out blast_results.txt

##for first top hit

#!/bin/bash
#SBATCH --job-name="tophit"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --partition amd
#SBATCH --time=50:00:00

blastn -query /gpfs/space/home/alihakim/analysis/2014_cluster/OTUs.fasta \
-db /gpfs/space/home/alihakim/analysis/databse/UNITE \
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
-out blast_tophit_results.txt

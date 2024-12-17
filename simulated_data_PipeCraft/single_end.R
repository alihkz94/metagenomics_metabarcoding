# This script generates realistic but reproducible single-end FASTQ simulation data.
# It ensures that DADA2 can consistently learn error rates and yield identical ASVs each run.
# Key points:
# - Fixed RNG seed for reproducibility.
# - Enough references, reads, and slight mutations ensure complexity.
# - Stable abundance distributions and controlled quality variation enable consistent DADA2 results.

set.seed(123)  # Fixed seed for deterministic results

library(Biostrings)
library(ShortRead)

# Parameters
num_ref_seqs <- 50       # Enough references to ensure complexity
seq_length <- 250
num_samples <- 20
reads_per_sample <- 10000
base_probs <- c(A=0.25, C=0.25, G=0.25, T=0.25)
mutation_rate <- 0.005    # Small mutation rate for slight variability
mean_abundance <- 1
abundance_sd <- 0.5
qual_mean_start <- 35
qual_mean_end <- 25
qual_sd <- 3
output_dir <- "simulated_fastq_test"

if(!dir.exists(output_dir)) dir.create(output_dir)

# Generate Reference Sequences
generate_reference_seqs <- function(num, length, probs) {
  refs <- replicate(num, paste0(sample(names(probs), size=length, replace=TRUE, prob=probs), collapse=""))
  DNAStringSet(refs)
}

ref_seqs <- generate_reference_seqs(num_ref_seqs, seq_length, base_probs)
names(ref_seqs) <- paste0("ref_", seq_len(num_ref_seqs))

# Generate Abundances
generate_abundances <- function(num_refs, mean_ab, sd_ab) {
  raw <- rnorm(num_refs, mean=mean_ab, sd=sd_ab)
  raw[raw < 0] <- runif(sum(raw<0), 0.1, 0.5) 
  rel_ab <- raw / sum(raw)
  rel_ab
}

# Mutate Sequence
mutate_sequence <- function(seq, mut_rate, probs) {
  bases <- unlist(strsplit(seq, ""))
  mutate_positions <- which(runif(length(bases)) < mut_rate)
  if(length(mutate_positions) > 0) {
    new_bases <- sample(names(probs), length(mutate_positions), replace=TRUE, prob=probs)
    for(i in seq_along(mutate_positions)) {
      while(new_bases[i] == bases[mutate_positions[i]]) {
        new_bases[i] <- sample(names(probs), 1, prob=probs)
      }
    }
    bases[mutate_positions] <- new_bases
  }
  paste0(bases, collapse="")
}

# Generate Qualities with a gradient
generate_qualities <- function(read_length, q_start_mean, q_end_mean, q_sd) {
  qslope <- (q_end_mean - q_start_mean)/read_length
  q_means <- q_start_mean + qslope*(0:(read_length-1))
  q_vals <- rnorm(read_length, mean=q_means, sd=q_sd)
  q_vals[q_vals < 2] <- 2
  q_vals[q_vals > 41] <- 41
  q_vals
}

simulate_sample_reads <- function(refs, reads_per_sample, mutation_rate, probs, 
                                  q_start_mean, q_end_mean, q_sd) {
  rel_ab <- generate_abundances(length(refs), mean_abundance, abundance_sd)
  counts <- round(rel_ab * reads_per_sample)
  
  all_reads <- DNAStringSet()
  all_quals <- BStringSet()
  
  for(i in seq_along(refs)) {
    if(counts[i] == 0) next
    base_seq <- as.character(refs[i])
    mutated_reads <- replicate(counts[i], mutate_sequence(base_seq, mutation_rate, base_probs))
    mutated_reads_set <- DNAStringSet(mutated_reads)
    
    q_matrix <- sapply(seq_len(counts[i]), function(x) {
      q_values <- generate_qualities(width(refs[i]), q_start_mean, q_end_mean, q_sd)
      rawToChar(as.raw(q_values+33))
    })
    
    all_reads <- c(all_reads, mutated_reads_set)
    all_quals <- c(all_quals, BStringSet(q_matrix))
  }
  
  ids <- BStringSet(sprintf("read_%05d", seq_len(length(all_reads))))
  srq <- ShortReadQ(sread=all_reads, quality=FastqQuality(all_quals), id=ids)
  srq
}

# Simulate all samples
for(s in 1:num_samples) {
  cat("Simulating sample", s, "...\n")
  srq <- simulate_sample_reads(
    refs = ref_seqs,
    reads_per_sample = reads_per_sample,
    mutation_rate = mutation_rate,
    probs = base_probs,
    q_start_mean = qual_mean_start,
    q_end_mean = qual_mean_end,
    q_sd = qual_sd
  )
  
  out_file <- file.path(output_dir, paste0("sample_", s, "_single_end.fastq"))
  writeFastq(srq, file=out_file, compress=FALSE)
  cat("Wrote:", out_file, "\n")
}

cat("Simulation complete. The generated data should support stable, repeatable DADA2 analysis.\n")

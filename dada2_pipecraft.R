# Load necessary libraries
library(BiocGenerics)
library(Biostrings)
library(Rcpp)
library(dada2)
library(ShortRead)
library(tibble)
library(seqinr)
library(Hmisc)
library(digest)
library(fs)

# Directory
setwd("/media/ali/data_store/test_dada2")

# List all fastq files in the working directory
fnFs <- list.files(pattern = ".fastq", full.names = TRUE)

# Read the original FASTQ files
original_reads <- lapply(fnFs, readFastq)

# Convert to DNAStringSet objects
original_dna <- lapply(original_reads, sread)

# Extract the headers
original_headers <- lapply(original_dna, names)

# Define file paths for filtered files
filtFs <- file.path("filtered", basename(fnFs))

# Quality filtering
out <- filterAndTrim(fnFs, filtFs, truncLen = 0, maxN = 0, maxEE = 2, truncQ = 2, minQ = 3, rm.phix = FALSE, multithread = TRUE)
cat("Filtered", sum(out), "reads from", length(out), "sample(s).\n")

# Check if all filtered files exist
missing_files <- filtFs[!file.exists(filtFs)]
if(length(missing_files) > 0) {
  warning("The following filtered files do not exist:\n", paste(missing_files, collapse = "\n"))
}

# Update filtFs only to include existing files
filtFs <- filtFs[file.exists(filtFs)]

# Determine the number of subsets to process
n_subsets <- 10  # Adjust this number based on your system's capacity

# Determine the number of files per subset
files_per_subset <- ceiling(length(filtFs) / n_subsets)

# Initialize an empty list to store the dereplicated data
derep_list <- vector("list", n_subsets)

# Loop through each subset of files
for(i in seq_len(n_subsets)){
  
  # Determine the file indices for this subset
  idx <- ((i - 1) * files_per_subset + 1):min(i * files_per_subset, length(filtFs))
  
  # Select the files for this subset
  subset_files <- filtFs[idx]
  
  # Ensure that subset_files is a character vector
  if (!is.character(subset_files)) {
    subset_files <- as.character(subset_files)
  }
  
  # Dereplicate this subset of files
  derep_list[[i]] <- derepFastq(subset_files, verbose = TRUE)
  
}

# Combine all the dereplicated data into a single list
derepFs <- do.call(c, derep_list)

# Create sequence table directly from dereplicated data
seqtab <- makeSequenceTable(derepFs)

# Inspect sequence table
print(dim(seqtab))

# Optionally, save the sequence table for later analysis
saveRDS(seqtab, "sequence_table.rds")

# Load the sequence table if needed
seqtab <- readRDS("sequence_table.rds")

# Existing code for chimera filtering
nonchimeric_seqtab <- removeBimeraDenovo(seqtab, minSampleFraction = 0.9, ignoreNNegatives = 1, method = "consensus",
                                         minFoldParentOverAbundance = 1.5, minParentAbundance = 2, 
                                         allowOneOff = TRUE, minOneOffParentDistance = 2, maxShift = 16, 
                                         multithread = TRUE, verbose = TRUE)
# Extract sequences
nonchimeric_sequences <- as.character(getSequences(nonchimeric_seqtab))
discarded_chimera <- seqtab[,colnames(seqtab) %nin% colnames(nonchimeric_seqtab)]
chimeric_sequences <- as.character(getSequences(discarded_chimera))

# Ensure directory existence
ensureDirExists <- function(dir_name) {
  if (!dir_exists(dir_name)) {
    dir_create(dir_name)
  }
}

# Write sequences to FASTA with counts, excluding sequences with zero counts
writeFastaWithCounts <- function(sequences, counts, headers, dir_name, file_name) {
  ensureDirExists(dir_name)
  file_path <- file.path(dir_name, file_name)
  
  # Initialize a vector to hold FASTA content
  fasta_lines <- c()
  
  # Loop through sequences and prepare FASTA content
  for (i in seq_along(sequences)) {
    if (counts[i] > 0) {  # Check if count is greater than zero
      # Modified header format to use ";size=" instead of "_count_"
      fasta_header <- paste(">", headers[i], ";size=", counts[i], sep = "")
      fasta_sequence <- sequences[i]
      fasta_lines <- c(fasta_lines, fasta_header, fasta_sequence)
    }
  }
  
  # Write to file only if there's content to avoid empty files
  if (length(fasta_lines) > 0) {
    con <- file(file_path, "w")
    writeLines(fasta_lines, con)
    close(con)
  }
}

# Process and organize sequences by sample, respecting counts
processAndWriteSequences <- function(seqtab, isChimeric = FALSE) {
  samples <- unique(gsub("\\.fastq\\.gz$", "", rownames(seqtab)))
  dir_name <- ifelse(isChimeric, "chimeras", "non_chimeric")
  
  # Loop through each sample
  for (sample in samples) {
    # Subset for the current sample; ensure drop = FALSE to keep data frame structure
    sample_seqtab <- seqtab[grep(paste0("^", sample), rownames(seqtab)), , drop = FALSE]
    
    # Skip processing for empty or all-zero-count samples
    if (ncol(sample_seqtab) == 0 || all(colSums(sample_seqtab) == 0)) next
    
    # Extracting sequences and counts
    sequences <- colnames(sample_seqtab)
    counts <- colSums(sample_seqtab)
    valid_indices <- counts > 0  # Indices where count > 0
    headers <- sapply(sequences[valid_indices], function(x) digest(x, algo = "sha1"))
    
    # File naming convention
    file_name <- paste0(sample, ifelse(isChimeric, "_chimeric.fasta", ".fasta"))
    
    # Write sequences with valid counts to FASTA
    writeFastaWithCounts(sequences[valid_indices], counts[valid_indices], headers, dir_name, file_name)
  }
}

# Assuming 'nonchimeric_seqtab' and 'discarded_chimera' are your data
processAndWriteSequences(nonchimeric_seqtab, isChimeric = FALSE)
processAndWriteSequences(discarded_chimera, isChimeric = TRUE)

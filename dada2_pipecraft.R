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


# Set working directory to where your fastq files are stored
setwd("/media/ali/data_store/test_dada2")

# List all fastq files in the working directory
fnFs <- list.files(pattern = ".fastq", full.names = TRUE)

# Read and preprocess the FASTQ files
original_reads <- lapply(fnFs, readFastq)
original_dna <- lapply(original_reads, sread)
original_headers <- lapply(original_dna, names)

# Filter reads and setup filtered file paths
filtFs <- file.path("filtered", basename(fnFs))
out <- filterAndTrim(fnFs, filtFs, truncLen = 0, maxN = 0, maxEE = 2, truncQ = 2, minQ = 3, rm.phix = FALSE, multithread = TRUE)

# Dereplication
derepFs <- derepFastq(filtFs, verbose = TRUE)

# Create sequence table
seqtab <- makeSequenceTable(derepFs)

# Chimera removal
nonchimeric_seqtab <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose=TRUE)

# Extract sequences
nonchimeric_sequences <- as.character(getSequences(nonchimeric_seqtab))
discarded_chimera <- seqtab[,colnames(seqtab) %nin% colnames(nonchimeric_seqtab)]
chimeric_sequences <- as.character(getSequences(discarded_chimera))

# Function to generate SHA1 hashes
generateSHA1 <- function(sequence) {
  digest(sequence, algo = "sha1", serialize = FALSE)
}

# Function to create directories if they don't already exist
ensureDirExists <- function(dir_name) {
  if (!dir_exists(dir_name)) {
    dir_create(dir_name)
  }
}

# Modified function to write sequences to FASTA files with directory support
writeFastaWithCounts <- function(sequences, counts, headers, dir_name, file_name) {
  ensureDirExists(dir_name) # Ensure the directory exists
  file_path <- file.path(dir_name, file_name)
  
  fasta_lines <- vector("list", length(sequences))
  for (i in seq_along(sequences)) {
    fasta_header <- paste(">", headers[i], "_count_", counts[i], sep = "")
    fasta_sequence <- sequences[i]
    fasta_lines[[i]] <- paste(fasta_header, fasta_sequence, sep = "\n")
  }
  
  con <- file(file_path, "w")
  writeLines(unlist(fasta_lines), con)
  close(con)
}

# Adjusted function to process and organize sequences by sample with directory support
processAndWriteSequences <- function(seqtab, isChimeric = FALSE) {
  samples <- unique(gsub("__.+$", "", rownames(seqtab)))
  dir_name <- ifelse(isChimeric, "chimeras", "non_chimeric")
  
  for (sample in samples) {
    sample_indices <- grepl(paste0("^", sample, "__"), rownames(seqtab))
    sample_seqtab <- seqtab[sample_indices, ]
    sequence_counts <- colSums(sample_seqtab)
    sequences <- names(sequence_counts)
    counts <- as.integer(sequence_counts)
    headers <- sapply(sequences, generateSHA1, USE.NAMES = FALSE)
    
    file_name <- paste0(sample, ifelse(isChimeric, "_chimeric.fasta", ".fasta"))
    
    writeFastaWithCounts(sequences, counts, headers, dir_name, file_name)
  }
}

# Create directories and process sequences
processAndWriteSequences(nonchimeric_seqtab, isChimeric = FALSE)
processAndWriteSequences(discarded_chimera, isChimeric = TRUE)

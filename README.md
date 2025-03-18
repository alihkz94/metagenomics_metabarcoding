# Metabarcoding and Metagenomics Analysis Tools

A comprehensive toolkit for environmental DNA analysis, containing scripts and workflows for both metabarcoding and metagenomics analysis on High-Performance Computing (HPC) environments.

## Overview

This repository contains a collection of scripts, pipelines, and utilities for processing and analyzing environmental DNA data. It supports both targeted metabarcoding approaches and whole-genome shotgun metagenomics analysis, with tools organized by function and technology.

## Repository Structure
Collecting workspace informationI'll create a comprehensive README file that covers both metabarcoding and metagenomics analysis components of your repository:


# Metabarcoding and Metagenomics Analysis Tools

A comprehensive toolkit for environmental DNA analysis, containing scripts and workflows for both metabarcoding and metagenomics analysis on High-Performance Computing (HPC) environments.

## Overview

This repository contains a collection of scripts, pipelines, and utilities for processing and analyzing environmental DNA data. It supports both targeted metabarcoding approaches and whole-genome shotgun metagenomics analysis, with tools organized by function and technology.

## Repository Structure

```plaintext
.
├── LICENSE                        # MIT License
├── README.md                      # This file
├── HPC_scripts/                   # Scripts for various analysis steps
│   ├── abundance.sh               # Calculate sequence abundances
│   ├── BlasCh_pipecraft.py        # Chimera detection and recovery
│   ├── concatanate.sh             # Sequence concatenation
│   ├── custom pipeline.sh         # Custom workflow for metabarcoding
│   ├── dada2_pipecraft.R          # DADA2 pipeline implementation
│   ├── demultiplexing_indexcheck.sh # Check duplicate headers and indices
│   ├── dereplication.sh           # Remove duplicate sequences
│   ├── fastq2fasta.sh             # Convert FASTQ to FASTA format
│   ├── ITSx.sh                    # Extract ITS regions
│   ├── mafft.sh                   # Sequence alignment
│   ├── makedb.sh                  # Create reference databases
│   ├── merge_same_fasta.py        # Combine FASTA files
│   ├── OTU_singelton_removal.sh   # Remove singleton OTUs
│   ├── parsing_issue.sh           # Handle FASTA parsing issues
│   ├── sintax.sh                  # Taxonomic classification with SINTAX
│   ├── total_abudances.sh         # Calculate total abundances
│   ├── usearch.sh                 # OTU clustering and chimera removal
│   └── blast/                     # BLAST utilities
│       ├── blast_PC.sh            # Parallel BLAST implementation
│       ├── blast_folder.sh        # Process multiple directories
│       ├── blast_general.sh       # General purpose BLAST script
│       └── concat.py              # Concatenate BLAST results
├── nextpac_pipeline/              # Nextflow pipelines for long-read metagenomics
│   ├── MAG_flye.nf               # MAG assembly with Flye 
│   ├── MetaFusion.nf             # Metagenome assembly and analysis
│   ├── nextflow_flye.config      # Configuration for Flye pipeline
│   └── nextflow_spades.config    # Configuration for SPAdes pipeline
├── simulated_data_PipeCraft/      # Tools for generating simulated datasets
│   ├── paired_end.py             # Generate paired-end reads
│   ├── Readme.md                 # Documentation for simulation tools
│   └── single_end.R              # Generate single-end reads
    └── simulate_DADA2.R          # Gnerate reads suitable for DADA2 testing
└── Snake_EUK/                    # Snakemake pipeline for taxonomy processing
    ├── readme.md                 # Documentation for Snakemake workflow
    ├── snakefile                 # Snakemake workflow definition
    ├── config.yaml               # Configuration parameters
    └── scripts/                  # Python scripts for different tools
        ├── dada2.py              # DADA2 format conversion
        ├── general.py            # General taxonomy processing
        ├── mothur.py             # Mothur format conversion
        ├── qiime.py              # QIIME2 format conversion
        ├── sintax.py             # SINTAX format conversion
        └── utils.py              # Shared utility functions
```

## Key Features

### Metabarcoding Analysis
- **Data Preprocessing**: Quality filtering, primer trimming, and format conversion
- **Sequence Processing**: Dereplication, chimera detection, and ITS extraction
- **OTU/ASV Generation**: DADA2 implementation and clustering tools
- **Taxonomy Assignment**: SINTAX and BLAST-based classification

### Metagenomics Analysis
- **Assembly Pipelines**: Nextflow workflows for short and long read assembly
- **MAG Recovery**: Tools for binning and refining metagenome-assembled genomes
- **Taxonomic Profiling**: Multiple approaches for community profiling
- **Functional Analysis**: Integration with annotation tools

### Data Simulation
- Generation of realistic sequencing data for benchmarking and testing
- Support for both single-end and paired-end read simulation
- Parameter documentation for customizing simulations

### Workflow Management
- Snakemake pipeline for consistent taxonomy handling across tools
- Nextflow workflows for complex metagenomics analysis
- Robust error handling and reporting

## Installation & Dependencies

### General Requirements
- **Bash** environment (Linux/MacOS/WSL)
- **Python 3.6+** with Biopython
- **R** with BiocManager for DADA2-related scripts

### Specialized Tools
- **BLAST+**: For sequence similarity searches
- **VSEARCH**: For various sequence operations
- **ITSx**: For ITS region extraction
- **Seqkit**: For FASTA/FASTQ manipulation

### Workflow Systems
- **Snakemake**: For taxonomy pipeline
- **Nextflow**: For metagenomics workflows

## Usage

### HPC Scripts

Most scripts in the `HPC_scripts/` directory can be run directly:

```bash
bash HPC_scripts/fastq2fasta.sh
```

For more complex scripts, check the header of each file for usage instructions.

### Nextflow Pipelines

Navigate to the nextpac_pipeline directory and run:

```bash
nextflow run MAG_flye.nf -c nextflow_flye.config
# or
nextflow run MetaFusion.nf -c nextflow_spades.config
```

### Snakemake Taxonomy Pipeline

Navigate to the Snake_EUK directory and run:

```bash
snakemake --cores <number_of_cores> --rerun-incomplete --keep-going
```

### Data Simulation

For generating test datasets:

```bash
# For single-end reads
Rscript simulated_data_PipeCraft/single_end.R

# For paired-end reads
python simulated_data_PipeCraft/paired_end.py

# For DADA2 specific reads
Rscript simulated_data_PipeCraft/simulate_DADA2.R
```

## Documentation

Each directory contains specific documentation:
- Simulation Parameters
- Taxonomy Pipeline

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author

Ali Hakimzadeh

## Citation

If you use these tools in your research, please cite:
```
Hakimzadeh A. (2025). Metabarcoding and Metagenomics Analysis Tools. GitHub repository. https://github.com/alihkz94/metagenomics_metabarcoding
```

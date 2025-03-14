# Nextpac Pipeline

A comprehensive Nextflow-based pipeline for metagenomic analysis, including assembly, binning, and taxonomic annotation of metagenomic sequences using various bioinformatics tools.

## Overview

The Nextpac Pipeline provides two main workflows:

1. **MAG_flye.nf**: A workflow for metagenome assembly and binning using Flye for long reads
2. **MetaFusion.nf**: A pipeline for quality trimming, assembly, binning, and taxonomic annotation of metagenomic sequences using metaSPAdes

Both pipelines utilize multiple binning approaches and integrate results using DAS Tool for improved metagenomic analysis.

## Repository Structure

```plaintext
.
├── MAG_flye.nf                 # Nextflow script for long-read metagenome assembly and binning
├── MetaFusion.nf               # Nextflow script for short-read metagenome assembly and binning
├── nextflow_flye.config        # Configuration file for the Flye-based pipeline
├── nextflow_spades.config      # Configuration file for the metaSPAdes-based pipeline
└── README.md                   # This README file
```

## Dependencies

The pipeline requires the following tools, which are managed through Conda environments:

- **Flye**: For long-read assembly
- **metaSPAdes**: For short-read assembly
- **Trimmomatic**: For quality trimming
- **minimap2**: For read mapping
- **SemiBin**: For binning with long-read data
- **MetaBAT2**: For binning
- **MaxBin2**: For binning
- **CONCOCT**: For binning
- **DAS Tool**: For integrating binning results
- **CAT/BAT**: For taxonomic annotation
- **ITSx**: For extracting fungal ITS regions
- **Kaiju**: For taxonomic classification of reads

## Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/alihkz94/metabarcoding_analysis.git
   cd metabarcoding_analysis/nextpac_pipeline
   ```

2. **Install Nextflow** (if not already installed):
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

3. **Install Conda** (if not already installed):
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

4. **Create Conda environments** for each tool as specified in the config files:
   ```bash
   # Example for creating Flye environment
   conda create -n flye -c bioconda flye
   
   # Example for creating metabat2 environment
   conda create -n metabat2 -c bioconda metabat2
   
   # Continue creating environments for other tools
   ```

## Usage

### Long-read Pipeline (MAG_flye.nf)

This pipeline is designed for processing long-read sequencing data (e.g., Oxford Nanopore).

```bash
nextflow run MAG_flye.nf -c nextflow_flye.config -profile conda \
  --reads "path/to/your/long_reads.fastq" \
  --outdir "path/to/output/directory"
```

### Short-read Pipeline (MetaFusion.nf)

This pipeline is designed for processing short-read sequencing data (e.g., Illumina).

```bash
nextflow run MetaFusion.nf -c nextflow_spades.config -profile conda \
  --reads "path/to/your/*_R{1,2}.fastq" \
  --outdir "path/to/output/directory"
```

## Pipeline Details

### MAG_flye.nf

A workflow for metagenome assembly and binning using tools for long-read sequencing data.

**Key Processes**:
- **Flye**: Assembles high-fidelity reads
- **minimap2**: Maps reads to the assembly
- **SemiBin**: Bins the assembly using long-read sequencing data
- **depth**: Summarizes contig depths
- **metabat**: Bins the assembly using MetaBAT2
- **table**: Generates contig-to-bin mappings
- **dastool**: Integrates binning results
- **BAT**: Taxonomically annotates bins
- **CAT**: Taxonomically annotates contigs
- **ITSx**: Extracts fungal ITS regions

### MetaFusion.nf

A comprehensive pipeline for quality trimming, assembly, binning, and taxonomic annotation of short-read metagenomic sequences.

**Key Processes**:
- **trimmomatic**: Quality trims the reads
- **assembly**: Assembles the reads using metaSPAdes
- **maxbin**: Bins the assembly using MaxBin2
- **minimap2**: Maps reads to the assembly
- **depth**: Estimates contig depths
- **metabat**: Bins the assembly using MetaBAT2
- **concoct**: Bins the assembly using CONCOCT
- **table**: Generates contig-to-bin mappings
- **dastool**: Integrates binning results
- **BAT**: Taxonomically annotates bins using CAT
- **CAT**: Taxonomically annotates contigs using CAT
- **ITSx**: Extracts fungal ITS regions
- **kaiju**: Classifies reads taxonomically

## Configuration Files

### nextflow_flye.config

Configuration file for running the Nextflow pipeline with Flye, specifying conda environments for various processes:

- **Flye**: Uses the Flye conda environment
- **SemiBin**: Uses the SemiBin conda environment
- **depth**: Uses the metabat2 conda environment
- **metabat**: Uses the metabat2 conda environment
- **BAT**: Uses the CAT conda environment
- **CAT**: Uses the CAT conda environment

### nextflow_spades.config

Configuration file for running the Nextflow pipeline with metaSPAdes, specifying conda environments for various processes:

- **metabat**: Uses the metabat2 conda environment
- **concoct**: Uses the metaWRAP conda environment
- **BAT**: Uses the CAT conda environment
- **CAT**: Uses the CAT conda environment
- **bins_abundance**: Uses the metaWRAP conda environment
- **metaxa**: Uses the metaxa conda environment

## Expected Outputs

The pipeline generates various outputs, including:

- Assembled contigs
- Binned metagenome-assembled genomes (MAGs)
- Taxonomic classifications of contigs and bins
- Depth profiles of contigs
- Extracted ITS regions for fungal identification

## Troubleshooting

- **Memory Issues**: If you encounter memory issues, try adjusting the memory allocation in the Nextflow configuration files.
- **Conda Environment Issues**: Ensure all required conda environments are properly created and activated.
- **Input Data Format**: Verify that your input data is in the correct format as expected by the pipeline.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

This pipeline was developed to streamline metagenomic analysis workflows. Contributions, suggestions, and bug reports are welcome.

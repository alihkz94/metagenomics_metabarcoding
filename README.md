# Metabarcoding and Metagenomics Analysis Tools

A comprehensive toolkit for environmental DNA analysis, containing scripts, workflows, and utilities for both metabarcoding and metagenomics on High-Performance Computing (HPC) systems.

## Overview

This repository includes tools for data preprocessing, chimera detection, OTU/ASV generation, taxonomic classification, assembly workflows, and simulated dataset generation. It is designed for environmental DNA and metabarcoding analyses, with a focus on HPC-friendly scripts and workflow automation.

## Repository Structure

```plaintext
.
├── BlasCh/                        # BlasCh package for chimera recovery
│   ├── BlasCh_HPC/                # Scripts for running BlasCh in HPC environments
│   │   ├── blasch.sh
│   │   ├── blasch_chimeric.py
│   │   └── blasch_nonchimeric.py
│   ├── __init__.py
│   ├── main.py
│   ├── readme.md
│   └── setup.py
├── Blasch_PipeCraft.py            # Chimera detection and recovery script
├── HPC_scripts/                   # General HPC scripts for metabarcoding workflows
│   ├── abundance.sh
│   ├── BlasCh_pipecraft.py
│   ├── concatanate.sh
│   ├── custom pipeline.sh
│   ├── dada2_pipecraft.R
│   ├── demultiplexing_indexcheck.sh
│   ├── dereplication.sh
│   ├── fastq2fasta.sh
│   ├── ITSx.sh
│   ├── mafft.sh
│   ├── makedb.sh
│   ├── merge_same_fasta.py
│   ├── OTU_singelton_removal.sh
│   ├── parsing_issue.sh
│   ├── sintax.sh
│   ├── total_abudances.sh
│   ├── usearch.sh
│   └── blast/
│       ├── blast_PC.sh
│       ├── blast_folder.sh
│       ├── blast_general.sh
│       └── concat.py
├── LICENSE
├── README.md
├── nextpac_pipeline/              # Nextflow-based metagenomics workflows
│   ├── MAG_flye.nf
│   ├── MetaFusion.nf
│   ├── nextflow_flye.config
│   └── nextflow_spades.config
└── simulated_data_PipeCraft/      # Simulated dataset generation tools
    ├── paired_end.py
    ├── Readme.md
    ├── simulate_DADA2.R
    └── single_end.R
```

## Key Features

### Metabarcoding Analysis
- Data preprocessing and format conversion
- Chimera detection and recovery
- Dereplication and OTU/ASV generation
- ITS extraction and SINTAX taxonomy assignment
- BLAST-based classification and reporting

### Metagenomics Analysis
- Nextflow assembly pipelines for long-read and short-read datasets
- MAG assembly support with Flye and SPAdes configurations
- Modular workflow configuration for HPC environments

### Chimera Recovery
- `Blasch_PipeCraft.py` for BLAST XML chimera classification and recovery
- `BlasCh/` package for false-positive chimera detection and recovery
- HPC wrapper scripts in `BlasCh/BlasCh_HPC/`

### Data Simulation
- Single-end and paired-end read generation
- DADA2-compatible simulated read sets
- Useful for benchmarking and pipeline development

## Installation & Dependencies

### General Requirements
- Linux or macOS shell environment
- Python 3.6+ with Biopython
- R (for DADA2 and simulation scripts)

### Recommended Tools
- BLAST+
- VSEARCH
- ITSx
- Seqkit
- Nextflow

## Usage

### Run HPC scripts

Example:

```bash
bash HPC_scripts/fastq2fasta.sh
```

### Run Nextflow pipelines

```bash
nextflow run nextpac_pipeline/MAG_flye.nf -c nextpac_pipeline/nextflow_flye.config
# or
nextflow run nextpac_pipeline/MetaFusion.nf -c nextpac_pipeline/nextflow_spades.config
```

### Run Blasch chimera recovery

```bash
python Blasch_PipeCraft.py --help
```

### Generate simulated data

```bash
Rscript simulated_data_PipeCraft/single_end.R
python simulated_data_PipeCraft/paired_end.py
Rscript simulated_data_PipeCraft/simulate_DADA2.R
```

## Documentation

- `BlasCh/readme.md` for the BlasCh package
- `simulated_data_PipeCraft/Readme.md` for simulation tool details
- Inline comments and usage headers in `HPC_scripts/`

## License

This project is licensed under the MIT License. See `LICENSE` for details.

## Author

Ali Hakimzadeh

## Citation

If you use these tools in your research, please cite:

```text
Hakimzadeh A. (2025). Metabarcoding and Metagenomics Analysis Tools. GitHub repository. https://github.com/alihkz94/metagenomics_metabarcoding
```

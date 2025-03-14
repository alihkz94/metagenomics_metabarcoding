# Eukaryome Taxonomy Pipeline

A comprehensive and reproducible pipeline for processing taxonomic FASTA files and converting them into formats suitable for multiple downstream tools including General taxonomy, DADA2, Mothur, QIIME2, and SINTAX.

## Overview

This pipeline is designed to process and clean FASTA files containing taxonomic information. It performs robust encoding conversion (from latin‑1 to ASCII) on the fly, filters and cleans taxonomy headers, and generates outputs tailored for various bioinformatics tools. The workflow is orchestrated via [Snakemake](https://snakemake.readthedocs.io/) to ensure reproducibility and parallel processing. Each tool-specific script is modular and can be adjusted to meet specific requirements.

## Repository Structure

```plaintext
.
├── Snakefile                 # Snakemake workflow file
├── config.yaml               # Pipeline configuration (versioning, etc.)
├── scripts
│   ├── utils.py              # Utility module for robust file handling (encoding conversion)
│   ├── general.py            # Generates General taxonomy formatted FASTA files
│   ├── dada2.py              # Converts FASTA headers for DADA2 (removes accession numbers)
│   ├── mothur.py             # Generates Mothur-compatible FASTA and TAX files
│   ├── qiime.py              # Generates QIIME2-compatible FASTA and TSV files with cleaned taxonomy fields
│   └── sintax.py             # Converts FASTA files to SINTAX format
└── README.md                 # This README file
```

## Installation & Dependencies

### Requirements

- **Python 3.6+**
- [Biopython](https://biopython.org/)  
  Install via: `pip install biopython`
- [Snakemake](https://snakemake.readthedocs.io/)  
  Install via: `pip install snakemake`
- [Seqkit](https://bioinf.shenwei.me/seqkit/)  
  Install via Conda: `conda install -c bioconda seqkit` or follow instructions on the [Seqkit website](https://bioinf.shenwei.me/seqkit/)

### Installation

Clone the repository and install Python dependencies:

```bash
git clone https://github.com/alihkz94/metagenomics_metabarcoding.git
cd metagenomics_metabarcoding/Snake_EUK
pip install biopython snakemake
```

Ensure that `seqkit` is installed and available in your system PATH.

## Usage

### Input Files

Place your input FASTA files (e.g., `ITS.fasta`, `LSU.fasta`, `SSU.fasta`, `longread.fasta`) in the repository root (or designated input folder).

### Configuration

Edit `config.yaml` to set the desired version string, for example:

```yaml
version: "1.9.4"
```

### Running the Pipeline

Run the entire pipeline with:

```bash
snakemake --cores <number_of_cores> --rerun-incomplete --keep-going
```

To run only for instance the DADA2 conversion over your original FASTA files, execute:

```bash
snakemake dada2/DADA2_EUK_ITS_v1.9.4.fasta dada2/DADA2_EUK_LSU_v1.9.4.fasta dada2/DADA2_EUK_SSU_v1.9.4.fasta dada2/DADA2_EUK_longread_v1.9.4.fasta --rerun-incomplete --keep-going --cores 8
```

> **Tip:** If you suspect output files are outdated or incorrect, you can force re-run of jobs using `--forceall` or `--forcerun <target>`.

### Output Files

- **General:** `general/General_EUK_{base}_v{version}.fasta`
- **DADA2:** `dada2/DADA2_EUK_{base}_v{version}.fasta`  
  *Note:* The DADA2 script removes the accession number and truncates the header at the first placeholder.
- **Mothur:** `mothur/mothur_EUK_{base}_v{version}.fasta` and `mothur/mothur_EUK_{base}_v{version}.tax`
- **QIIME2:** `qiime2/QIIME2_EUK_{base}_v{version}.fasta` and `qiime2/QIIME2_EUK_{base}_v{version}.tsv`  
  *Note:* Taxonomy in the TSV files is cleaned to remove trailing dot placeholders.
- **SINTAX:** `sintax/SINTAX_EUK_{base}_v{version}.fasta`

## Pipeline Details

- **Robust Encoding Conversion:**  
  All scripts utilize a custom file-like wrapper (implemented in `utils.py`) that reads files using `latin-1` decoding and converts them on the fly to ASCII. This ensures that all non-ASCII characters are handled gracefully.

- **Taxonomy Header Cleaning:**  
  The scripts are designed to remove unwanted placeholders (`.`) and extra taxonomic levels.  
  - For DADA2, the script removes the accession number and retains taxonomy fields only until the first placeholder.
  - For QIIME2 and Mothur, the scripts clean the TSV taxonomy output by stripping trailing dot placeholders.

- **Reproducible Workflow:**  
  The entire process is managed by Snakemake, ensuring that jobs run in the correct order with proper dependencies, even when running in parallel.

## Customization & Troubleshooting

- **Modifying Header Formatting:**  
  The header processing logic is contained within each script (e.g., `dada2.py`, `qiime.py`). You can modify the functions `transform_header` or `clean_taxonomy` as needed.

- **Resource Management:**  
  If you encounter memory or process issues (e.g., SIGKILLs), try reducing the number of cores with `--cores 2` or increasing system resources.

- **Debugging:**  
  Use Snakemake's verbose and print shell command options (`--printshellcmds`) for detailed execution logs.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

This pipeline was developed to streamline taxonomic FASTA file processing for downstream bioinformatics applications. Contributions, suggestions, and bug reports are welcome.

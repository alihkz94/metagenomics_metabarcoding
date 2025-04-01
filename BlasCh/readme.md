# BlasCh

![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)
![Python](https://img.shields.io/badge/python-3.6%2B-blue.svg)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## FALSE POSITIVE Chimera detection and recovery for metabarcoding & environmental DNA

BlasCh is a specialized tool designed to process BLAST XML results for identifying, classifying, and recovering false positive chimeric sequences from metabarcoding or environmental DNA (eDNA) datasets.

### Key Features

- **Automated Chimera Detection**: Identify and classify chimeric sequences
- **Sequence Recovery**: Rescue borderline sequences that might be mistakenly flagged as chimeras
- **Multiple Classification Categories**: Non-chimeric, Chimeric, Borderline, and Multiple Alignment
- **Comprehensive Reporting**: Detailed reports with classification statistics

## Installation

```bash
# Install from PyPI
pip install blasch

# OR install from GitHub
git clone https://github.com/yourusername/blasch.git
cd blasch
pip install -e .

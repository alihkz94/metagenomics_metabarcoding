configfile: "config.yaml"

# Define the base names of the FASTA files (without extension)
BASE_FILES = ["ITS", "LSU", "SSU", "longread"]
VERSION = config["version"]

rule all:
    input:
        expand("general/General_EUK_{base}_v{version}.fasta", base=BASE_FILES, version=VERSION),
        expand("mothur/mothur_EUK_{base}_v{version}.fasta", base=BASE_FILES, version=VERSION),
        expand("mothur/mothur_EUK_{base}_v{version}.tax", base=BASE_FILES, version=VERSION),
        expand("qiime2/QIIME2_EUK_{base}_v{version}.fasta", base=BASE_FILES, version=VERSION),
        expand("qiime2/QIIME2_EUK_{base}_v{version}.tsv", base=BASE_FILES, version=VERSION),
        expand("sintax/SINTAX_EUK_{base}_v{version}.fasta", base=BASE_FILES, version=VERSION),
        expand("dada2/DADA2_EUK_{base}_v{version}.fasta", base=BASE_FILES, version=VERSION)

rule general:
    input:
        "{base}.fasta"
    output:
        "general/General_EUK_{base}_v{version}.fasta"
    shell:
        "python scripts/general.py --input {input} --output {output}"

rule mothur:
    input:
        "{base}.fasta"
    output:
        fasta="mothur/mothur_EUK_{base}_v{version}.fasta",
        tax="mothur/mothur_EUK_{base}_v{version}.tax"
    shell:
        "python scripts/mothur.py --input {input} --output_fasta {output.fasta} --output_tax {output.tax}"

rule qiime:
    input:
        "{base}.fasta"
    output:
        fasta="qiime2/QIIME2_EUK_{base}_v{version}.fasta",
        tsv="qiime2/QIIME2_EUK_{base}_v{version}.tsv"
    shell:
        "python scripts/qiime.py --input {input} --output_fasta {output.fasta} --output_tsv {output.tsv}"

rule sintax:
    input:
        "general/General_EUK_{base}_v{version}.fasta"
    output:
        "sintax/SINTAX_EUK_{base}_v{version}.fasta"
    shell:
        "python scripts/sintax.py --input {input} --output {output}"

rule dada2:
    input:
        "{base}.fasta"
    output:
        "dada2/DADA2_EUK_{base}_v{version}.fasta"
    shell:
        "python scripts/dada2.py --input {input} --output {output}"

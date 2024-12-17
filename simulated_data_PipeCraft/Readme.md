## Parameter Documentation

This section details the parameters used for simulating sequencing data.

**Parameters:**

*   **`num_ref_seqs`**: Number of unique reference sequences to simulate.
    *   *Default*: 50
    *   *Impact*: A higher value increases sequence diversity for more realistic testing.

*   **`seq_length`**: Length of each reference sequence.
    *   *Default*: 250
    *   *Impact*: Matches target sequencing technology (e.g., 250 bp reads for Illumina platforms).

*   **`num_samples`**: Total number of samples to generate.
    *   *Default*: 20
    *   *Impact*: Produces multiple paired-end FASTQ files for testing.

*   **`reads_per_sample`**: Number of reads per sample.
    *   *Default*: 10000
    *   *Impact*: Defines sequencing depth, influencing downstream error rate estimation.

*   **`base_probs`**: Probabilities of nucleotide composition (A, C, G, T).
    *   *Default*: `c(A=0.25, C=0.25, G=0.25, T=0.25)`
    *   *Impact*: Default: Equal distribution of all bases.

*   **`mutation_rate`**: Rate of random mutations introduced in sequences.
    *   *Default*: 0.005
    *   *Impact*: Adds small variability while retaining overall sequence structure.

*   **`mean_abundance`**: Average abundance of sequences within a sample.
    *   *Default*: 1
    *   *Impact*: Used to model realistic variation in sequence coverage.

*   **`abundance_sd`**: Standard deviation for sequence abundances.
    *   *Default*: 0.5
    *   *Impact*: Introduces variability in the relative abundances of sequences.

*   **`qual_mean_start`**: Average PHRED quality score at the start of reads.
    *   *Default*: 35
    *   *Impact*: High scores simulate high-quality Illumina sequencing.

*   **`qual_mean_end`**: Average PHRED quality score at the end of reads.
    *   *Default*: 25
    *   *Impact*: Simulates real-world quality degradation in sequencing runs.

*   **`qual_sd`**: Variability of quality scores along each read.
    *   *Default*: 3
    *   *Impact*: Adds slight randomness around mean quality scores.

*   **`output_dir`**: Directory where the generated FASTQ files will be saved.
    *   *Default*: `"simulated_fastq_test"`
    *   *Impact*: Ensure this directory exists and the script has write permissions.

**Note:** You can modify these parameters within the script to customize your simulations.
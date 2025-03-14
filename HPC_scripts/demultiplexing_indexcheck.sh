#!/bin/bash

# Check if file was provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <fasta_file>"
    exit 1
fi

# File paths
FASTA_FILE=$1
HEADER_COUNTS="header_counts.txt"
INDEX_REPORT="index_report.txt"

# Function to check headers
check_headers() {
    # Get header counts
    echo "Checking headers..."
    grep '^>' "$FASTA_FILE" | sort | uniq -c | awk '{printf("%s\t%s\n", $1, $2)}' > "$HEADER_COUNTS"
    
    # Find duplicates
    echo "\nDuplicate headers found:"
    grep -v '^ *1\t' "$HEADER_COUNTS"
}

# Function to check indexes
check_indexes() {
    # Extract sequences and count repetitive patterns
    echo "Checking index sequences..."
    awk '
    BEGIN {RS=">"; ORS=""}
    NR>1 {
        seq=$0
        # Look for repeated patterns of length 8-20 characters
        for(len=8; len<=20; len++) {
            for(i=1; i+len<=length(seq); i++) {
                pattern=substr(seq,i,len)
                count=gsub(pattern, "&", seq)
                if(count > 1 && length(pattern) >= 8) {
                    print "Index pattern found:", pattern
                    print "Count:", count
                    print "Sequence:", seq
                    print "---"
                }
            }
        }
    }
    ' "$FASTA_FILE" > "$INDEX_REPORT"
}

# Main script
echo "Starting validation of $FASTA_FILE..."

# Check headers
check_headers

# Check indexes
check_indexes

echo "\nValidation complete!"
echo "Results saved in:"
echo "- $HEADER_COUNTS (header analysis)"
echo "- $INDEX_REPORT (index pattern analysis)"

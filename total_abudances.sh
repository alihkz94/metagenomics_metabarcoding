grep -o 'size=[0-9]*' your_file.fasta | awk -F'=' '{sum += $2} END {print sum}'

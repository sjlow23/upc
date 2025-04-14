#!/bin/bash

# Directory containing the fasta files
fasta_dir="$1"
output_dir="$2"
keyfile="$3"

# Check if the directory exists
if [ ! -d "$fasta_dir" ]; then
    echo "Directory $fasta_dir not found!"
    exit 1
fi

# Initialize the segment files by clearing or creating them
for i in {1..8}; do
    > "$output_dir"/"segment${i}.fna"  # Empty each segment file
done

# Loop over each fasta file in the directory
for input_file in "$fasta_dir"/*.fna; do
    # Ensure the file exists
    if [ ! -f "$input_file" ]; then
        echo "File $input_file not found!"
        continue
    fi

    echo "Processing file: $input_file"

    # Use awk to process each fasta file and split it into 8 parts
    awk '
    BEGIN { 
        seq_count = 0;
    }
    /^>/ {
        # When a header line is found, increment sequence count
        seq_count++;
        # Only keep first field of header
        $0 = $1;
    }
    # Print each sequence to the corresponding segment file
    {
        if (seq_count >= 1 && seq_count <= 8) {
            print $0 >> "'"$output_dir"'"/"segment" seq_count ".fna";  # Append to the correct segment file
        }
    }
    ' "$input_file"
done


for i in "$output_dir"/segment?.fna; do seqkit replace -p '(.+)' -r '{kv}' -k "$keyfile" $i > $i.renamed; \
    mv $i.renamed $i;
done
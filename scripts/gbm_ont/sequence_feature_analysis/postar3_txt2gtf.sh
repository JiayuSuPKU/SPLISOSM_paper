#!/bin/bash

# This script converts the RBS data downloaded from CLIPdb POSTAR3 in TXT format to GTF format.
# POSTAR3: http://111.198.139.65/index.html

# Input and output file
POSTAR3_DIR="/Users/jysumac/reference/POSTAR3"
INPUT_FILE="$POSTAR3_DIR/human.txt"
RBP_NAME_LIST="$POSTAR3_DIR/human.rbp_with_clip.txt"
OUTPUT_FILE="$POSTAR3_DIR/human.gtf"

# Check if input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: Input file '$INPUT_FILE' not found!"
  exit 1
fi

# sort unique RBP names in a separate file
awk -F"\t" '
BEGIN {
  # Define the GTF field separator as tab
  OFS="\t"
}
{
  # Skip empty lines or lines starting with a comment
  if (NF == 0 || $0 ~ /^#/) next;

  # Extract columns from the input file
  rbp = $6;            # RNA-binding protein

  # Print the RBP name
  print rbp;
}' "$INPUT_FILE" | sort | uniq > "$RBP_NAME_LIST"

# Create or overwrite the output file
echo "Converting $INPUT_FILE to GTF format..."
true > "$OUTPUT_FILE"

# Process each line of the input file
awk -F"\t" '
BEGIN {
  # Define the GTF field separator as tab
  OFS="\t"
}
{
  # Skip empty lines or lines starting with a comment
  if (NF == 0 || $0 ~ /^#/) next;

  # Extract columns from the input file
  seqname = $1;        # Chromosome
  start = $2;          # Start position
  end = $3;            # End position
  feature_id = $4;     # Feature ID
  strand = $5;         # Strand
  rbp = $6;            # RNA-binding protein
  source = $7;         # Source (e.g., method)
  tissue = $8;         # Tissue type
  dataset = $9;       # Dataset IDs
  score = $10;         # Confidence score

  # Construct the attributes field in GTF format
  attributes = "gene_id \"" rbp "\"; transcript_id \"" feature_id "\"; tissue \"" tissue "\"; dataset \"" dataset "\";"

  # From 0-based to 1-based coordinates
  start += 1;
  end += 1;

  # Print the GTF-formatted line
  print seqname, source, "binding_site", start, end, score, strand, ".", attributes;
}' "$INPUT_FILE" >> "$OUTPUT_FILE"

echo "Conversion complete. GTF file saved as $OUTPUT_FILE"

# sort, compress and index the GTF file
bedtools sort -i "$OUTPUT_FILE" | bgzip > "$(dirname "$OUTPUT_FILE")/$(basename "$OUTPUT_FILE" .gtf).sorted.gtf.gz"
tabix -p gff "$(dirname "$OUTPUT_FILE")/$(basename "$OUTPUT_FILE" .gtf).sorted.gtf.gz"


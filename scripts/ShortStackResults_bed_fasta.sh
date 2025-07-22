#!/bin/bash

# Usage: bash prep_bed_fasta.sh <input_gff3> <genome_fa> <output_prefix>
# Example: bash prep_bed_fasta.sh Results.gff3 genome.fa output_dir/Results

set -e

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <input_gff3> <genome_fasta> <output_prefix>"
  exit 1
fi

GFF3="$1"
GENOME="$2"
OUT_PREFIX="$3"

BED="${OUT_PREFIX}.bed"
FA="${OUT_PREFIX}.fa"

# Generate BED from GFF3
awk 'BEGIN{OFS="\t"} 
  $1 ~ /^Chr/ && $4 ~ /^[0-9]+$/ {
    print $1, $4 - 1, $5, $9, $8, $7
}' "$GFF3" > "$BED"

# Extract FASTA from genome using BED
bedtools getfasta \
  -fi "$GENOME" \
  -bed "$BED" \
  -s -name \
  -fo "$FA"

echo "[INFO] BED and FASTA files created:"
echo "BED: $BED"
echo "FA : $FA"

#!/bin/bash

# Script to analyze enrichment of known CIS-BP motifs using SEA from the MEME suite
# Usage: ./peak_enrichment.sh <peak_file> <output_dir>
# Example: ./peak_enrichment.sh /data_vol/analysis_outputs/peak_clustering/K_means_16_L2_permissive_peaks_L2_cell_types.bed /data_vol/analysis_outputs/peak_clustering

# Check if correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <peak_file> <output_dir>"
    exit 1
fi

# Input arguments
PEAK_FILE=$1
OUTPUT_DIR=$2
PEAK_FILE_PREFIX=$(basename "$PEAK_FILE" .bed)

# Ensure MEME suite is installed and SEA is available
if ! command -v sea &> /dev/null; then
    echo "Error: SEA (from MEME suite) is not installed or not in PATH."
    exit 1
fi

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Define the CIS-BP motif database (update the path if necessary)
CIS_BP_DB="/data_vol/motif_databases/CIS-BP_2.00/Homo_sapiens.meme"

# Check if the CIS-BP database exists
if [ ! -f "$CIS_BP_DB" ]; then
    echo "Error: CIS-BP motif database not found at $CIS_BP_DB."
    exit 1
fi

fasta="/data_vol/References/hg38.fa"
genome="/data_vol/References/human.hg38.genome"

# create random genomic background sequences 
if [ ! -f "$OUTPUT_DIR/random_genomic_background_seqs.fasta" ]; then
    bedtools random -l 500 -n 10000 -seed 42 -g ${genome} > ${OUTPUT_DIR}/random_genomic_background_seqs.bed
    bedtools getfasta -fi ${fasta} -bed ${OUTPUT_DIR}/random_genomic_background_seqs.bed > ${OUTPUT_DIR}/random_genomic_background_seqs.fasta
fi

# Extract unique cluster values from the fourth column of PEAK_FILE,
clusters=$(awk '{print $4}' "$PEAK_FILE" | sort -nu)

# Loop through each unique value, create a subset BED file, and run SEA
# for cluster in $clusters; do
for cluster in $clusters; do
    awk -v val="$cluster" '$4 == val' "$PEAK_FILE" | shuf -n 10000 > ${OUTPUT_DIR}/tmp.bed
    bedtools getfasta -fi ${fasta} -bed ${OUTPUT_DIR}/tmp.bed > ${OUTPUT_DIR}/tmp.fasta
    sea --text --noseqs --n ${OUTPUT_DIR}/random_genomic_background_seqs.fasta --pvalue 0.01 --p ${OUTPUT_DIR}/tmp.fasta --m ${CIS_BP_DB} > ${OUTPUT_DIR}/${PEAK_FILE_PREFIX}_cluster_${cluster}_vs_random_genomic_background_seqs_SEA_enrichment.tsv
done
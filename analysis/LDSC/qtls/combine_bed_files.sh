#!/bin/bash

# Usage: ./combine_bed_files.sh output.bed file1.bed file2.bed ... fileN.bed

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 output.bed file1.bed [file2.bed ... fileN.bed]"
  exit 1
fi

output="$1"
shift

cat "$@" > "$output"
echo "Combined BED files saved to $output"
#!/bin/bash

atac_tsv="/data_vol/FINAL/Pseudobulks/atac/stringent_L2_peaks/L2_cell_type_pseudobulk_permissive_peaks_lcpm.tsv"
outdir="/data_vol/analysis_outputs/peak_clustering"
file_prefix="L2_permissive_peaks_L2_cell_types"
for k in 4 8 12 16 20 24 28 32;
do
    echo $k
    python cluster_peaks.py -k $k --method kmeans --atac_tsv $atac_tsv --outdir $outdir --prefix $file_prefix
done
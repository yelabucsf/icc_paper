#!/bin/bash

atac_tsv="/data_vol/FINAL/Pseudobulks/atac/stringent_L2_peaks/L2_cell_type_pseudobulk_permissive_peaks_lcpm.tsv"
outdir="/data_vol/analysis_outputs/peak_clustering"
file_prefix="L2_permissive_peaks_L2_cell_types"

# cluster
# for k in 4 8 12 16 20 24 28 32;
# do
#     echo $k
#     python cluster_peaks.py -k $k --method kmeans --atac_tsv $atac_tsv --outdir $outdir --prefix $file_prefix
# done

# plot
for method in kmeans NMF;
do
    echo $method
    for k in 16;
    do
        cluster_bed="${outdir}/clustering_results/${method}_${k}_${file_prefix}.bed"
        outfile="${outdir}/figures/${method}_${k}_${file_prefix}_vertical.png"
        python plot_heatmap.py --cluster_bed $cluster_bed --atac_tsv $atac_tsv --outfile $outfile --orientation vertical
    done
done

method="kmeans"
k="16"
cluster_order="5,11,13,2,14,12,4,6,3,8,0,10,9,15,7,1"
cluster_bed="${outdir}/clustering_results/${method}_${k}_${file_prefix}.bed"
outfile="${outdir}/figures/${method}_${k}_${file_prefix}_vertical_manual_cluster_order.png"
python plot_heatmap.py --cluster_bed $cluster_bed --atac_tsv $atac_tsv --outfile $outfile --orientation vertical --cluster_order $cluster_order
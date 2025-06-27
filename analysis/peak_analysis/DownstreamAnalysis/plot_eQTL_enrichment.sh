#!/bin/bash

outdir="/data_vol/analysis_outputs/peak_clustering"
k="16"
cluster_order="5,11,13,2,14,12,4,6,3,8,0,10,9,15,7,1"
for method in kmeans NMF;
do
    cell_type_order="B,CD4T,CD8T,gdT,NK,cDC,cM,ncM"
    enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_caQTL_enrichment_all_snps.tsv"
    out_file="${outdir}/figures/${method}_${k}_L1_caQTL_enrichment_ordered_clusters_all_snps.pdf"
    bed_file="${outdir}/clustering_results/${method}_${k}_L2_permissive_peaks_L2_cell_types.bed"
    python peak_eQTL_enrichment.py -q /data_vol/FINAL/QTL/atac/permissive/L1_manual/ -b $bed_file -o $enrichment_file -m atac
    python plot_eQTL_enrichment.py -e $enrichment_file -o $out_file -c $cell_type_order -x $cluster_order

    cell_type_order="B_Naive,B_Memory,Treg,CD4T_Naive,CD4T_Memory,CD8T_Naive,CD8T_Memory,gdT,MAIT,NK_dim,cDC2,cM,ncM"
    enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_L2_manual_caQTL_enrichment_all_snps.tsv"
    out_file="${outdir}/figures/${method}_${k}_L2_caQTL_enrichment_ordered_clusters_all_snps.pdf"
    bed_file="${outdir}/clustering_results/${method}_${k}_L2_permissive_peaks_L2_cell_types.bed"
    python peak_eQTL_enrichment.py -q /data_vol/FINAL/QTL/atac/permissive/L2_manual/ -b $bed_file -o $enrichment_file -m atac
    python plot_eQTL_enrichment.py -e $enrichment_file -o $out_file -c $cell_type_order -x $cluster_order
done
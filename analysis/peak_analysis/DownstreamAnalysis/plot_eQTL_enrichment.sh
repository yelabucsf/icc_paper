#!/bin/bash

outdir="/data_vol/analysis_outputs/peak_clustering"
k="16"
for method in kmeans NMF;
do
    cell_type_order="B,CD4T,CD8T,gdT,NK,cDC,cM,ncM"
    enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_eQTL_enrichment.tsv"
    out_file="${outdir}/figures/${method}_${k}_L1_eQTL_enrichment.pdf"
    python plot_eQTL_enrichment.py -e $enrichment_file -o $out_file -c $cell_type_order

    cell_type_order="B_Naive,B_Memory,Treg,CD4T_Naive,CD4T_Memory,CD8T_Naive,CD8T_Memory,gdT,MAIT,NK_dim,cDC2,cM,ncM"
    enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_L2_manual_eQTL_enrichment.tsv"
    out_file="${outdir}/figures/${method}_${k}_L2_eQTL_enrichment.pdf"
    python plot_eQTL_enrichment.py -e $enrichment_file -o $out_file -c $cell_type_order
done
#!/bin/bash

outdir="/data_vol/analysis_outputs/peak_clustering"
k="16"
cluster_order="5,11,14,13,2,10,9,3,12,4,8,0,6,1"
for method in kmeans NMF;
do
    bed_file="${outdir}/clustering_results/exploratory/${method}_${k}_L2_permissive_peaks_L2_cell_types.bed"

    enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_eQTL_enrichment_all_snps.tsv"
    python peak_eQTL_enrichment.py -q /data_vol/FINAL/QTL/rna/L1_manual/ -b $bed_file -o $enrichment_file -m rna
    enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_caQTL_enrichment_all_snps.tsv"
    python peak_eQTL_enrichment.py -q /data_vol/FINAL/QTL/atac/permissive/L1_manual/ -b $bed_file -o $enrichment_file -m atac

    enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_L2_manual_eQTL_enrichment_all_snps.tsv"
    python peak_eQTL_enrichment.py -q /data_vol/FINAL/QTL/rna/L2_manual/ -b $bed_file -o $enrichment_file -m rna
    enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_L2_manual_caQTL_enrichment_all_snps.tsv"
    python peak_eQTL_enrichment.py -q /data_vol/FINAL/QTL/atac/permissive/L2_manual/ -b $bed_file -o $enrichment_file -m atac


    for modality in eQTL caQTL;
    do
        cell_type_order="B,CD4T,CD8T,gdT,NK,cDC,cM,ncM"
        enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_${modality}_enrichment_all_snps.tsv"
        out_file="${outdir}/figures/${method}_${k}_L1_${modality}_enrichment_ordered_clusters_all_snps.pdf"
        python plot_eQTL_enrichment.py -e $enrichment_file -o $out_file -c $cell_type_order -x $cluster_order

        cell_type_order="Treg,CD4T_Naive,CD4T_Memory,CD8T_Naive,CD8T_Memory,gdT,MAIT,NK_dim,B_Naive,B_Memory,cDC2,cM,ncM"
        enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_L2_manual_${modality}_enrichment_all_snps.tsv"
        out_file="${outdir}/figures/${method}_${k}_L2_${modality}_enrichment_ordered_clusters_all_snps.pdf"
        python plot_eQTL_enrichment.py -e $enrichment_file -o $out_file -c $cell_type_order -x $cluster_order
    done

    cell_type_order="Treg,CD4T_Naive,CD4T_Memory,CD8T_Naive,CD8T_Memory,gdT,MAIT,NK_dim,B_Naive,B_Memory,cDC2,cM,ncM"
    eqtl_enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_L2_manual_eQTL_enrichment_all_snps.tsv"
    caqtl_enrichment_file="${outdir}/eqtl_enrichment/${method}_${k}_L2_permissive_peaks_L2_cell_types_L2_manual_caQTL_enrichment_all_snps.tsv"
    
    cluster_order="5,11,14,2,10,9,3,12,4,8,0"
    out_file="${outdir}/figures/${method}_${k}_L2_enrichment_ordered_clusters_all_snps_combined_separate_colorbars_vertical_without_lineage_shared_clusters.pdf"
    python plot_QTL_enrichment_combined.py -e $eqtl_enrichment_file -ca $caqtl_enrichment_file -o $out_file -c $cell_type_order -x $cluster_order

    cluster_order="5,11,14,2,10,9,3,12,4,8,0;13,6,1"
    out_file="${outdir}/figures/${method}_${k}_L2_enrichment_ordered_clusters_all_snps_combined_separate_colorbars_vertical_separate_lineage_shared_clusters.pdf"
    python plot_QTL_enrichment_combined.py -e $eqtl_enrichment_file -ca $caqtl_enrichment_file -o $out_file -c $cell_type_order -x $cluster_order
done
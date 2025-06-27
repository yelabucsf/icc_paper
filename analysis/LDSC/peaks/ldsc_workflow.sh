#!/bin/bash

###
# 1. Comparison of different peak calling parameters
### 

# make annots for each peak set
L2_stringent="/data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/PeakCalls/stringent_peaks/stringent_L2_peaks_spm5e250_rulen12/All_Samples.fwp.filter.non_overlapping"
L2_permissive="/data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/XHLT1-POOL-XX1-SC1/peaks"
encode_cCREs="/data_vol/References/ENCODE_cCREs/GRCh38-cCREs_non_TF"

./make_annots.sh /data_vol/analysis_outputs/ldsc ${L2_stringent} L2_stringent
./make_annots.sh /data_vol/analysis_outputs/ldsc ${L2_permissive} L2_permissive
./make_annots.sh /data_vol/analysis_outputs/ldsc ${encode_cCREs} encode_cCREs

data_dir="/data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/PeakCalls/stringent_peaks"
rule="2"
for spm in 2 0 5;
do
    file_name="${data_dir}/stringent_L2_peaks_spm${spm}e250_rule${rule}/All_Samples.fwp.filter.non_overlapping"
    prefix="L2_stringent_spm${spm}e250_rule${rule}"
    ./make_annots.sh /data_vol/analysis_outputs/ldsc ${file_name} ${prefix}
done

rule="n12"
spm="0"
file_name="${data_dir}/stringent_L2_peaks_spm${spm}e250_rule${rule}/All_Samples.fwp.filter.non_overlapping"
prefix="L2_stringent_spm${spm}e250_rule${rule}"
./make_annots.sh /data_vol/analysis_outputs/ldsc ${file_name} ${prefix}

# run LDSC
ldsc_dir="/data_vol/analysis_outputs/ldsc"
for peak_set in L2_stringent L2_permissive encode_cCREs L2_stringent_spm2e250_rule2 L2_stringent_spm2e250_rulen12 L2_stringent_spm5e250_rule2 L2_stringent_spm0e250_rule2 L2_stringent_spm0e250_rulen12;
do
    ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/${peak_set}/${peak_set}. ${peak_set}
done

###
# 2. Enrichment in cell type specific peaks
### 
ldsc_dir="/data_vol/analysis_outputs/ldsc"
peak_dir="/data_vol/analysis_outputs/peak_clustering/clustering_results"
for k in {0..15}
do
    ./make_annots.sh ${ldsc_dir} ${peak_dir}/kmeans_16_L2_permissive_peaks_L2_cell_types_cluster_${k} kmeans_16/kmeans_16_L2_permissive_peaks_L2_cell_types_cluster_${k}
done

# create annots for cell type shared vs non-shared peaks
# cluster 1 is cell type shared peaks, all other clusters are non-shared
ldsc_dir="/data_vol/analysis_outputs/ldsc"
peak_dir="/data_vol/analysis_outputs/peak_clustering/clustering_results"
cat_cmd=""
for k in {0..15}
do
    if [ "$k" -ne 1 ]
    then
        cat_cmd+=" ${peak_dir}/kmeans_16_L2_permissive_peaks_L2_cell_types_cluster_${k}.bed"
    fi
done
cat ${cat_cmd} > ${peak_dir}/kmeans_16_L2_permissive_peaks_L2_cell_types_non_shared.bed
./make_annots.sh ${ldsc_dir} ${peak_dir}/kmeans_16_L2_permissive_peaks_L2_cell_types_non_shared kmeans_16/kmeans_16_L2_permissive_peaks_L2_cell_types_non_shared

# run LDSC for all cell type specific peak clusters altogether
annot_list=""
for k in {0..15}
do
    annot_list="${annot_list}${ldsc_dir}/annots/kmeans_16/kmeans_16_L2_permissive_peaks_L2_cell_types_cluster_${k}."
    if [ $k -ne 15 ]
    then
        annot_list="${annot_list},"
    fi
done
./run_ldsc.sh ${ldsc_dir} ${annot_list} kmeans_16_L2_permissive_peaks_L2_cell_types

# run LDSC for each cell type specific peak cluster separately
for k in {0..15}
do
    ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/kmeans_16/kmeans_16_L2_permissive_peaks_L2_cell_types_cluster_${k}. kmeans_16_L2_permissive_peaks_L2_cell_types_cluster_${k}
done

# run LDSC for non-shared combined annotation
 ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/kmeans_16/kmeans_16_L2_permissive_peaks_L2_cell_types_non_shared. kmeans_16_L2_permissive_peaks_L2_cell_types_non_shared
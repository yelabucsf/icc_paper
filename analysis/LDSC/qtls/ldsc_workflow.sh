#!/bin/bash

ldsc_dir="/data_vol/analysis_outputs/ldsc"

# cell types for each level
declare -A cell_types_dict
cell_types_dict["Cohort"]="ICC"
cell_types_dict["L1_manual"]="B CD4T CD8T NK cDC cM gdT ncM"
cell_types_dict["L2_manual"]="B_Memory B_Naive CD4T_Memory CD4T_Naive CD8T_Memory CD8T_Naive MAIT NK_dim Treg cDC2 cM gdT ncM"

###
# RNA
###
QTL_dir="/data_vol/FINAL/QTL/rna"

# mamba activate single_cell
# create variant and gene bed files
for level in "${!cell_types_dict[@]}"; do
    echo "Processing level: $level"
    
    combined_bed_files_cmd="./qtls/combine_bed_files.sh ${QTL_dir}/${level}/all_cell_types_significant_genes_slop_100kb.bed"
    # Convert the space-separated string into an array and iterate
    read -ra cell_types <<< "${cell_types_dict[$level]}"
    for cell_type in "${cell_types[@]}"; do
        echo "Processing cell type: $cell_type"
        
        python qtls/create_variant_beds.py -q $QTL_dir/$level --celltype $cell_type --finemap
        python qtls/create_gene_beds.py -q $QTL_dir/$level --celltype $cell_type
        combined_bed_files_cmd+=" ${QTL_dir}/${level}/${cell_type}/${cell_type}.significant_genes_slop_100kb.bed"
    done
    ${combined_bed_files_cmd}

    python qtls/create_variant_beds_all_cell_types.py -q $QTL_dir --level $level
    python qtls/create_variant_beds_modality_specific.py --level $level
done

# mamba activate py27
# make annotations for each bed file
for level in "${!cell_types_dict[@]}"; do
    echo "Processing level: $level"

    read -ra cell_types <<< "${cell_types_dict[$level]}"
    for cell_type in "${cell_types[@]}"; do
        mkdir -p ${ldsc_dir}/annots/qtls/rna/${level}/${cell_type}

        bed_file="${QTL_dir}/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant"
        ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/rna/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant

        # bed_file="${QTL_dir}/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_CS"
        # ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/rna/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_CS

        bed_file="${QTL_dir}/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_PIP_0.3"
        ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/rna/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_PIP_0.3

        bed_file="${QTL_dir}/${level}/${cell_type}/${cell_type}.significant_genes_slop_100kb"
        ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/rna/${level}/${cell_type}/${cell_type}.significant_genes_slop_100kb
    done

    mkdir -p ${ldsc_dir}/annots/qtls/rna/${level}/all_cell_types
    bed_file="${QTL_dir}/${level}/all_cell_types_significant_genes_slop_100kb"
    ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/rna/${level}/all_cell_types/all_cell_types_significant_genes_slop_100kb

    for prefix in "significant_qtls_pooled" "significant_qtl_in_one_cell_type" "significant_qtl_in_multiple_cell_types"; do
        bed_file="${QTL_dir}/${level}/all_cell_types.${prefix}"
        ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/rna/${level}/all_cell_types/${prefix}
    done

    mkdir -p ${ldsc_dir}/annots/qtls/${level}/all_cell_types
    for prefix in "significant_qtls_pooled_rna_and_atac" "significant_qtls_pooled_rna_only" "significant_qtls_pooled_atac_only"; do
        bed_file="${QTL_dir}/${level}/all_cell_types.${prefix}"
        ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/${level}/all_cell_types/${prefix}
    done
done

###
# ATAC
###
QTL_dir="/data_vol/FINAL/QTL/atac/permissive"

# mamba activate single_cell
# create variant and gene bed files
for level in "${!cell_types_dict[@]}"; do
    echo "Processing level: $level"
    
    combined_bed_files_cmd="./qtls/combine_bed_files.sh ${QTL_dir}/${level}/all_cell_types_significant_peaks_1000bp.bed"
    # Convert the space-separated string into an array and iterate
    read -ra cell_types <<< "${cell_types_dict[$level]}"
    for cell_type in "${cell_types[@]}"; do
        echo "Processing cell type: $cell_type"
        
        python qtls/create_variant_beds.py -q $QTL_dir/$level --celltype $cell_type --finemap
        python qtls/create_peak_beds.py -q $QTL_dir/$level --celltype $cell_type
        combined_bed_files_cmd+=" ${QTL_dir}/${level}/${cell_type}/${cell_type}.significant_peaks_1000bp.bed"
    done
    ${combined_bed_files_cmd}
    python qtls/create_variant_beds_all_cell_types.py -q $QTL_dir --level $level
done

# mamba activate py27
# make annotations for each bed file
for level in "${!cell_types_dict[@]}"; do
    echo "Processing level: $level"

    read -ra cell_types <<< "${cell_types_dict[$level]}"
    for cell_type in "${cell_types[@]}"; do
        mkdir -p ${ldsc_dir}/annots/qtls/atac/${level}/${cell_type}

        bed_file="${QTL_dir}/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant"
        ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/atac/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant

        # bed_file="${QTL_dir}/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_CS"
        # ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/atac/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_CS

        bed_file="${QTL_dir}/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_PIP_0.3"
        ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/atac/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_PIP_0.3

        bed_file="${QTL_dir}/${level}/${cell_type}/${cell_type}.significant_peaks_1000bp"
        ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/atac/${level}/${cell_type}/${cell_type}.significant_peaks_1000bp
    done

    mkdir -p ${ldsc_dir}/annots/qtls/atac/${level}/all_cell_types
    bed_file="${QTL_dir}/${level}/all_cell_types_significant_peaks_1000bp"
    ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/atac/${level}/all_cell_types/all_cell_types_significant_peaks_1000bp

    for prefix in "significant_qtls_pooled" "significant_qtl_in_one_cell_type" "significant_qtl_in_multiple_cell_types"; do
        bed_file="${QTL_dir}/${level}/all_cell_types.${prefix}"
        ./make_annots.sh ${ldsc_dir} ${bed_file} qtls/atac/${level}/all_cell_types/${prefix}
    done
done

###
# Run LDSC
###
for modality in rna atac; do
    echo "Processing modality: $modality" 

    for level in "${!cell_types_dict[@]}"; do
        echo "Processing level: $level" 

        read -ra cell_types <<< "${cell_types_dict[$level]}"
        for cell_type in "${cell_types[@]}"; do
            echo "Processing cell type: $cell_type"

            ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/qtls/${modality}/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant. ${modality}_${cell_type}.cis_qtl_pairs.significant
            # ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/qtls/${modality}/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_CS. ${modality}_${cell_type}.cis_qtl_pairs.significant.SuSiE_CS
            # ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/qtls/${modality}/${level}/${cell_type}/${cell_type}.cis_qtl_pairs.significant.SuSiE_PIP_0.3. ${modality}_${cell_type}.cis_qtl_pairs.significant.SuSiE_PIP_0.3

            # ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/qtls/${modality}/${level}/${cell_type}/${cell_type}.significant_genes_slop_100kb. ${modality}_${cell_type}.significant_genes_slop_100kb
            # ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/qtls/${modality}/${level}/${cell_type}/${cell_type}.significant_peaks_1000bp. ${modality}_${cell_type}.significant_peaks_1000bp
        done
        ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/qtls/${modality}/${level}/all_cell_types/all_cell_types_significant_genes_slop_100kb. ${modality}_all_cell_types.significant_genes_slop_100kb
        ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/qtls/${modality}/${level}/all_cell_types/all_cell_types_significant_peaks_1000bp. ${modality}_all_cell_types.significant_peaks_1000bp

        for prefix in "significant_qtls_pooled" "significant_qtl_in_one_cell_type" "significant_qtl_in_multiple_cell_types" "significant_qtls_pooled_rna_and_atac" "significant_qtls_pooled_rna_only" "significant_qtls_pooled_atac_only"; do
            ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/qtls/${modality}/${level}/all_cell_types/${prefix}. ${modality}_${level}_all_cell_types.${prefix}
        done
    done    
done

for level in "${!cell_types_dict[@]}"; do
    echo "Processing level: $level" 

    for prefix in "significant_qtls_pooled_rna_and_atac" "significant_qtls_pooled_rna_only" "significant_qtls_pooled_atac_only"; do
        ./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/qtls/${level}/all_cell_types/${prefix}. ${level}_all_cell_types.${prefix}
    done
done 


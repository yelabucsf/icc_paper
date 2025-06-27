#!/bin/bash

# cell types for each level
declare -A cell_types_dict
cell_types_dict["Cohort"]="ICC"
cell_types_dict["L1_manual"]="B CD4T CD8T NK cDC cM gdT ncM"
cell_types_dict["L2_manual"]="B_Memory B_Naive CD4T_Memory CD4T_Naive CD8T_Memory CD8T_Naive MAIT NK_dim Treg cDC2 cM gdT ncM"

out_dir="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi"

for level in "${!cell_types_dict[@]}"; do
    echo "Processing level: $level"

    read -ra cell_types <<< "${cell_types_dict[$level]}"

    # eQTLs
    parallel -j 4 ./icc_eqtl_process/01_LD_clump_per_eGene_hg38.sh $level {} ::: "${cell_types[@]}"
    Rscript icc_qtl_props/01_filter_eqtls.R $level

    # caQTLs
    parallel -j 4 ./icc_caqtl_process/01_LD_clump_per_caPeak_hg38.sh $level {} ::: "${cell_types[@]}"
    Rscript icc_qtl_props/01_filter_caqtls.R $level
done

for modality in rna atac; do
    for level in "${!cell_types_dict[@]}"; do
        parallel -j 12 Rscript icc_qtl_props/02_bootstrap_hits.R {} $modality $level ::: {1..1000}
        Rscript icc_qtl_props/03_match_SNPs.R $modality $level

        # basic properties
        Rscript icc_qtl_props/04_calc_properties.R $modality/$level/filter_qtls.hg19.assoc $modality/$level/out/qtl.basic_props
        parallel -j 12 Rscript icc_qtl_props/04_calc_properties.R $modality/$level/boot/boot_{}.txt $modality/$level/out/qtl_boot.{}.basic_props ::: {1..1000}
        parallel -j 12 Rscript icc_qtl_props/04_calc_properties.R $modality/$level/match/qtl.match_snps_{}.txt $modality/$level/out/qtl_matched.{}.basic_props ::: {1..1000}
        # combine files for bootstrap/matched samples
        cat $out_dir/$modality/$level/out/qtl_boot.1.basic_props >  $out_dir/$modality/$level/out/qtl_boot.basic_props
        cat $out_dir/$modality/$level/out/qtl_matched.1.basic_props >  $out_dir/$modality/$level/out/qtl_matched.basic_props
        for i in {2..1000}; do
            tail -n 1 $out_dir/$modality/$level/out/qtl_boot.$i.basic_props >> $out_dir/$modality/$level/out/qtl_boot.basic_props
            tail -n 1 $out_dir/$modality/$level/out/qtl_matched.$i.basic_props >> $out_dir/$modality/$level/out/qtl_matched.basic_props
        done
        Rscript icc_qtl_props/04_extract_props_mean_CI.R $modality $level

        # TSS clustering
        Rscript icc_qtl_props/05_calc_TSS_clustering.R $modality/$level/filter_qtls.hg19.assoc $modality/$level/out/qtl.TSS_props
        parallel -j 12 Rscript icc_qtl_props/05_calc_TSS_clustering.R $modality/$level/boot/boot_{}.txt $modality/$level/out/qtl_boot.{}.TSS_props ::: {1..1000}
        parallel -j 12 Rscript icc_qtl_props/05_calc_TSS_clustering.R $modality/$level/match/qtl.match_snps_{}.txt $modality/$level/out/qtl_matched.{}.TSS_props ::: {1..1000}
        # combine files for bootstrap/matched samples
        cat $out_dir/$modality/$level/out/qtl_boot.1.TSS_props >  $out_dir/$modality/$level/out/qtl_boot.TSS_props
        cat $out_dir/$modality/$level/out/qtl_matched.1.TSS_props >  $out_dir/$modality/$level/out/qtl_matched.TSS_props
        for i in {2..1000}; do
            tail -n 1 $out_dir/$modality/$level/out/qtl_boot.$i.TSS_props >> $out_dir/$modality/$level/out/qtl_boot.TSS_props
            tail -n 1 $out_dir/$modality/$level/out/qtl_matched.$i.TSS_props >> $out_dir/$modality/$level/out/qtl_matched.TSS_props
        done
        Rscript icc_qtl_props/05_extract_TSS_clustering_mean_CI.R $modality $level

        # promoter/enhancer properties
        Rscript icc_qtl_props/06_calc_promoter_enhancer_props.R $modality/$level/filter_qtls.hg19.assoc $modality/$level/out/qtl.promoter_enhancer.enrichments
        parallel -j 12 Rscript icc_qtl_props/06_calc_promoter_enhancer_props.R $modality/$level/boot/boot_{}.txt $modality/$level/out/qtl_boot.{}.promoter_enhancer.enrichments ::: {1..1000}
        parallel -j 12 Rscript icc_qtl_props/06_calc_promoter_enhancer_props.R $modality/$level/match/qtl.match_snps_{}.txt $modality/$level/out/qtl_matched.{}.promoter_enhancer.enrichments ::: {1..1000}
        # combine files for bootstrap/matched samples
        cat $out_dir/$modality/$level/out/qtl_boot.1.promoter_enhancer.enrichments >  $out_dir/$modality/$level/out/qtl_boot.promoter_enhancer.enrichments
        cat $out_dir/$modality/$level/out/qtl_matched.1.promoter_enhancer.enrichments >  $out_dir/$modality/$level/out/qtl_matched.promoter_enhancer.enrichments
        for i in {2..1000}; do
            tail -n 1 $out_dir/$modality/$level/out/qtl_boot.$i.promoter_enhancer.enrichments >> $out_dir/$modality/$level/out/qtl_boot.promoter_enhancer.enrichments
            tail -n 1 $out_dir/$modality/$level/out/qtl_matched.$i.promoter_enhancer.enrichments >> $out_dir/$modality/$level/out/qtl_matched.promoter_enhancer.enrichments
        done
        Rscript icc_qtl_props/06_extract_promoter_enhancer_enrichments_CI.R $modality $level
    done
done
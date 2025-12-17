#!/bin/bash


out_dir="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi"

# GWAS
parallel -j 20 Rscript gwas_props/02_bootstrap_hits.R {} ::: {1..1000}
Rscript gwas_props/03_match_SNPs.R

# basic properties
Rscript gwas_props/04_calc_properties.R gwas/filter_indep_gwas.assoc gwas/out/gwas.basic_props
parallel -j 20 Rscript gwas_props/04_calc_properties.R gwas/boot/boot_{}.txt gwas/out/gwas_boot.{}.basic_props ::: {1..1000}
parallel -j 10 Rscript gwas_props/04_calc_properties.R gwas/match/gwas.match_snps_{}.txt gwas/out/gwas_matched.{}.basic_props ::: {1..1000}
# combine files for bootstrap/matched samples
cat $out_dir/gwas/out/gwas_boot.1.basic_props >  $out_dir/gwas/out/gwas_boot.basic_props
cat $out_dir/gwas/out/gwas_matched.1.basic_props >  $out_dir/gwas/out/gwas_matched.basic_props
for i in {2..1000}; do
    tail -n 1 $out_dir/gwas/out/gwas_boot.$i.basic_props >> $out_dir/gwas/out/gwas_boot.basic_props
    tail -n 1 $out_dir/gwas/out/gwas_matched.$i.basic_props >> $out_dir/gwas/out/gwas_matched.basic_props
done
Rscript gwas_props/04_extract_props_mean_CI.R

# TSS clustering
Rscript gwas_props/05_calc_TSS_clustering.R gwas/filter_indep_gwas.assoc gwas/out/gwas.TSS_props
parallel -j 20 Rscript gwas_props/05_calc_TSS_clustering.R gwas/boot/boot_{}.txt gwas/out/gwas_boot.{}.TSS_props ::: {1..1000}
parallel -j 10 Rscript gwas_props/05_calc_TSS_clustering.R gwas/match/gwas.match_snps_{}.txt gwas/out/gwas_matched.{}.TSS_props ::: {1..1000}
# combine files for bootstrap/matched samples
cat $out_dir/gwas/out/gwas_boot.1.TSS_props >  $out_dir/gwas/out/gwas_boot.TSS_props
cat $out_dir/gwas/out/gwas_matched.1.TSS_props >  $out_dir/gwas/out/gwas_matched.TSS_props
for i in {2..1000}; do
    tail -n 1 $out_dir/gwas/out/gwas_boot.$i.TSS_props >> $out_dir/gwas/out/gwas_boot.TSS_props
    tail -n 1 $out_dir/gwas/out/gwas_matched.$i.TSS_props >> $out_dir/gwas/out/gwas_matched.TSS_props
done
Rscript gwas_props/05_extract_TSS_clustering_mean_CI.R

# promoter/enhancer properties
Rscript gwas_props/06_calc_promoter_enhancer_props.R gwas/filter_indep_gwas.assoc gwas/out/gwas.promoter_enhancer.enrichments
parallel -j 20 Rscript gwas_props/06_calc_promoter_enhancer_props.R gwas/boot/boot_{}.txt gwas/out/gwas_boot.{}.promoter_enhancer.enrichments ::: {1..1000}
parallel -j 10 Rscript gwas_props/06_calc_promoter_enhancer_props.R gwas/match/gwas.match_snps_{}.txt gwas/out/gwas_matched.{}.promoter_enhancer.enrichments ::: {1..1000}
# combine files for bootstrap/matched samples
cat $out_dir/gwas/out/gwas_boot.1.promoter_enhancer.enrichments >  $out_dir/gwas/out/gwas_boot.promoter_enhancer.enrichments
cat $out_dir/gwas/out/gwas_matched.1.promoter_enhancer.enrichments >  $out_dir/gwas/out/gwas_matched.promoter_enhancer.enrichments
for i in {2..1000}; do
    tail -n 1 $out_dir/gwas/out/gwas_boot.$i.promoter_enhancer.enrichments >> $out_dir/gwas/out/gwas_boot.promoter_enhancer.enrichments
    tail -n 1 $out_dir/gwas/out/gwas_matched.$i.promoter_enhancer.enrichments >> $out_dir/gwas/out/gwas_matched.promoter_enhancer.enrichments
done
Rscript gwas_props/06_extract_promoter_enhancer_enrichments_CI.R
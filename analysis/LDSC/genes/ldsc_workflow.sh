#!/bin/bash

# create bed file for 100kb window around expressed genes
expressed_genes_bed="/data_vol/FINAL/QTL/rna/Cohort/ICC/pheno/ICC_phenos_tensorqtl_sorted.bed"
out_dir="/data_vol/analysis_outputs/ldsc/annots"
mkdir ${out_dir}/expressed_genes
cut -f1-4 ${expressed_genes_bed} | tail -n +2 | sed 's/^/chr/' | bedtools slop -b 100000 -g /data_vol/References/human.hg38.genome > ${out_dir}/expressed_genes/ICC_all_expressed_genes_slop_100kb.bed

# make annots for expressed genes
mamba activate py27
ldsc_dir="/data_vol/analysis_outputs/ldsc"
bed_file="${ldsc_dir}/annots/expressed_genes/ICC_all_expressed_genes_slop_100kb"
./make_annots.sh ${ldsc_dir} ${bed_file} expressed_genes/ICC_all_expressed_genes_slop_100kb

ldsc_dir="/data_vol/analysis_outputs/ldsc"
./run_ldsc.sh ${ldsc_dir} ${ldsc_dir}/annots/expressed_genes/ICC_all_expressed_genes_slop_100kb. ICC_all_expressed_genes_slop_100kb
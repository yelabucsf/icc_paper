
### Run after pseudobulking pipeline


### Steps for prepping final QTL files from pseudobulks and covariate files

## Run pheno PCs for each cell type
data_dir="/data/processed/pseudobulk/rna"
out_dir="/data/processed/QTL_calling/rna/"
genos_dir="/data/processed/genotypes/plink/"
geno_pref='combined_cleaned'
cell_annot_col='L2_manual'
pref='L3'
#chr=5

conda activate tensorQTL
## Generate Expression PCS
for cell_type in `ls -1 ${data_dir}/${cell_annot_col}/*.exp | cut -d / -f7 | cut -d . -f1`; do python make_phenoPC.py -p 40 -e ${data_dir}/${cell_annot_col}/${cell_type}_expr.lcpm.qnorm -o ${data_dir}/${cell_annot_col}/${cell_type}_expr_lcpm_qnorm_pcs.txt ; done

## LOOK AT ELBOW PLOT HERE TO DETERMINE THE NUMBER OF PCS

PC_num=10

## merge covariates
conda activate archr
for cell_type in `ls -1 ${data_dir}/${cell_annot_col}/*.exp | cut -d / -f7 | cut -d . -f1`; do mkdir -p  ${out_dir}/${cell_annot_col}/${cell_type}/cov/ ; Rscript combine_covariates.R -Z ${data_dir}/${cell_annot_col}/${cell_type}_expr_lcpm_qnorm_pcs.txt -M /data/processed/meta/person_demo_anc.2025-03-14.txt -C ${cell_type} -P ${cell_type} -B ${PC_num} -O ${out_dir}/${cell_annot_col}/${cell_type}/cov/ ; done 


## filter and order matrices for tensorQTL
conda activate archr
for cell_type in `ls -1 ${data_dir}/${cell_annot_col}/*.exp | cut -d / -f7 | cut -d . -f1`; do mkdir -p ${out_dir}/${cell_annot_col}/${cell_type}/pheno ; Rscript filter_order_matrices_tensorQTL.R -X ${data_dir}/${cell_annot_col}/${cell_type}_expr.lcpm.qnorm -V ${out_dir}/${cell_annot_col}/${cell_type}/cov/${cell_type}_${cell_type}.covs  -P ${cell_type}  -O ${out_dir}/${cell_annot_col}/${cell_type}/pheno -L /data/references/gene_locs_autosome_clean.txt ;done

## sort beds
conda activate genetics
for cell_type in `ls -1 ${data_dir}/${cell_annot_col}/*.exp | cut -d / -f7 | cut -d . -f1`; do head -n 1 ${out_dir}/${cell_annot_col}/${cell_type}/pheno/${cell_type}_phenos_tensorqtl.bed > ${out_dir}/${cell_annot_col}/${cell_type}/pheno/${cell_type}_phenos_tensorqtl_sorted.bed; tail -n +2 ${out_dir}/${cell_annot_col}/${cell_type}/pheno/${cell_type}_phenos_tensorqtl.bed | sort -k1,1V -k2,2n -k3,3n >> ${out_dir}/${cell_annot_col}/${cell_type}/pheno/${cell_type}_phenos_tensorqtl_sorted.bed ; done

#bedtools sort -i -g ${out_dir}/${pref}_phenos_tensorqtl.bed >> ${out_dir}/${pref}_phenos_tensorqtl_sorted.bed


## run tensor QTL
conda activate tensorQTL
for cell_type in `ls -1 ${data_dir}/${cell_annot_col}/*.exp | cut -d / -f7 | cut -d . -f1`; do python run_tensorqtl.py -g ${genos_dir}/${geno_pref} -p ${out_dir}/${cell_annot_col}/${cell_type}/pheno/${cell_type}_phenos_tensorqtl_sorted.bed -c ${out_dir}/${cell_annot_col}/${cell_type}/cov/${cell_type}_${cell_type}.covs -n ${out_dir}/${cell_annot_col}/${cell_type}/${cell_type} ; done


### run tensor QTL
#conda activate tensorQTL
#for cell_type in `ls -1 ${data_dir}/*.exp | cut -d / -f6 | cut -d . -f1`; do python run_tensorqtl.py -g ${genos_dir}/chr${chr}.dose -p ${out_dir}/${cell_type}/pheno/${pref}_phenos_tensorqtl_sorted.bed -c ${out_dir}/${cell_type}/cov/${pref}_${cell_type}.covs -n ${out_dir}/${cell_type}/${pref} --chrom ${chr} ; done
#

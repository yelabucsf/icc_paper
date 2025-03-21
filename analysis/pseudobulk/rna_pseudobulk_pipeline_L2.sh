

cell_annot_col='L2_manual'
h5_file='/data/march2025/MultiVI_posttraining_RNA_Celltypist_VaryingClusteringResolution_manual_annots.h5ad'
outdir='/data/processed/pseudobulk/rna/'
donor_col='best_sample'
cell_num=50

# make pseudobulk counts
conda activate scanpy

mkdir -p ${outdir}/${cell_annot_col}
python make_basic_counts.py -R ${h5_file} -O ${outdir}/${cell_annot_col}/ -S ${donor_col} ${cell_annot_col} -G ${cell_annot_col}  


# lcpm normalize
conda activate archr
for cell_type in `ls -1 ${outdir}${cell_annot_col}/*.exp | cut -d / -f7 | cut -d . -f1`;
do echo $cell_type;  Rscript normalized_lcpm.R -E ${outdir}/${cell_annot_col}/${cell_type}.exp -G ${outdir}/${cell_annot_col}/gene_counts_exp.txt -D ${outdir}/${cell_annot_col}/donor_counts_exp.txt -C ${cell_type} -O  ${outdir}/${cell_annot_col} -N ${cell_num} -F ${outdir}/${cell_annot_col}/donor_cell_counts_exp.txt  ; done





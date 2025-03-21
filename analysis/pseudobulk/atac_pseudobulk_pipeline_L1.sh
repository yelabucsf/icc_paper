
# Run stringent peak pipeline first
# takes in peak h5ad file
# We are using stringent peaks from L2 for L1 QTLs since these are almost identical 

cell_annot_col='L1_manual'
#h5_file='/data/processed/pseudobulk/atac/stringent_L1_peaks/aggregate_phase2_ICC_stringent_L1_peaks_counts_mtx.h5'
h5_file='/data/processed/pseudobulk/atac/stringent_L2_peaks/aggregate_phase2_ICC_stringent_L2_peaks_counts_mtx_annot.h5'
outdir='/data/processed/pseudobulk/atac/'
donor_col='best_sample'
cell_num=50

# make pseudobulk counts
conda activate scanpy

# add metadata to atac h5 file, note that this is hardcoded and should only be run if annotations have not already been added
python add_meta_atac_h5_L1.py
python add_meta_atac_h5_L2.py

# make basic counts
mkdir -p ${outdir}/${cell_annot_col}
python make_basic_counts.py -A ${h5_file} -O ${outdir}/${cell_annot_col}/ -C ${donor_col} ${cell_annot_col} -G ${cell_annot_col}  


# normalize
conda activate archr
for cell_type in `ls -1 ${outdir}${cell_annot_col}/*.atac | cut -d / -f7 | cut -d . -f1`;
do echo $cell_type;  Rscript normalized_lcpm.R -E ${outdir}/${cell_annot_col}/${cell_type}.atac -G ${outdir}/${cell_annot_col}/gene_counts_atac.txt -D ${outdir}/${cell_annot_col}/donor_counts_atac.txt -C ${cell_type} -O  ${outdir}/${cell_annot_col} -N ${cell_num} -F ${outdir}/${cell_annot_col}/donor_cell_counts_atac.txt -M atac ; done





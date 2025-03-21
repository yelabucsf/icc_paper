import scanpy as sc
import pandas as pd

adata=sc.read_h5ad('/data/processed/pseudobulk/atac/stringent_L2_peaks/aggregate_phase2_ICC_stringent_L2_peaks_counts_mtx.h5')
meta=pd.read_csv('/data/march2025/MultiVI_posttraining_RNA_Celltypist_VaryingClusteringResolution_manual_annots_filtered_meta.csv')
print(meta.head())
meta.index=meta['Unnamed: 0']
meta.index.name = ''
print(meta.head())
tmp=adata.obs.merge(meta, left_index=True, right_index=True, how='left')
#tmp=adata.obs.merge(meta, by.x=0, by.y=0, how='left')
adata.obs=tmp
adata.write('/data/processed/pseudobulk/atac/stringent_L2_peaks/aggregate_phase2_ICC_stringent_L2_peaks_counts_mtx_annot.h5')

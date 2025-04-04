### create a peak location file (like gen location file) from h5ad


#### make gene loc files for peaks
import scanpy as sc
import pandas as pd

# read in stringent peak file 
adata=sc.read_h5ad('/data/processed/pseudobulk/atac/permissive_peaks/aggregate_phase2_ICC_permissive_L2_peaks_counts_mtx_annot.h5')

peaks=adata.var.index.tolist()
chr=[ i.split(':')[0] for i in peaks ]
locs=[ i.split(':')[1] for i in peaks ]
s1=[ i.split('-')[0] for i in locs ]
s2=[ i.split('-')[1] for i in locs ]

dict={'gene': peaks, 'chr': chr, 's1':s1,'s2':s2}

df=pd.DataFrame(dict)
df.to_csv('/data/references/peak_locs_L2_permissive.txt', sep='\t', index=False)


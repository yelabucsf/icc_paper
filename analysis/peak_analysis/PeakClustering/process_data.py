import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import qnorm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--atac_h5",
                      action = "store", type=str, required=False,
                      dest = "atac_h5", 
                      default = "/data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/aggregate_phase2_ICC_permissive_L2_peaks_counts_mtx.h5",
                      help = "ATAC h5 file") 
    parser.add_argument("-o", "--outdir",
                      action = "store", type=str, required=True,
                      dest = "outdir",
                      help = "Output directory to save files to.") 
    parser.add_argument("-p", "--prefix",
                      action = "store", type=str, required=True,
                      dest = "file_prefix",
                      help = "File prefix to use too save files.") 
    args = parser.parse_args()

    # load obs from RNA h5
    adata_rna = sc.read_h5ad("/data_vol/FINAL/SingleCellObjects/rna/MultiVI_posttraining_RNA_Celltypist_VaryingClusteringResolution_manual_annots.h5ad",
                            backed="r")
    rna_obs = adata_rna.obs
    rna_obs.index = rna_obs["Well"].astype(str) + "#" + rna_obs["Barcode"].astype(str)
    del adata_rna

    # load atac h5
    adata = sc.read_h5ad(args.atac_h5)
    adata.obs = adata.obs.join(rna_obs, how="left")

    norm_factors = np.array(adata.X.sum(axis=0)).flatten()
    norm_factors = (10^6)/norm_factors
    cell_types = np.unique(adata.obs["L2_manual"])

    cellbypeaks_unnorm  = pd.DataFrame([], index=adata.var.index)
    for cell_type in tqdm(cell_types):
        cellbypeaks_unnorm[cell_type] = np.array(adata[adata.obs["L2_manual"] == cell_type].X.sum(axis=0)).flatten()
    cellbypeaks_unnorm.to_csv(f"{args.outdir}/{args.file_prefix}_cell_type_pseudobulk_unnorm.tsv",
                            sep="\t", header=True, index=True)

    # log cpm and quantile normalize
    cellbypeaks_logcpm = np.log2((cellbypeaks_unnorm.T*norm_factors).T + 1)
    cellbypeaks_logcpm = qnorm.quantile_normalize(cellbypeaks_logcpm)  
    cellbypeaks_L2_logcpm = qnorm.quantile_normalize(cellbypeaks_L2_logcpm)  
    cellbypeaks_logcpm.to_csv(f"{args.outdir}/{args.file_prefix}_cell_type_pseudobulk_lcpm.tsv",
                            sep="\t", header=True, index=True)


if __name__ == '__main__':
    main()
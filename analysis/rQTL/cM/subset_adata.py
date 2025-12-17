#!/usr/bin/env python3
# subset anndata RNA and ATAC files to cM + IL1B cells

import anndata as ad
import pandas as pd
from pathlib import Path
from subprocess import run

dir_dat = Path("./dat/processed/cm_il1b/reqtl/")
dir_dat.mkdir(parents=True, exist_ok=True)

dir_sc = Path("./dat/processed/sc_objects/")
dir_sc_rna = dir_sc / "rna"
dir_sc_atac = dir_sc / "atac"


s3_rna = "s3://gracie-icc/production_data/analysis_outputs/FINAL/SingleCellObjects/rna/MultiVI_posttraining_RNA_Celltypist_VaryingClusteringResolution_manual_annots_filtered.h5ad"
rna_flt_fn = dir_sc_rna / "rna_flt.h5ad"
ret = run(
    [
        "aws", "s3", "cp",
        s3_rna, rna_flt_fn
    ]
)

s3_atac = "s3://gracie-icc/production_data/analysis_outputs/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/aggregate_phase2_ICC_permissive_L2_peaks_counts_mtx.h5"

atac_fn = dir_sc_atac / "atac_flt.h5ad"
ret = run(
    [
        "aws", "s3", "cp",
        s3_atac, atac_fn
    ]
)


rna_flt = ad.read_h5ad(rna_flt_fn, backed='r')
atac = ad.read_h5ad(atac_fn, backed='r')

atac_bcs = atac.obs_names.str.replace("#", ";").to_series()
atac_bcs.index = atac_bcs.values
atac.obs_names = atac_bcs.values

keep_cts = ["cM", "cM_IL1B"]

k_bcs = rna_flt.obs_names[rna_flt.obs["L3_manual"].isin(keep_cts)]

rna_flt_s = rna_flt[k_bcs].to_memory()

flt_meta = rna_flt_s.obs.copy()
flt_meta_fn = dir_dat / "meta.tsv"
flt_meta.to_csv(
    flt_meta_fn, header=True, index=True, sep='\t'
)

obs_k = [
    "Well", "Pool", "best_sample", "rna_pct_mt",
    "L0_manual", "L1_manual", "L2_manual", "L3_manual"
]

rna_flt_s.obs = rna_flt_s.obs[obs_k]

uns_keys = list(rna_flt_s.uns.keys())
for k in uns_keys:
    del rna_flt_s.uns[k]

obsm_keys = list(rna_flt_s.obsm.keys())
for k in obsm_keys:
    del rna_flt_s.obsm[k]

obsp_keys = list(rna_flt_s.obsp.keys())
for k in obsp_keys:
    del rna_flt_s.obsp[k]

rna_flt_s.X = rna_flt_s.layers["counts"]
del rna_flt_s.layers["counts"]

rna_out_fn = dir_dat / "rna.h5ad"
rna_flt_s.write_h5ad(rna_out_fn)

k_bcs_atac = k_bcs[k_bcs.isin(atac_bcs)]
atac_s = atac[k_bcs_atac].to_memory()

atac_s.obs = rna_flt_s.obs.loc[atac_s.obs_names]

atac_out_fn = dir_dat / "atac.h5ad"
atac_s.write_h5ad(atac_out_fn)

# write sample IDs present
rna_flt_sam = rna_flt.obs["best_sample"].unique()
sam_out_fn = dir_dat / "sample_ids.txt"
rna_flt_sam = pd.Series(rna_flt_sam)
rna_flt_sam.to_csv(sam_out_fn, header=False, index=False)


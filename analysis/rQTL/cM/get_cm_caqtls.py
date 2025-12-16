#!/usr/bin/env python3
# Get caQTLs in cM

import pandas as pd
import numpy as np
from pathlib import Path
from subprocess import run

dir_dat = Path("./dat/processed/cm_il1b/reqtl/caqtls/")
dir_dat.mkdir(exist_ok=True, parents=True)

caqtl_fn = Path("./dat/processed/qtl/atac/permissive/L2_manual/")

cell_type = "cM"

# pull lead caQTL-caPeak pairs from S3
s3_dir = f"s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/L2_manual/{cell_type}"

caqtl_l = []
for chr in range(1, 23):
    s3_atac_chr = f"{s3_dir}/{cell_type}.cis_qtl_egene.{chr}.txt"
    df = pd.read_csv(s3_atac_chr, sep='\t', header=0)
    df["cell_type"] = cell_type
    caqtl_l.append(df)
caqtl_df = pd.concat(caqtl_l)

sk = caqtl_df["qval"] < .05
caqtl_df = caqtl_df.loc[sk]

pairs = caqtl_df[["phenotype_id", "variant_id"]]
pair_fn = dir_dat / "cM_caqtls.txt"
pairs.to_csv(pair_fn, header=True, index=False, sep="\t")

egenes = caqtl_df["phenotype_id"].unique()
esnps = caqtl_df["variant_id"].unique()

egenes_fn = dir_dat / "cM_caPeaks.txt"
np.savetxt(egenes_fn, egenes, fmt='%s')

esnps_fn = dir_dat / "cM_casnps.txt"
np.savetxt(esnps_fn, esnps, fmt='%s')

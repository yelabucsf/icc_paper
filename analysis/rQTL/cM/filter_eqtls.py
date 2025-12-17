#!/usr/bin/env python3
# Get candidate SNPs for testing

import pandas as pd
import numpy as np
from pathlib import Path

dir_dat = Path("./dat/processed/cm_il1b/reqtl/eqtls/")
dir_dat.mkdir(exist_ok=True, parents=True)

eqtl_fn = Path("./dat/processed/qtl/rna/L2_manual/all.cis_eqtl.tsv.gz")
eqtl = pd.read_csv(eqtl_fn, header=0, sep='\t')

k = eqtl["cell_type"] == "cM"
eqtl_k = eqtl.loc[k]

sk = eqtl_k["qval"] < .05
eqtl_k = eqtl_k.loc[sk]

pairs = eqtl_k[["phenotype_id", "variant_id"]]
pair_fn = dir_dat / "cM_eqtls.txt"
pairs.to_csv(pair_fn, header=True, index=False, sep="\t")

egenes = eqtl_k["phenotype_id"].unique()
esnps = eqtl_k["variant_id"].unique()

egenes_fn = dir_dat / "cM_egenes.txt"
np.savetxt(egenes_fn, egenes, fmt='%s')

esnps_fn = dir_dat / "cM_esnps.txt"
np.savetxt(esnps_fn, esnps, fmt='%s')

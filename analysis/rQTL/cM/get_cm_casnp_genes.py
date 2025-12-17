#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pyranges as pr
import anndata as ad

WINDOW = 100_000

# read gene locations
genes_fn = "./dat/processed/cm_il1b/reqtl/gene_tss_loc.bed"
genes_pr = pr.read_bed(genes_fn)

# subset to those in RNA data
rna_fn = "./dat/processed/cm_il1b/reqtl/rna.h5ad"
rna = ad.read_h5ad(rna_fn, backed='r')

rna_vars = rna.var_names
k = genes_pr.Name.isin(rna_vars)
genes_pr = genes_pr[k].copy()

rna.file.close()
del rna

# read SNP locations
snps_fn = "./dat/processed/cm_il1b/reqtl/vcfs/chrALL.info.bed"
snps_pr = pr.read_bed(snps_fn)
# remove duplicates, usually multiallelic
k = ~snps_pr.Name.duplicated()
snps_pr = snps_pr[k]

# read in caSNPs
casnps_fn = "./dat/processed/cm_il1b/reqtl/caqtls/cM_casnps.txt"
casnps = np.loadtxt(casnps_fn, dtype=str)

# subset snps_pr to caSNPs
k = snps_pr.Name.isin(pd.Series(casnps))
ca_snps_pr = snps_pr[k]

# get all genes
ca_snps_win_pr = ca_snps_pr.extend(ext=WINDOW)
ca_snps_win_gene_pr = ca_snps_win_pr.join(genes_pr)

# get distance between snps and genes
snp_loc = pd.DataFrame(
    {
        "Chromosome": ca_snps_pr.Chromosome.values,
        "Position": ca_snps_pr.End.values,
    },
    index=ca_snps_pr.Name
)
dist = ca_snps_win_gene_pr.End_b.values - snp_loc.loc[ca_snps_win_gene_pr.Name, "Position"].values

# join pairs
ca_snp_gene_pairs = pd.DataFrame(
    {
        "variant_id": ca_snps_win_gene_pr.Name,
        "gene_id": ca_snps_win_gene_pr.Name_b,
        "distance": dist
    }
)

snp_gene_pairs_fn = "./dat/processed/cm_il1b/reqtl/caqtls/cM_casnp_gene_pairs.txt"
ca_snp_gene_pairs.to_csv(snp_gene_pairs_fn, header=True, index=False, sep='\t')

#!/usr/bin/env python3

from tqdm import tqdm
from cyvcf2 import VCF
import pandas as pd
import numpy as np
import gpboost as gpb
import anndata as ad
from pathlib import Path
from scipy.stats import chi2
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

# ==============================================================================
# centering fx for covaraites
# ==============================================================================


def center_and_scale(x):
    """
    Center and scale a 1-D numpy array.
    Returns:
        x_scaled : transformed array with mean 0 and sd 1
        mean     : original mean of x
        sd       : original standard deviation of x
    """
    x = np.asarray(x, dtype=float)
    
    mean = np.nanmean(x)
    sd = np.nanstd(x, ddof=0)
    
    if sd == 0:
        raise ValueError("SD of variable is zero; cannot scale.")
    
    x_scaled = (x - mean) / sd
    return x_scaled, mean, sd


# ==============================================================================
# read snps, genes, samples, and vcf
# ==============================================================================

# read caQTL-gene pairs
pairs_fn = "./dat/processed/cm_il1b/reqtl/caqtls/cM_casnp_gene_pairs.txt"
pairs = pd.read_csv(pairs_fn, header=0, sep='\t')

# read SNP locations
snp_loc_fn = "./dat/processed/cm_il1b/reqtl/vcfs/chrALL.info.vcf"
snp_loc = pd.read_csv(snp_loc_fn, header=None, sep='\t')
# deduplicate, use first occurrence for pos
k = ~snp_loc[2].duplicated()
snp_loc = snp_loc.loc[k]
snp_loc.index = snp_loc[2].values

# get sample IDs
samples_fn = Path("./dat/processed/cm_il1b/reqtl/sample_ids.txt")
sample_ids = np.loadtxt(samples_fn, dtype=str)

# vcf file
vcf_fn = "./dat/processed/cm_il1b/reqtl/vcfs/chrALL.dose.flt.bcf"
vcf = VCF(vcf_fn, samples=sample_ids.tolist(), gts012=True)
vcf_sams = vcf.samples

# ==============================================================================
# Parse RNA data and metadata/covariates
# ==============================================================================

# RNA adata
rna_fn = Path("./dat/processed/cm_il1b/reqtl/rna.h5ad")
rna = ad.read_h5ad(rna_fn)

# PC1 condition
pc1_fn = Path("./dat/processed/cm_il1b/reqtl/pc1.txt")
pc1 = pd.read_csv(pc1_fn, index_col=0, sep='\t', header=None)

meta_fn = Path("./dat/processed/cm_il1b/reqtl/meta.tsv")
meta = pd.read_csv(meta_fn, sep='\t', header=0, index_col=0)

# map sex
meta["sex_code"] = 0
meta.loc[meta["Sex"] == "Female", "sex_code"] = 1
# map site
meta["site_cumc"] = 0
meta.loc[meta["Site"] == "CUMC", "site_cumc"] = 1
# map genotype array
meta["geno_gda"] = 0
meta.loc[meta["GenotypingArray"] == "GDA", "geno_gda"] = 1
meta["geno_other"] = 0
meta.loc[meta["GenotypingArray"] == "Other", "geno_other"] = 1
# map geno samples to meta samples
geno_sam_idx_map = {str(s): i for i, s in enumerate(sample_ids)}
meta_geno_sam = np.array(
    [geno_sam_idx_map[a] for a in meta['best_sample']]
)
meta['geno_sam_ix'] = meta_geno_sam

# subset RNA to cMs
rna = rna[meta.index.values].to_memory()

Xcsc = rna.X.tocsc()

log_umis = np.log(rna.X.sum(-1).A[:, 0])
meta["log_umi"] = log_umis

meta["expr_PC1"] = pc1.loc[meta.index]

donor = meta["best_sample"].astype("category").cat.codes.values
batch = meta["Well"].astype("category").cat.codes.values
group = np.vstack([donor, batch]).T

# rna gene to int index
gene2ix = {
    s: i for i, s in enumerate(rna.var_names)
}

# build base covariates
covs_dict = {
    "intercept": np.ones(meta.shape[0]),
    "log_umi": meta["log_umi"].values,
    "age": center_and_scale(meta["Age"].values)[0],
    "sex": meta["sex_code"].values,
    "rna_pct_mt": center_and_scale(meta["rna_pct_mt"].values)[0],
    "site": meta["site_cumc"].values,
    "geno_gda": meta["geno_gda"].values,
    "geno_other": meta["geno_other"].values,
    "context": center_and_scale(meta["expr_PC1"].values)[0],
}
for i in range(1, 4):
    covs_dict[f"geno_pc{i}"] = center_and_scale(meta[f"PC{i}"].values)[0]


# ==============================================================================
# Run interaction QTL analysis
# ==============================================================================

snps_uniq = pairs["variant_id"].unique()
res_l = []

for snp in tqdm(snps_uniq):
    snp_k = pairs["variant_id"] == snp
    snp_pairs = pairs.loc[snp_k]
    var_chr = snp_loc.loc[snp, 0]
    var_pos0 = int(snp_loc.loc[snp, 1] - 1)
    var_pos = int(var_pos0 + 1)

    hits = []
    for v in vcf(f"{var_chr}:{var_pos0}-{var_pos}"):
        if v.ID == snp:
            hits.append(v)

    if not hits:
        print(f"No genotypes found for {snp}")
        continue
    if len(hits) > 1:
        print(f"Multiple hits found for {snp}")

    hit = hits[0]

    if len(hit.ALT) > 1:
        print(f"snp {snp} is multiallelic ({len(hit.ALT)})")

    alt_counts = hit.gt_types.copy().astype(np.float32)
    alt_counts[alt_counts == 3] = np.nan  # set missing to nan

    # ensure samples and genotyeps are correctly matched
    alt_counts_a = pd.Series(
        alt_counts,
        index=vcf_sams
    )
    alt_counts = alt_counts_a.loc[sample_ids].values
    cell_geno = alt_counts_a.loc[meta["best_sample"]].values

    ref_allele = hit.REF
    alt_alleles = ",".join(hit.ALT)

    alt_af = np.nanmean(alt_counts) / 2.0
    alt_var = np.nanvar(cell_geno)
    if np.isnan(alt_var) or alt_var <= 0.0:
        print(f"Variant is monomorphic")
        continue

    # get genotypes
    for gene_id in snp_pairs['gene_id']:
        gene_ix = gene2ix.get(gene_id)
        if gene_ix is None:
            print(f"Gene '{gene_id}' not found in anndata")
            continue
        gene_expr = Xcsc[:, gene_ix].toarray()[:, 0]

        covs_dict["geno"] = center_and_scale(cell_geno)[0]
        covs_dict["geno*context"] = covs_dict["context"] * covs_dict["geno"]

        cell_valid = np.isfinite(cell_geno)

        # fit full model
        covs = pd.DataFrame(covs_dict, index=meta.index)

        gp_full = gpb.GPModel(group_data=donor[cell_valid], likelihood="poisson")
        gp_full.fit(
            y=gene_expr[cell_valid], X=covs.values[cell_valid]
        )

        # fit reduced model
        covs_reduced = covs.iloc[:, :-1]
        gp_reduced = gpb.GPModel(group_data=donor[cell_valid], likelihood="poisson")
        gp_reduced.fit(
            y=gene_expr[cell_valid], X=covs_reduced.values[cell_valid]
        )

        ll_full = -gp_full.get_current_neg_log_likelihood()
        ll_reduced = -gp_reduced.get_current_neg_log_likelihood()

        chisq = 2 * (ll_full - ll_reduced)
        df = 1
        lrt_p = chi2.sf(chisq, df)

        coefs = gp_full.get_coef(std_err=True, format_pandas=False)
        coef = coefs[0, -1]
        se = coefs[1, -1]
        z = coef / se
        wald_p = 2 * norm.sf(abs(z))

        g_coef, g_se = coefs[0, -2], coefs[1, -2]
        g_z = g_coef / g_se
        g_wald_p = 2 * norm.sf(abs(g_z))

        vg_k = (snp_pairs["variant_id"] == snp) & (snp_pairs["gene_id"] == gene_id)
        var_gene_dist = int(snp_pairs.loc[vg_k]["distance"].iloc[0])
        res = pd.DataFrame(
            {
                "gene": gene_id,
                "variant": snp,
                "distance": var_gene_dist,
                "ref": ref_allele,
                "alt": alt_alleles,
                "alt_af": alt_af,
                "geno_coef": float(g_coef),
                "geno_se": float(g_se),
                "geno_wald_p": float(g_wald_p),
                "interaction_coef": float(coef),
                "interaction_wald_p": float(wald_p),
                "interaction_lrt_p": float(lrt_p),
            },
            index=pd.Index([f"{gene_id}:{snp}"])
        )

        res_l.append(res)

if len(res_l) == 0:
    raise RuntimeError("No interactions QTLs run")

rna.file.close()
vcf.close()

res_df = pd.concat(res_l)

reject, qvals, _, _ = multipletests(
    res_df["interaction_lrt_p"].values, alpha=0.05, method='fdr_bh'
)
res_df["interaction_lrt_q"] = qvals

i_wald_p = res_df["interaction_wald_p"].values
kna = ~np.isnan(i_wald_p)
reject, qvals, _, _ = multipletests(
    i_wald_p[kna], alpha=0.05, method='fdr_bh'
)
res_df["interaction_wald_q"] = np.nan
res_df.loc[kna, "interaction_wald_q"] = qvals

res_df = res_df.sort_values("interaction_wald_p")

# ==============================================================================
# Save
# ==============================================================================

res_df_fn = Path("./dat/processed/cm_il1b/reqtl/cm.caqtl.pc1_reqtl.txt")
res_df.to_csv(res_df_fn, header=True, index=False, sep='\t')


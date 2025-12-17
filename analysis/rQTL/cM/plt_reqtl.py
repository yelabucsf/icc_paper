#!/usr/bin/env python3
# plot some reQTL results

from cyvcf2 import VCF
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import matplotlib.pyplot as plt


def get_geno(vcf, var_chr, var_pos0, var_pos, var_id) -> pd.Series | None:
    hits = []
    for v in vcf(f"{var_chr}:{var_pos0}-{var_pos}"):
        if v.ID == var_id:
            hits.append(v)

    if not hits:
        print(f"No genotypes found for {var_id}")
        return None
    if len(hits) > 1:
        print(f"Multiple hits found for {var_id}")

    hit = hits[0]

    if len(hit.ALT) > 1:
        print(f"snp {var_id} is multiallelic ({len(hit.ALT)})")

    alt_counts = hit.gt_types.copy().astype(np.float32)
    alt_counts[alt_counts == 3] = np.nan  # set missing to nan

    # ensure samples and genotyeps are correctly matched
    alt_counts_a = pd.Series(
        alt_counts,
        index=vcf.samples
    )

    return alt_counts_a


dir_plt = Path("./res/cm_il1b/reqtl/cm.eqtl.pc1_reqtl/")
dir_plt.mkdir(exist_ok=True, parents=True)


rna_fn = Path("./dat/processed/cm_il1b/reqtl/rna.h5ad")
rna = ad.read_h5ad(rna_fn)

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

pairs_fn = Path("./dat/processed/cm_il1b/reqtl/eqtls/cM_eqtls.txt")
pairs = pd.read_csv(pairs_fn, header=0, sep='\t')

# PC1 condition
pc1_fn = Path("./dat/processed/cm_il1b/reqtl/pc1.txt")
pc1 = pd.read_csv(pc1_fn, index_col=0, sep='\t', header=None)

meta_fn = Path("./dat/processed/cm_il1b/reqtl/meta.tsv")
meta = pd.read_csv(meta_fn, sep='\t', header=0, index_col=0)

# map geno samples to meta samples
geno_sam_idx_map = {str(s): i for i, s in enumerate(sample_ids)}
meta_geno_sam = np.array(
    [geno_sam_idx_map[a] for a in meta['best_sample']]
)
meta['geno_sam_ix'] = meta_geno_sam

rna = rna[meta.index.values].to_memory()

# normalize values for plotting
sc.pp.normalize_total(rna, target_sum=1e3)
Xcsc = rna.X.tocsc()

meta["expr_PC1"] = pc1.loc[meta.index]

donor = meta["best_sample"].astype("category").cat.codes.values

# rna gene to index
gene2ix = {
    s: i for i, s in enumerate(rna.var_names)
}

res_df_fn = Path("./dat/processed/cm_il1b/reqtl/cm.eqtl.pc1_reqtl.txt")
res_df = pd.read_csv(res_df_fn, header=0, sep='\t')

res_df = res_df.sort_values("interaction_wald_p")


n_plt = min(50, res_df.shape[0])
for plt_ix in range(n_plt):
    gene_id, var_id = res_df["gene"].iloc[plt_ix], res_df["variant"].iloc[plt_ix]
    gene_name = rna.var.loc[gene_id, "name"]
    ibeta = res_df["interaction_coef"].iloc[plt_ix]
    ip = res_df["interaction_wald_q"].iloc[plt_ix]

    var_chr = snp_loc.loc[var_id, 0]
    var_pos0 = int(snp_loc.loc[var_id, 1] - 1)
    var_pos = int(var_pos0 + 1)

    alt_counts_a = get_geno(
        vcf, var_chr, var_pos0, var_pos, var_id
    )
    if alt_counts_a is None:
        continue

    cell_geno = alt_counts_a.loc[meta["best_sample"]].values

    gene_ix = gene2ix.get(gene_id)
    if gene_ix is None:
        print(f"Gene {gene_id} not found")
        continue

    gene_expr = Xcsc[:, gene_ix].toarray()[:, 0]
    context = meta["expr_PC1"].values

    quartiles = pd.qcut(context, q=4, labels=False, duplicates="drop") + 1

    plt_df = pd.DataFrame(
        {
            "geno": cell_geno,
            "gene_expr": gene_expr,
            "quartile": quartiles,
            "donor": donor,
        }
    )

    geno_lvls = [0, 1, 2]

    pdf_fn = dir_plt / f"{gene_id}.{var_id}.iqtl.pdf"
    png_fn = dir_plt / f"{gene_id}.{var_id}.iqtl.png"

    fig, axes = plt.subplots(1, 4, figsize=(8, 4), sharey=True)
    rng = np.random.default_rng(0)

    for q, ax in zip(range(1, 5), axes):
        df_q = plt_df.loc[plt_df["quartile"] == q]
        # donor-level means within each genotype for this quartile
        df_q = df_q.assign(log_expr=np.log1p(df_q["gene_expr"]))
        donor_mean = (
            df_q.groupby(["donor", "geno"], observed=True)["log_expr"]
            .mean()
            .reset_index()
        )
        group_vals = [donor_mean.loc[donor_mean["geno"] == g, "log_expr"] for g in geno_lvls]

        ax.boxplot(group_vals, positions=geno_lvls, widths=0.6, showfliers=False)
        # jittered donor means per genotype
        for g, vals in zip(geno_lvls, group_vals):
            if len(vals) == 0:
                continue
            x_jitter = g + rng.normal(0, 0.08, size=len(vals))
            ax.scatter(x_jitter, vals, s=6, alpha=0.2, color="black", edgecolors="none", zorder=2)

        # add linear trend across genotype means
        x = np.array(geno_lvls)
        y = np.array([d.mean() if len(d) else np.nan for d in group_vals])
        ok = ~np.isnan(y)
        if ok.sum() >= 2:
            coef = np.polyfit(x[ok], y[ok], 1)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax.plot(x_line, np.polyval(coef, x_line), color="blue", lw=2, zorder=3)

        ax.set_title(f"Quartile {q}")
        ax.set_xticks(geno_lvls)
        ax.set_xlabel("geno")

    fig.suptitle(f"{gene_name} - {var_id} beta={ibeta:.2e} q={ip:.2e}")
    axes[0].set_ylabel("log1p(gene_expr)")
    plt.tight_layout()
    plt.savefig(pdf_fn)
    plt.savefig(png_fn)
    plt.close(fig)

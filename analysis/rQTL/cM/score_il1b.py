#!/usr/bin/env python3
# Characterize cM IL1B cells

import anndata as ad
import random
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm, ListedColormap
from matplotlib.ticker import PercentFormatter
import numpy as np
from pathlib import Path
import boto3
import json
import scanpy as sc
from pandas.api import types as pdtypes


def create_umap_plot(
    data,
    *,
    color_by=None,
    x_col=0,
    y_col=1,
    figsize=(6, 6),
    point_size=6,
    alpha=0.8,
    rasterized=True,
    cmap="viridis",
    palette=None,
    discrete=None,
    legend=True,
    legend_title=None,
    legend_kwargs=None,
    colorbar=False,
    colorbar_kwargs=None,
    axis_labels=("UMAP 1", "UMAP 2"),
    title=None,
    ax=None,
    scatter_kwargs=None,
):
    if isinstance(x_col, int):
        x = data.iloc[:, x_col]
    else:
        x = data[x_col]

    if isinstance(y_col, int):
        y = data.iloc[:, y_col]
    else:
        y = data[y_col]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    plot_kwargs = {
        "s": point_size,
        "linewidths": 0,
        "alpha": alpha,
        "rasterized": rasterized,
    }
    if scatter_kwargs:
        plot_kwargs.update(scatter_kwargs)

    color_series = None
    if color_by is not None:
        if isinstance(color_by, str):
            color_series = data[color_by]
        elif isinstance(color_by, pd.Series):
            color_series = color_by.reindex(data.index)
        else:
            color_series = pd.Series(color_by, index=data.index)

        if len(color_series) != len(data):
            raise ValueError("color_by must have the same length as the input data")

    if color_series is None:
        ax.scatter(x, y, **plot_kwargs)
    else:
        if discrete is None:
            discrete = not pdtypes.is_numeric_dtype(color_series)

        if discrete:
            if palette is not None:
                lookup = dict(palette)
                category_order = [cat for cat in palette if cat in color_series.values]
            else:
                lookup = {}
                category_order = []

            seen = set()
            for value in color_series:
                if pd.isna(value):
                    continue
                if value in seen:
                    continue
                seen.add(value)
                if palette is None:
                    category_order.append(value)

            color_lookup = lookup.copy()
            if palette is None:
                base_cmap = plt.cm.tab20
                for idx, category in enumerate(category_order):
                    color_lookup.setdefault(category, base_cmap(idx % base_cmap.N))

            mapped_colors = color_series.map(color_lookup)
            mapped_colors = mapped_colors.fillna("lightgray")
            ax.scatter(x, y, c=mapped_colors, **plot_kwargs)

            if legend and category_order:
                handles = [
                    plt.Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        markerfacecolor=color_lookup.get(category, "lightgray"),
                        markersize=point_size,
                        label=str(category),
                    )
                    for category in category_order
                ]
                legend_options = {
                    "title": legend_title,
                    "bbox_to_anchor": (1.05, 1),
                    "loc": "upper left",
                    "fontsize": "small",
                }
                if legend_kwargs:
                    legend_options.update(legend_kwargs)
                ax.legend(handles=handles, **legend_options)
        else:
            scatter = ax.scatter(x, y, c=color_series, cmap=cmap, **plot_kwargs)
            if colorbar:
                cb_kwargs = dict(colorbar_kwargs or {})
                if "label" not in cb_kwargs:
                    cb_kwargs["label"] = legend_title or (color_by if isinstance(color_by, str) else None)
                fig.colorbar(scatter, ax=ax, **{k: v for k, v in cb_kwargs.items() if v is not None})

    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    ax.set_title(title or "")
    ax.grid(False)
    return fig, ax


dir_plt = Path("./res/cm_il1b/reqtl/")
dir_plt.mkdir(exist_ok=True, parents=True)

# Set random seeds for reproducibility
random.seed(42)
np.random.seed(42)


adat_fn = Path("./dat/processed/cm_il1b/reqtl/rna.h5ad")
adat = ad.read_h5ad(adat_fn)

sc.pp.normalize_total(adat)
sc.pp.log1p(adat)

# gene id map
ens_fn = "ref/cellranger/ens_name.txt"
ens_map = pd.read_csv(ens_fn, sep='\t', header=None)
gene2ens = pd.Series(
    ens_map[0].values,
    index=ens_map[1].values
)

ens2gene = pd.Series(
    ens_map[1].values,
    index=ens_map[0].values
)

# read in DE results
de_fn = "./res/de_da/2025-03-27/cm_il1b_rna_de.csv"
rna_cm_il1b_de = pd.read_csv(de_fn, sep=',')

# score genes
k = (
    (rna_cm_il1b_de["pval_adj"] < 0.05) & (rna_cm_il1b_de["log2fc"] > 0)
)
il1b_genes = rna_cm_il1b_de.loc[k, "gene_id"]

il1b_genes = rna_cm_il1b_de["gene_id"].iloc[:50]

sc.tl.score_genes(
    adat, il1b_genes.to_list(), score_name="il1b_score",
    ctrl_as_ref=False
)

# get IL1B expression
kg = adat.var["name"] == "IL1B"
ens_il1b = adat.var_names.values[kg][0]
kvar = np.where(adat.var_names == ens_il1b)[0]
il1b_expr = adat.X[:, kvar].todense().A[:, 0]
adat.obs["il1b_expr"] = il1b_expr

# ==============================================================================
# Run PCA on RNA subset
# ==============================================================================

sc.pp.highly_variable_genes(adat, n_top_genes=2000)
sc.tl.pca(adat)
sc.pp.neighbors(adat)
sc.tl.umap(adat, key_added="X_umap_sub")

# ==============================================================================
# Plot all PCs
# ==============================================================================

s_umap_df = pd.DataFrame(adat.obsm["X_umap_sub"], index=adat.obs_names)
s_umap_df["cell_type"] = adat.obs["L3_manual"].astype(str)
s_umap_df["il1b_score"] = adat.obs["il1b_score"]
s_umap_df["il1b_expr"] = adat.obs["il1b_expr"]
for i in range(10):
    s_umap_df[f"PC{i+1}"] = adat.obsm["X_pca"][:, i]

cell_type_vals = adat.obs["L3_manual"].astype(str)

# scatter plot of UMAP colored by cell type and IL1B score
fig, axs = plt.subplots(6, 2, figsize=(12, 24))
cell_types = sorted(cell_type_vals.unique())
color_map = {ct: plt.cm.tab20(i % plt.cm.tab20.N) for i, ct in enumerate(cell_types)}

axsf = axs.flatten()

create_umap_plot(
    s_umap_df,
    color_by="cell_type",
    palette=color_map,
    legend_title="cell type",
    title="RNA (cM) UMAP — cell type",
    ax=axsf[0],
    scatter_kwargs={"s": 6},
)

create_umap_plot(
    s_umap_df,
    color_by="il1b_expr",
    colorbar=True,
    colorbar_kwargs={"label": "IL1B\nexpression"},
    title="RNA (cM) UMAP — IL1B expression",
    ax=axsf[1],
    scatter_kwargs={"s": 6},
)

for i in range(10):
    create_umap_plot(
        s_umap_df,
        color_by=f"PC{i+1}",
        colorbar=True,
        colorbar_kwargs={"label": "PC"},
        title=f"RNA (cM) UMAP — PC{i+1}",
        ax=axsf[i+2],
        scatter_kwargs={"s": 6},
    )

fig.tight_layout()
fig.savefig(dir_plt / "rna_cM_umap_il1b_scatter.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ==============================================================================
# Plot PC1 on RNA
# ==============================================================================

s_umap_df = pd.DataFrame(adat.obsm["X_umap_sub"], index=adat.obs_names)
s_umap_df["cell_type"] = adat.obs["L3_manual"].astype(str)
s_umap_df["il1b_score"] = adat.obs["il1b_score"]
s_umap_df["il1b_expr"] = adat.obs["il1b_expr"]
for i in range(10):
    s_umap_df[f"PC{i+1}"] = adat.obsm["X_pca"][:, i]

cell_type_vals = adat.obs["L3_manual"].astype(str)

# scatter plot of UMAP colored by cell type and IL1B score
fig, axs = plt.subplots(1, 3, figsize=(12, 3))
cell_types = sorted(cell_type_vals.unique())
color_map = {ct: plt.cm.tab20(i % plt.cm.tab20.N) for i, ct in enumerate(cell_types)}

axsf = axs.flatten()

create_umap_plot(
    s_umap_df,
    color_by="cell_type",
    palette=color_map,
    legend_title="cell type",
    title="RNA (cM) UMAP — cell type",
    ax=axsf[0],
    scatter_kwargs={"s": 6},
)

create_umap_plot(
    s_umap_df,
    color_by="il1b_expr",
    colorbar=True,
    colorbar_kwargs={"label": "IL1B\nexpression"},
    title="RNA (cM) UMAP — IL1B expression",
    ax=axsf[1],
    scatter_kwargs={"s": 6},
)

create_umap_plot(
    s_umap_df,
    color_by="PC1",
    colorbar=True,
    colorbar_kwargs={"label": "PC"},
    title=f"RNA (cM) UMAP — PC1",
    ax=axsf[2],
    scatter_kwargs={"s": 6},
)

fig.tight_layout()
fig.savefig(dir_plt / "rna_cM_umap_il1b_PC1_scatter.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ==============================================================================
# save PCs
# ==============================================================================

dir_dat = Path("./dat/processed/cm_il1b/reqtl")

pc1_dat = pd.DataFrame(
    {"PC1": adat.obsm["X_pca"][:, 0]},
    index=adat.obs_names
)

pc1_fn = dir_dat / "pc1.txt"
pc1_dat.to_csv(pc1_fn, header=False, index=True, sep='\t')

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

cluster_names = {
    5: "Classical monocyte",
    11: "Monocyte",
    13: "Myeloid",
    2: "Myeloid & B",
    14: "cDC",
    12: "Innate-like",
    4: "T",
    6: "Lymphoid",
    3: "NK",
    8: "T Naive",
    0: "Treg & proliferating",
    10: "B & Plasmablast",
    9: "B",
    15: "pDC",
    7: "Progenitor",
    1: "Shared"
}
colors = {"eQTLs": "Reds", 
          "caQTLs": "Blues"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--enrichment_results_eqtl",
                        action="store", type=str, required=True,
                        dest="enrichment_results_eqtl")
    parser.add_argument("-ca", "--enrichment_results_caqtl",
                        action="store", type=str, required=True,
                        dest="enrichment_results_caqtl")
    parser.add_argument("-c", "--cell_type_order",
                        action="store", type=str, required=False,
                        default=None, dest="cell_type_order")
    parser.add_argument("-x", "--cluster_order",
                        action="store", type=str, required=False,
                        default=None, dest="cluster_order")
    parser.add_argument("-o", "--outfile",
                        action="store", type=str, required=True,
                        dest="outfile",
                        help="Output file to save plot.")
    args = parser.parse_args()

    enrichment_dfs = {"eQTLs": pd.read_csv(args.enrichment_results_eqtl,
                                          sep="\t", index_col=0),
                      "caQTLs": pd.read_csv(args.enrichment_results_caqtl,
                                           sep="\t", index_col=0)}

    if args.cell_type_order is not None:
        cell_type_order = args.cell_type_order.split(",")
        for k, v in enrichment_dfs.items():
            enrichment_dfs[k] = v.loc[cell_type_order]
    if args.cluster_order is not None:
        cluster_order = args.cluster_order.split(";")
        cluster_order = [c.split(",") for c in cluster_order]

    enrichment_df_sigs = {}
    for k, v in enrichment_dfs.items():
        _, pvals, _, _ = multipletests(v.astype(np.float64).values.flatten(),
                                       alpha=0.01, method="fdr_bh")
        enrichment_dfs[k] = pd.DataFrame(pvals.reshape(v.shape),
                                         index=v.index,
                                         columns=v.columns)
        enrichment_df_sigs[k] = np.full(enrichment_dfs[k].shape, "").astype("<U4")
        for i, thresh in enumerate([1e-2, 1e-3, 1e-4]):
            print(thresh, "*"*(i+1))
            enrichment_df_sigs[k][np.where(enrichment_dfs[k] < thresh)] = "*"*(i+1)
        enrichment_df_sigs[k] = pd.DataFrame(enrichment_df_sigs[k],
                                             index=enrichment_dfs[k].index,
                                             columns=enrichment_dfs[k].columns)

    fig, ax = plt.subplots(2, len(cluster_order), figsize=(6, 10), gridspec_kw={'width_ratios': [2.9, 1]})
    plt.subplots_adjust(wspace=-0.05, hspace=0.15)
    x_center = (ax[0, 0].get_position().x0 + ax[0, 1].get_position().x1) / 2
    # cbar_ax = fig.add_axes([.99, .3, .03, .5])

    for i, k in enumerate(enrichment_dfs.keys()):
        for j, c in enumerate(cluster_order):
            sns.heatmap(-np.log(enrichment_dfs[k][c]), ax=ax[i, j],
                        annot=enrichment_df_sigs[k][c].values.astype(str), 
                        cmap=colors[k],
                        vmin=0, 
                        vmax=30 if i == 0 else 15, 
                        fmt='', square=True,
                        xticklabels=False if i == 0 else True,
                        yticklabels=False if j == 1 else True,
                        cbar_kws={"shrink": 0.65},
                        cbar=False if j == 0 else True)
                        # cbar_ax=None if i == 0 else cbar_ax)
            if i == 1:
                ax[i, j].set_xticklabels([cluster_names[int(c.get_text())] for c in ax[i, j].get_xticklabels()], 
                                         rotation=45, ha='right', fontsize=13)
            if j == 0:
                ax[i, j].set_ylabel("QTL cell type", fontsize=20)
                ax[i, j].set_yticklabels(ax[i, j].get_yticklabels(), rotation=0, ha='right', fontsize=14)
                ax[i, j].set_title(k, fontsize=20, x=x_center)

    ax[1, 0].set_xlabel("Open chromatin clusters", x=x_center, fontsize=20)     
    plt.savefig(args.outfile, bbox_inches="tight", dpi=300)


if __name__ == '__main__':
    main()

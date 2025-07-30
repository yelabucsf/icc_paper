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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--enrichment_results",
                        action="store", type=str, required=True,
                        dest="enrichment_results")
    parser.add_argument("-c", "--cell_type_order",
                        action="store", type=str, required=False,
                        default=None, dest="cell_type_order")
    parser.add_argument("-x", "--cluster_order",
                        action="store", type=str, required=False,
                        default=None, dest="cluster_order")
    parser.add_argument("-v", "--vmax",
                        action="store", type=int, required=False,
                        default=25, dest="vmax")
    parser.add_argument("-o", "--outfile",
                        action="store", type=str, required=True,
                        dest="outfile",
                        help="Output file to save plot.")
    args = parser.parse_args()

    enrichment_df = pd.read_csv(args.enrichment_results,
                                sep="\t", index_col=0)

    if args.cell_type_order is not None:
        cell_type_order = args.cell_type_order.split(",")
        enrichment_df = enrichment_df.loc[cell_type_order]
    if args.cluster_order is not None:
        cluster_order = args.cluster_order.split(",")
        enrichment_df = enrichment_df[cluster_order]
    enrichment_df.columns = [cluster_names[int(c)] for c in enrichment_df.columns]

    _, pvals, _, _ = multipletests(enrichment_df.astype(np.float64).values.flatten(),
                                   alpha=0.01, method="fdr_bh")
    enrichment_df_corrected = pd.DataFrame(pvals.reshape(enrichment_df.shape),
                                           index=enrichment_df.index,
                                           columns=enrichment_df.columns)
    enrichment_df_sigs = np.full(enrichment_df.shape, "").astype("<U4")
    for i, thresh in enumerate([1e-2, 1e-3, 1e-4]):
        print(thresh, "*"*(i+1))
        enrichment_df_sigs[np.where(enrichment_df_corrected < thresh)] = "*"*(i+1)

    if args.cluster_order is None:
        fg = sns.clustermap(-np.log(enrichment_df_corrected), 
                            row_cluster=False,
                            annot=enrichment_df_sigs.astype(str), cmap="rocket",
                            figsize=(7,5), vmin=0, vmax=args.vmax, fmt='')
    else:
        fg = sns.clustermap(-np.log(enrichment_df_corrected), 
                            row_cluster=False,
                            col_cluster=False, square=True,
                            annot=enrichment_df_sigs.astype(str), cmap="rocket",
                            figsize=(6,6), vmin=0, vmax=args.vmax, fmt='',
                            cbar_pos=(0.08, 0.31, 0.03, 0.48))
    fg.ax_heatmap.set_title("caQTL enrichment in open chromatin clusters", fontsize=14)
    fg.ax_heatmap.set_xlabel("Open chromatin clusters", fontsize=14)
    fg.ax_heatmap.set_ylabel("QTL cell type", fontsize=14)
    plt.savefig(args.outfile, bbox_inches="tight", dpi=300)


if __name__ == '__main__':
    main()

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--enrichment_results",
                        action="store", type=str, required=True,
                        dest="enrichment_results")
    parser.add_argument("-c", "--cell_type_order",
                        action="store", type=str, required=False,
                        default=None, dest="cell_type_order")
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

    _, pvals, _, _ = multipletests(enrichment_df.astype(np.float64).values.flatten(),
                                   alpha=0.01, method="fdr_bh")
    enrichment_df_corrected = pd.DataFrame(pvals.reshape(enrichment_df.shape),
                                           index=enrichment_df.index,
                                           columns=enrichment_df.columns)
    enrichment_df_sigs = np.full(enrichment_df.shape, "").astype("<U4")
    for i, thresh in enumerate([1e-2, 1e-3, 1e-4]):
        print(thresh, "*"*(i+1))
        enrichment_df_sigs[np.where(enrichment_df_corrected < thresh)] = "*"*(i+1)

    fg = sns.clustermap(-np.log(enrichment_df_corrected), row_cluster=False,
                        annot=enrichment_df_sigs.astype(str), cmap="coolwarm",
                        figsize=(12,5), vmin=0, vmax=25, fmt='')
    plt.savefig(args.outfile, dpi=300)


if __name__ == '__main__':
    main()

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

CMAP=plt.get_cmap("tab20")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--cluster_bed",
                        action="store", type=str, required=True,
                        dest="cluster_bed",
                        help="Input bed file with cluster assignments.")
    parser.add_argument("-a", "--atac_tsv",
                      action = "store", type=str, required=True,
                      dest = "atac_tsv", 
                      help = "TSV file contain per cell type pseudobulked ATAC data")
    parser.add_argument("--orientation",
                        action="store", type=str, required=False,
                        default="horizontal",
                        dest="orientation")
    parser.add_argument("--cell_type_order",
                        action="store", type=str, required=False,
                        default=None, dest="cell_type_order")
    parser.add_argument("--cluster_order",
                        action="store", type=str, required=False,
                        default=None, dest="cluster_order")
    parser.add_argument("-o", "--outfile",
                        action="store", type=str, required=True,
                        dest="outfile",
                        help="Output file to save plot.")
    args = parser.parse_args()

    cellbypeaks_logcpm = pd.read_csv(args.atac_tsv, sep="\t", index_col=0) 
    bed_df = pd.read_csv(args.cluster_bed,
                         sep="\t", names=["chr", "start", "end", "k"])
    
    if args.cluster_order is not None:
        cluster_order = [int(i) for i in args.cluster_order.split(",")]
        
    else:
        cluster_order = np.unique(bed_df["k"].values)
    cluster_labels = pd.Series(pd.Categorical(bed_df["k"].values, 
                                              categories=cluster_order, ordered=True),
                                index=cellbypeaks_logcpm.index)
    cellbypeaks_logcpm = cellbypeaks_logcpm.loc[cluster_labels.sort_values().index]
    cluster_labels = cluster_labels.loc[cellbypeaks_logcpm.index]
    cluster_color_mapping = {c: CMAP(i) for i, c in enumerate(cluster_labels.unique())}
    cluster_colors = cluster_labels.astype(int).apply(lambda x: cluster_color_mapping[int(x)])

    if args.cell_type_order is not None:
        cell_type_order = args.cell_type_order.split(",")
        cellbypeaks_logcpm = cellbypeaks_logcpm[cell_type_order]

    if args.orientation == "horizontal":
        sns.clustermap(cellbypeaks_logcpm.T, 
                    row_cluster=True if args.cell_type_order is None else False,
                    col_cluster=False, 
                    col_colors=cluster_colors, 
                    xticklabels=False, cmap="Blues", figsize=(16,12))
    else:
        sns.clustermap(cellbypeaks_logcpm, 
                    row_cluster=False, 
                    col_cluster=True if args.cell_type_order is None else False, 
                    row_colors=cluster_colors, 
                    yticklabels=False, cmap="Blues", figsize=(12,16))
                    
    patches = [mpatches.Patch(color=v, label=k) for k, v in cluster_color_mapping.items()]
    if args.orientation == "horizontal":
        plt.legend(handles=patches, loc=(2,0.1))
    else:
        plt.legend(handles=patches, loc="lower left")
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=300)


if __name__ == '__main__':
    main()
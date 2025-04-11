import argparse
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
    parser.add_argument("-o", "--outfile",
                        action="store", type=str, required=True,
                        dest="outfile",
                        help="Output file to save plot.")
    args = parser.parse_args()

    cellbypeaks_logcpm = pd.read_csv(args.atac_tsv, sep="\t", index_col=0) 
    bed_df = pd.read_csv(args.cluster_bed,
                         sep="\t", names=["chr", "start", "end", "k"])
    cluster_labels = pd.Series(bed_df["k"].values, index=cellbypeaks_logcpm.index)
    cluster_color_mapping = {c: CMAP(i) for i, c in enumerate(cluster_labels.unique())}
    cluster_colors = cluster_labels.apply(lambda x: cluster_color_mapping[x])
    
    sns.clustermap(cellbypeaks_logcpm.loc[cluster_labels.sort_values().index].T, 
                   row_cluster=True, col_cluster=False, col_colors=cluster_colors, 
                   xticklabels=False, cmap="Blues", figsize=(16,12))
    patches = [mpatches.Patch(color=v, label=k) for k, v in cluster_color_mapping.items()]
    plt.legend(handles=patches, loc=(2,0.1))
    plt.savefig(args.outfile, dpi=300)


if __name__ == '__main__':
    main()
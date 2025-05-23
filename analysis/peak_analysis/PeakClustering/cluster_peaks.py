import argparse
import re
import json
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--atac_tsv",
                      action = "store", type=str, required=True,
                      dest = "atac_tsv", 
                      help = "TSV file contain per cell type pseudobulked ATAC data")
    parser.add_argument("-o", "--outdir",
                      action = "store", type=str, required=True,
                      dest = "outdir",
                      help = "Output directory to save files to.")
    parser.add_argument("-p", "--prefix",  
                        action = "store", type=str, required=True,
                        dest = "file_prefix",
                        help = "File prefix to use too save files.")
    parser.add_argument("-k",
                        action = "store", type=int, required=True,
                        dest = "k",
                        help = "Number of clusters to use for clustering.")
    parser.add_argument("-m", "--method",
                        action = "store", type=str, required=False,
                        default = "kmeans",
                        dest = "method",
                        help = "Clustering method to use.")
    args = parser.parse_args()
    
    cellbypeaks_logcpm = pd.read_csv(args.atac_tsv, sep="\t", index_col=0)
    if args.method == "kmeans":
        kmeans = KMeans(n_clusters=args.k, random_state=0, n_init="auto").fit(cellbypeaks_logcpm)
        cluster_labels = kmeans.labels_
        with open(f"{args.outdir}/{args.method}_{args.k}_{args.file_prefix}_WSS_error.json", 'w') as f:
            json.dump({"Method": args.method,
                       "K": args.k,
                       "WSS": kmeans.inertia_}, f)
    else:
        nmf = NMF(n_components=args.k, init="nndsvd", random_state=0, verbose=True, max_iter=1_000).fit(cellbypeaks_logcpm)
        W = nmf.transform(cellbypeaks_logcpm)
        cluster_labels = W.argmax(axis=1)
        with open(f"{args.outdir}/{args.method}_{args.k}_{args.file_prefix}_WSS_error.json", 'w') as f:
            json.dump({"Method": args.method,
                       "K": args.k,
                       "Reconstructioin error": nmf.reconstruction_err_}, f)

    # save results to bed file
    bed_df = pd.DataFrame([re.split("[:-]", i) for i in cellbypeaks_logcpm.index],
                          columns=["chr", "start", "end"])
    bed_df["k"] = cluster_labels
    bed_df.to_csv(f"{args.outdir}/{args.method}_{args.k}_{args.file_prefix}.bed", 
                      sep="\t", header=False, index=False)


if __name__ == '__main__':
    main()
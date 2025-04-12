import argparse
import os
from tqdm import tqdm
import pandas as pd
import pyranges as pr
from scipy import stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--qtldir",
                      action = "store", type=str, required=False,
                      dest = "qtl_dir", default = "/data_vol/FINAL/QTL/rna/L1_manual",
                      help = "QTL dir") 
    parser.add_argument("-b", "--bed",
                      action = "store", type=str, required=False,
                      dest = "bed",
                      help = "Cell type specific peak bed file") 
    parser.add_argument("-o", "--outfile",
                      action = "store", type=str, required=False,
                      dest = "outfile",
                      help = "outfile") 
    (options, args) = parser.parse_args()
    
    peak_bed = pd.read_csv(args.bed, sep="\t", names=["Chromosome", "Start", "End", "Cluster"])
    clusters = peak_bed["Cluster"].unique()
    cell_types = os.listdir(args.qtl_dir)
    cell_types = [c for c in cell_types 
                  if os.path.exists(f"{args.qtl_dir}/{c}/{c}.cis_qtl_pairs.1.parquet")]
    
    enrichment_df = pd.DataFrame([], index=cell_types, columns=clusters)
    for cell_type in cell_types:
        print(cell_type)
        qtl_df = []
        for chr in tqdm(range(1, 23)):
            qtl_df.append(pd.read_parquet(f"{args.qtl_dir}/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet", 
            columns=["pval_nominal", "chr", "pos"]))
        qtl_df = pd.concat(qtl_df)
        qtl_df.rename(columns={"chr": "Chromosome", "pos": "End"}, inplace=True)
        qtl_df["Start"] = qtl_df["End"] - 1
        qtl_df_pr = pr.PyRanges(qtl_df)
        
        for cluster in clusters:
            cluster_pr = pr.PyRanges(peak_bed[peak_bed["Cluster"] == cluster])
            qtl_df[cluster] = False
            overlap_inds = qtl_df_pr.overlap(cluster_pr).df["index"].values
            qtl_df.loc[overlap_inds, cluster] = True
    
            enrichment_df.loc[cell_type, cluster] = stats.mannwhitneyu(
                qtl_df[qtl_df[cluster] == True]["pval_nominal"].values,
                qtl_df[qtl_df[cluster] == False]["pval_nominal"].values,
                alternative="less")[1]
            
    enrichment_df.to_csv(args.outfile, sep="\t", header=True, index=True)


if __name__ == '__main__':
    main()
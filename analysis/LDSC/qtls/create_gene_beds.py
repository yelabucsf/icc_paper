import argparse
import os
from tqdm import tqdm
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--qtldir",
                      action = "store", type=str, required=False,
                      dest = "qtl_dir", default = "/data_vol/FINAL/QTL/rna/L2_manual",
                      help = "QTL dir") 
    parser.add_argument("-c", "--celltype",
                      action = "store", type=str, required=False,
                      dest = "cell_type", default = "B_Memory",
                      help = "cell type") 
    parser.add_argument("-g", "--gene_locs",
                      action = "store", type=str, required=False,
                      dest = "gene_locs", default = "/data_vol/References/gene_locs_autosome_clean.txt",
                      help = "gene locs file") 
    args = parser.parse_args()

    # read in QTLs
    sig_genes = []
    for chr in tqdm(range(1, 23)):
        tmp_df = pd.read_csv(f"{args.qtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_egene.{chr}.txt", sep="\t")
        tmp_df = tmp_df.loc[tmp_df["qval"] < 0.05][["phenotype_id"]]
        sig_genes.append(tmp_df)
    sig_genes = pd.concat(sig_genes, ignore_index=True)

    gene_locs = pd.read_csv(args.gene_locs, sep="\t")
    gene_locs = gene_locs.loc[gene_locs["gene"].isin(sig_genes["phenotype_id"].values)]
    gene_locs["s1"]  = gene_locs["s1"] - 100_000
    gene_locs["s2"]  = gene_locs["s2"] + 100_000
    gene_locs[["chr", "s1", "s2"]].to_csv(f"{args.qtl_dir}/{args.cell_type}/{args.cell_type}.significant_genes_slop_100kb.bed", sep="\t", header=False, index=False)

            
if __name__ == "__main__":
    main()
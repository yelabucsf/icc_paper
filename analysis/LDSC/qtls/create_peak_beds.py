import argparse
import os
from tqdm import tqdm
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--qtldir",
                      action = "store", type=str, required=False,
                      dest = "qtl_dir", default = "/data_vol/FINAL/QTL/atac/permissive/L2_manual",
                      help = "QTL dir") 
    parser.add_argument("-c", "--celltype",
                      action = "store", type=str, required=False,
                      dest = "cell_type", default = "B_Memory",
                      help = "cell type") 
    args = parser.parse_args()

    # read in QTLs
    sig_peaks = []
    for chr in tqdm(range(1, 23)):
        tmp_df = pd.read_csv(f"{args.qtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_egene.{chr}.txt", sep="\t")
        tmp_df = tmp_df.loc[tmp_df["qval"] < 0.05][["phenotype_id"]]
        sig_peaks.append(tmp_df)
    sig_peaks = pd.concat(sig_peaks, ignore_index=True)
    sig_peaks["chr"] = sig_peaks["phenotype_id"].str.split(":").str[0]
    sig_peaks["start"] = sig_peaks["phenotype_id"].str.split(":").str[1].str.split("-").str[0].astype(int) - 250
    sig_peaks["end"] = sig_peaks["phenotype_id"].str.split(":").str[1].str.split("-").str[1].astype(int) + 250
    sig_peaks[["chr", "start", "end"]].to_csv(f"{args.qtl_dir}/{args.cell_type}/{args.cell_type}.significant_peaks_1000bp.bed", sep="\t", header=False, index=False)

            
if __name__ == "__main__":
    main()
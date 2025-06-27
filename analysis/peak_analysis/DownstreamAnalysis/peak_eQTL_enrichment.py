import argparse
import os
from tqdm import tqdm
import numpy as np
import pandas as pd
import pyranges as pr
from scipy import stats
import boto3


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
    parser.add_argument("-m", "--modality",
                      action = "store", type=str, required=False,
                      dest = "modality", default="rna",
                      help = "modality") 
    args = parser.parse_args()
    
    s3 = boto3.client('s3')
    s3_bucket = "gracie-icc"

    peak_bed = pd.read_csv(args.bed, sep="\t", names=["Chromosome", "Start", "End", "Cluster"])
    clusters = peak_bed["Cluster"].unique()
    cell_types = os.listdir(args.qtl_dir)
    cell_types = [c for c in cell_types 
                  if os.path.exists(f"{args.qtl_dir}/{c}/{c}.cis_qtl_egene.1.txt")]
    print(cell_types)
    
    enrichment_df = pd.DataFrame([], index=cell_types, columns=clusters)
    for cell_type in cell_types:
        print(cell_type)
        qtl_df = []
        for chr in tqdm(range(1, 23)):
            if args.modality == "rna":
                qtl_df.append(pd.read_parquet(f"{args.qtl_dir}/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet", 
                columns=["pval_nominal", "chr", "pos"]))
            elif args.modality == "atac":
                level = "L1" if "L1" in args.qtl_dir else "L2"
                rna_dir = f"/data_vol/FINAL/QTL/rna/{level}_manual"
                if not os.path.exists(f"{rna_dir}/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet"):
                    s3.download_file(s3_bucket, 
                                     f"production_data/analysis_outputs/FINAL/QTL/rna/{level}_manual/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet", 
                                     f"{rna_dir}/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet")
                num_rows = pd.read_parquet(f"{rna_dir}/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet", 
                                           columns=["pval_nominal"]).shape[0]
                os.remove(f"{rna_dir}/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet")

                if not os.path.exists(f"{args.qtl_dir}/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet"):
                    s3.download_file(s3_bucket, 
                                    f"production_data/analysis_outputs/FINAL/QTL/atac/permissive/{level}_manual/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet", 
                                    f"{args.qtl_dir}/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet")

                tmp_df = pd.read_parquet(f"{args.qtl_dir}/{cell_type}/{cell_type}.cis_qtl_pairs.{chr}.parquet", 
                                         columns=["phenotype_id", "start_distance", "pval_nominal"])
                tmp_df = tmp_df.sample(num_rows, random_state=42)
                tmp_df["chr"] = f"chr{chr}"
                tmp_df["pos"] = tmp_df["phenotype_id"].str.split(":").str[1].str.split("-").str[0].astype(int)
                tmp_df["pos"] = tmp_df["pos"] + tmp_df["start_distance"]
                tmp_df.drop(columns=["phenotype_id", "start_distance"], inplace=True)
                qtl_df.append(tmp_df)

        qtl_df = pd.concat(qtl_df, ignore_index=True)
        qtl_df.rename(columns={"chr": "Chromosome", "pos": "End"}, inplace=True)
        qtl_df["Start"] = qtl_df["End"] - 1
        qtl_df.reset_index(inplace=True)
        qtl_df_pr = pr.PyRanges(qtl_df)
        print("Done loading qtl data")

        for cluster in tqdm(clusters):
            cluster_pr = pr.PyRanges(peak_bed[peak_bed["Cluster"] == cluster])
            qtl_df[cluster] = False
            overlap_inds = qtl_df_pr.overlap(cluster_pr).df["index"].values
            qtl_df.loc[overlap_inds, cluster] = True
    
            enrichment_df.loc[cell_type, cluster] = stats.mannwhitneyu(
                qtl_df[qtl_df[cluster] == True]["pval_nominal"].values,
                qtl_df[qtl_df[cluster] == False]["pval_nominal"].values,
                alternative="less")[1]
            enrichment_df.to_csv(args.outfile, sep="\t", header=True, index=True)  # save intermediate progress
            
    enrichment_df.to_csv(args.outfile, sep="\t", header=True, index=True)


if __name__ == '__main__':
    main()

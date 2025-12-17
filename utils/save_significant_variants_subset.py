import pandas as pd
import boto3
import os
from tqdm import tqdm


def main():
    data_dir_dict = {
        "rna": {
            "s3": "production_data/analysis_outputs/FINAL/QTL/rna",
            "local": "/data_vol/FINAL/QTL/rna"
        },
        "atac": {
            "s3": "production_data/analysis_outputs/FINAL/QTL/atac/permissive",
            "local": "/data_vol/FINAL/QTL/atac/permissive"
        }
    }
    cell_types_dict = {
        "Cohort": ["ICC"],
        "L1_manual": ["B", "CD4T", "CD8T", "NK", "cDC", "cM", "gdT", "ncM"],
        "L2_manual": ["B_Memory", "B_Naive", "CD4T_Memory", "CD4T_Naive", "CD8T_Memory", "CD8T_Naive", "MAIT", "NK_dim", "Treg", "cDC2", "cM", "gdT", "ncM"]
    }
    
    s3 = boto3.client('s3')
    s3_bucket = "gracie-icc"

    for modality, data_dirs in data_dir_dict.items():
        for level, cell_types in cell_types_dict.items():
            for cell_type in cell_types:
                print(f"Processing {modality} {level} {cell_type}")
                sig_df = []
                for chr in tqdm(range(1, 23)):
                    if not os.path.exists(f"{data_dirs['local']}/{level}/{cell_type}/{cell_type}.cis_qtl_egene.{chr}.txt"):
                        s3.download_file(s3_bucket, 
                                         f"{data_dirs['s3']}/{level}/{cell_type}/{cell_type}.cis_qtl_egene.{chr}.txt", 
                                         f"{data_dirs['local']}/{level}/{cell_type}/{cell_type}.cis_qtl_egene.{chr}.txt")
                    egenes = pd.read_csv(f"{data_dirs['local']}/{level}/{cell_type}/{cell_type}.cis_qtl_egene.{chr}.txt", sep="\t")
                    egenes = egenes.loc[egenes["qval"] < 0.05]
                    
                    if not os.path.exists(f"{data_dirs['local']}/{level}/{cell_type}/{cell_type}.cis_qtl_snp.{chr}.annot.txt"):
                        s3.download_file(s3_bucket, 
                                         f"{data_dirs['s3']}/{level}/{cell_type}/{cell_type}.cis_qtl_snp.{chr}.annot.txt", 
                                         f"{data_dirs['local']}/{level}/{cell_type}/{cell_type}.cis_qtl_snp.{chr}.annot.txt")

                    df = pd.read_csv(f"{data_dirs['local']}/{level}/{cell_type}/{cell_type}.cis_qtl_snp.{chr}.annot.txt", sep="\t")
                    df = df.loc[df["phenotype_id"].isin(egenes["phenotype_id"].values)]
                    df = df.merge(egenes[["phenotype_id", "pval_nominal_threshold"]], on="phenotype_id", how="left")
                    df = df.loc[df["pval_nominal"] < df["pval_nominal_threshold"]]
                    sig_df.append(df)
                    os.remove(f"{data_dirs['local']}/{level}/{cell_type}/{cell_type}.cis_qtl_snp.{chr}.annot.txt")
                sig_df = pd.concat(sig_df, ignore_index=True)
                sig_df.to_csv(f"{data_dirs['local']}/{level}/{cell_type}/{cell_type}.cis_qtl_snp.all_significant.annot.txt", sep="\t", index=False)
                s3.upload_file(f"{data_dirs['local']}/{level}/{cell_type}/{cell_type}.cis_qtl_snp.all_significant.annot.txt", 
                               s3_bucket, 
                               f"{data_dirs['s3']}/{level}/{cell_type}/{cell_type}.cis_qtl_snp.all_significant.annot.txt")


if __name__ == "__main__":
    main()
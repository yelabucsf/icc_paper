import argparse
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
    parser.add_argument("-f", "--finemap",
                      action = "store_true", required=False,
                      dest = "finemap", default = False,
                      help = "Whether to create bed file for finemapping results") 
    args = parser.parse_args()

    # read in QTLs
    qtls = pd.read_csv(f"{args.qtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_snp.all_significant.annot.txt", sep="\t")
    if "chr" not in qtls.columns:
        if "chrom" in qtls.columns:
            qtls["chr"] = qtls["chrom"]
        else:
            raise ValueError("chr or chrom not in qtls")
        if not qtls["chr"].astype(str).str.contains("chr").all():
            qtls["chr"] = "chr" + qtls["chr"].astype(str)
        
    qtls["start"] = qtls["pos"] - 1
    qtls[["chr", "start", "pos"]].drop_duplicates().to_csv(f"{args.qtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_pairs.significant.bed", sep="\t", header=False, index=False)

    if args.finemap:
        if "PIP" not in qtls.columns:
            egenes = qtls["phenotype_id"].unique()
            susie_results = []
            for gene in egenes:
                try:
                    tmp_df = pd.read_csv(f"{args.qtl_dir}/{args.cell_type}/SuSiE_new/{args.cell_type}_{gene}_PIP.csv.gz", compression="gzip")
                    tmp_df.rename(columns={"SNP": "variant_id"}, inplace=True)
                    tmp_df["phenotype_id"] = gene
                    susie_results.append(tmp_df)
                except Exception as e:
                    print(e)
                    print(f"No results for {gene}")
            susie_results = pd.concat(susie_results)
            qtls = qtls.merge(susie_results, on=["phenotype_id", "variant_id"], how="left")
        qtls.loc[qtls["CS"] == 1, ["chr", "start", "pos"]].drop_duplicates().to_csv(f"{args.qtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_pairs.significant.SuSiE_CS.bed", sep="\t", header=False, index=False)
        for pip in [0.3]:
            qtls.loc[qtls["PIP"] >= pip, ["chr", "start", "pos"]].drop_duplicates().to_csv(f"{args.qtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_pairs.significant.SuSiE_PIP_{pip}.bed", sep="\t", header=False, index=False)

            
if __name__ == "__main__":
    main()
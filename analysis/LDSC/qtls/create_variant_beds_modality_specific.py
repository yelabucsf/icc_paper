import argparse
from tqdm import tqdm
import pandas as pd


cell_types_dict = {"Cohort": ["ICC"],
                   "L1_manual": ["B", "CD4T", "CD8T", "NK", "cDC", "cM", "gdT", "ncM"],
                   "L2_manual": ["B_Memory", "B_Naive", "CD4T_Memory", "CD4T_Naive", "CD8T_Memory", "CD8T_Naive", "MAIT", "NK_dim", "Treg", "cDC2", "cM", "gdT", "ncM"]}

RNA_QTL_dir = "/data_vol/FINAL/QTL/rna"
ATAC_QTL_dir = "/data_vol/FINAL/QTL/atac/permissive"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--level",
                      action = "store", type=str, required=False,
                      dest = "level", default = "L2_manual",
                      help = "Level of cell type resolution") 
    args = parser.parse_args()

    rna_qtls = pd.read_csv(f"{RNA_QTL_dir}/{args.level}/all_cell_types.significant_qtls_pooled.bed", sep="\t", names=["chrom", "start", "end"])
    atac_qtls = pd.read_csv(f"{ATAC_QTL_dir}/{args.level}/all_cell_types.significant_qtls_pooled.bed", sep="\t", names=["chrom", "start", "end"])
    
    df_all = rna_qtls.merge(atac_qtls, how='outer', indicator=True)

    df_all[df_all["_merge"] == "both"].to_csv(f"{RNA_QTL_dir}/{args.level}/all_cell_types.significant_qtls_pooled_rna_and_atac.bed", sep="\t", header=False, index=False)
    df_all[df_all["_merge"] == "left_only"].to_csv(f"{RNA_QTL_dir}/{args.level}/all_cell_types.significant_qtls_pooled_rna_only.bed", sep="\t", header=False, index=False)
    df_all[df_all["_merge"] == "right_only"].to_csv(f"{RNA_QTL_dir}/{args.level}/all_cell_types.significant_qtls_pooled_atac_only.bed", sep="\t", header=False, index=False)

            
if __name__ == "__main__":
    main()
import argparse
from tqdm import tqdm
import pandas as pd


cell_types_dict = {"Cohort": ["ICC"],
                   "L1_manual": ["B", "CD4T", "CD8T", "NK", "cDC", "cM", "gdT", "ncM"],
                   "L2_manual": ["B_Memory", "B_Naive", "CD4T_Memory", "CD4T_Naive", "CD8T_Memory", "CD8T_Naive", "MAIT", "NK_dim", "Treg", "cDC2", "cM", "gdT", "ncM"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--qtldir",
                      action = "store", type=str, required=False,
                      dest = "qtl_dir", default = "/data_vol/FINAL/QTL/rna",
                      help = "QTL dir") 

    parser.add_argument("-l", "--level",
                      action = "store", type=str, required=False,
                      dest = "level", default = "L2_manual",
                      help = "Level of cell type resolution") 
    args = parser.parse_args()

    # read in QTLs
    qtls = {}
    for cell_type in tqdm(cell_types_dict[args.level]):
        qtls[cell_type] = pd.read_csv(f"{args.qtl_dir}/{args.level}/{cell_type}/{cell_type}.cis_qtl_snp.all_significant.annot.txt", sep="\t")
        qtls[cell_type] = qtls[cell_type].drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
        qtls[cell_type]["Cell type"] = cell_type

    qtls = pd.concat(qtls.values())
    qtls["start"] = qtls["pos"] - 1
    if not qtls["chrom"].astype(str).str.contains("chr").all():
        qtls["chrom"] = "chr" + qtls["chrom"].astype(str)

    qtl_counts = qtls[["chrom", "start", "pos", "ref", "alt"]].value_counts()

    # save bed file with all QTLs across all cell types
    qtls[["chrom", "start", "pos"]].drop_duplicates().to_csv(f"{args.qtl_dir}/{args.level}/all_cell_types.significant_qtls_pooled.bed", 
                                                             sep="\t", header=False, index=False)

    # save bed file with all QTLs specific to one cell type
    qtl_counts[qtl_counts == 1].index.to_frame()[["chrom", "start", "pos"]].drop_duplicates().to_csv(f"{args.qtl_dir}/{args.level}/all_cell_types.significant_qtl_in_one_cell_type.bed", 
                                                                                                     sep="\t", header=False, index=False)

    # save bed file with all QTLs shared in two or more cell types
    qtl_counts[qtl_counts > 1].index.to_frame()[["chrom", "start", "pos"]].drop_duplicates().to_csv(f"{args.qtl_dir}/{args.level}/all_cell_types.significant_qtl_in_multiple_cell_types.bed", 
                                                                                                     sep="\t", header=False, index=False)

            
if __name__ == "__main__":
    main()
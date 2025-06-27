from argparse import ArgumentParser
from tqdm import tqdm
import pandas as pd

def main():
    parser = ArgumentParser()
    parser.add_argument("--maffile", type=str, required=True)
    parser.add_argument("--hg38_snps", type=str, required=True)
    parser.add_argument("--infile1", type=str, required=True)
    parser.add_argument("--infile2", type=str, required=True)
    parser.add_argument("--infile3", type=str, required=True)

    args = parser.parse_args()

    d_var = pd.read_csv(args.infile1, sep=" ", header=None)
    d_pos1 = pd.read_csv(args.infile2, sep=" ", header=None)
    d_pos2 = pd.read_csv(args.infile3, sep="\t", header=None)
    d_maf = pd.read_csv(args.maffile, sep="\t")
    with open(args.hg38_snps, "r") as f:
        d_hg38_snps = f.read().splitlines() 

    d_var.columns=["var","gene","P","snp"]
    d_var["chr"] = d_var["var"].str.split("_", expand=True)[0]
    d_var["pos_38"] = d_var["var"].str.split("_", expand=True)[1].astype(int)
    d_var["A1"] = d_var["var"].str.split("_", expand=True)[2]
    d_var["A2"] = d_var["var"].str.split("_", expand=True)[3]
    d_var["build"] = d_var["var"].str.split("_", expand=True)[4]

    d_pos1 = d_pos1[[0,2,3]]
    d_pos1.columns = ["chr", "pos_38", "snp"]

    d_pos2 = d_pos2[[0,2,3]]
    d_pos2.columns = ["chr", "pos_19", "snp"]
    d_pos = d_pos1.merge(d_pos2, on=["chr", "snp"])

    df = d_var.merge(d_pos, on=["chr", "pos_38", "snp"])
    df["chr_num"] = df["chr"].str.replace("[a-z]", "", regex=True).astype(int)
    df = df[(df["chr_num"] >= 1) & (df["chr_num"] <= 22)]

    df["allele1_hg38"] = df["chr_num"].astype(str) + ":" + df["pos_38"].astype(str) + ":" + df["A1"] + ":" + df["A2"]
    df["allele2_hg38"] = df["chr_num"].astype(str) + ":" + df["pos_38"].astype(str) + ":" + df["A2"] + ":" + df["A1"]
    dg1 = df[df["allele1_hg38"].isin(d_hg38_snps)]
    dg1["allele_hg38"] = dg1["allele1_hg38"]
    dg2 = df[df["allele2_hg38"].isin(d_hg38_snps)]
    dg2["allele_hg38"] = dg2["allele2_hg38"]
    dg = pd.concat([dg1, dg2])

    dg["allele1_hg19"] = dg["chr_num"].astype(str) + ":" + dg["pos_19"].astype(str) + ":" + dg["A1"] + ":" + dg["A2"]
    dg["allele2_hg19"] = dg["chr_num"].astype(str) + ":" + dg["pos_19"].astype(str) + ":" + dg["A2"] + ":" + dg["A1"]
    
    dh1 = dg[dg["allele1_hg19"].isin(d_maf["SNP"])]
    dh1["allele_hg19"] = dh1["allele1_hg19"]
    dh2 = dg[dg["allele2_hg19"].isin(d_maf["SNP"])]
    dh2["allele_hg19"] = dh2["allele2_hg19"]
    dh = pd.concat([dh1, dh2])

    dh = dh[["chr_num", "gene", "allele_hg38", "P"]]
    dh.columns = ["CHR", "GENE", "SNP", "P"]
    dh.sort_values(by="GENE", inplace=True)
    dh.to_csv(args.infile1, sep=" ", index=False, header=False)

if __name__ == "__main__":
    main()
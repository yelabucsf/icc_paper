import argparse
import os
from tqdm import tqdm
from dataclasses import dataclass, astuple
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt

# chromHMM states in PBMCs
# downloaded from https://www.encodeproject.org/annotations/ENCSR333EIL/
chromhmm_bed = "/data_vol/References/chromHMM/ENCFF863IVL.bed.gz"

#chromHMM state mapping
chromhmm_state_mapping_coarse = {"TssA": "TssActive",
    "TssAFlnk": "TssActive",
    "TssFlnk": "TssActive",
    "TssFlnkD": "TssActive",
    "TssFlnkU": "TssActive",
    "TxFlnk": "Transcript",
    "Tx":"Transcript",
    "TxWk": "Transcript",
    "EnhG1":"Enhancer",
    "EnhG2":"Enhancer",
    "EnhA1":"Enhancer",
    "EnhA2":"Enhancer",
    "Enh": "Enhancer",
    "EnhWk": "Enhancer",
    "ZNF/Rpts":"Heterochromatin",
    "Het":"Heterochromatin",
    "TssBiv": "TssBivalent",
    "BivFlnk": "TssBivalent",
    "EnhBiv": "Enhancer",
    "ReprPC": "Repressive",
    "ReprPCWk": "Repressive",
    "Quies": "Quiescent"}


chromhmm_state_mapping_granular = {"TssA": "TSS",
    "TssAFlnk": "Flanking TSS",
    "TssFlnk": "Flanking TSS",
    "TssFlnkD": "Flanking TSS",
    "TssFlnkU": "Flanking TSS",
    "TxFlnk": "Flanking Transcript",
    "Tx":"Transcript",
    "TxWk": "Weak Transcript",
    "EnhG1":"Genic Enhancer",
    "EnhG2":"Genic Enhancer",
    "EnhA1":"Active Enhancer",
    "EnhA2":"Active Enhancer",
    "Enh": "Enhancer",
    "EnhWk": "Weak Enhancer",
    "ZNF/Rpts":"Heterochromatin",
    "Het":"Heterochromatin",
    "TssBiv": "Bivalent TSS",
    "BivFlnk": "Bivalent TSS",
    "EnhBiv": "Bivalent Enhancer",
    "ReprPC": "Repressive",
    "ReprPCWk": "Repressive",
    "Quies": "Quiescent"}

@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: int
    alt: int

    @classmethod
    def create_from_row(cls, row):
        chrom = f"chr{row[0]}"
        pos = int(row[1])
        ref = row[3]
        alt = row[4]
        assert "," not in alt, alt
        return cls(chrom, pos, ref, alt)


def read_variants_from_pvar(path: str) -> dict[Variant]:
    open_fn = gzip.open if path.endswith(".gz") else open
    variants = dict()
    with open_fn(path, "rt") as f:
        for line in tqdm(f, desc="Reading variants"):
            if line.startswith("#"):
                continue
            row = line.strip().split("\t")
            variants[row[2]] = Variant.create_from_row(row)
    return variants


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--eqtldir",
                      action = "store", type=str, required=False,
                      dest = "eqtl_dir", default = "/data_vol/FINAL/QTL/rna/L2_manual",
                      help = "eQTL dir") 
    parser.add_argument("-c", "--caqtldir",
                      action = "store", type=str, required=False,
                      dest = "caqtl_dir", default = "/data_vol/FINAL/QTL/atac/permissive/L2_manual",
                      help = "caQTL dir") 
    parser.add_argument("--celltype",
                      action = "store", type=str, required=False,
                      dest = "cell_type", default = "B_Memory",
                      help = "cell type") 
    parser.add_argument("-o", "--outfile",
                      action = "store", type=str, required=False,
                      dest = "outfile",
                      help = "outfile") 
    args = parser.parse_args()

    # read in chromHMM sttes
    chromhmm_states = pr.read_bed(chromhmm_bed).df
    chromhmm_states["chromHMM State"] = chromhmm_states["Name"].apply(lambda x: chromhmm_state_mapping_granular[x])
    chromhmm_states = pr.PyRanges(chromhmm_states) 

    # read in caQTLs
    caqtls = pd.read_csv(f"{args.caqtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_snp.all_significant.annot.txt", sep="\t")

    # read in eQTLs
    eqtls = pd.read_csv(f"{args.eqtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_snp.all_significant.annot.txt", sep="\t")

    # add chromHMM annotations to QTLs
    qtls = {"eQTLs": eqtls, "caQTLs": caqtls}
    for k in qtls.keys():
        qtls[k].rename(columns={"chr": "Chromosome", "pos": "End"}, inplace=True)
        qtls[k]["Start"] = qtls[k]["End"] - 1

        # Add chromHMM annotations
        qtls[k] = pr.PyRanges(qtls[k]).join(chromhmm_states, how="left", suffix="_chromHMM").df
    
    # compute enrichment
    prop_caqtls_chromhmm = (qtls["caQTLs"]["chromHMM State"].value_counts()/len(qtls["caQTLs"]))
    prop_eqtls_chromhmm = (qtls["eQTLs"]["chromHMM State"].value_counts()/len(qtls["eQTLs"]))
    enrichment_df = pd.DataFrame({"caQTL enrichment relative to eQTLs":prop_caqtls_chromhmm/prop_eqtls_chromhmm, 
                                    "eQTL enrichment relative to caQTLs":prop_eqtls_chromhmm/prop_caqtls_chromhmm})
    enrichment_df["Enrichment"] = enrichment_df.max(axis=1) - 1
    enrichment_df.loc[enrichment_df["eQTL enrichment relative to caQTLs"] > enrichment_df["caQTL enrichment relative to eQTLs"], "Enrichment"] *= -1
    enrichment_df = enrichment_df[enrichment_df.index != "-1"]
    enrichment_df.to_csv(f"{args.outfile}.tsv", sep="\t", header=True, index=True) 

    # plot
    fig, ax = plt.subplots(figsize=(4,4))
    enrichment_df["Enrichment"].sort_values(ascending=False).plot.bar()
    locs, labels = plt.yticks()
    locs, labels = plt.yticks(locs, [round(l - 1, 1) if l < 0 else round(l + 1, 1) for l in locs])
    plt.plot([-1,20], [0,0], c="black")
    plt.ylabel("Enrichment", fontsize=12)
    plt.xlabel("Enriched in caQTLs --> Enriched in eQTLs", fontsize=12)
    plt.title(f"chromHMM State Enrichment in caQTLs vs. eQTLs - {cell_type}")
    plt.savefig(f"{args.outfile}.png", dpi=300)
    plt.savefig(f"{args.outfile}.pdf", dpi=300)


if __name__ == '__main__':
    main()
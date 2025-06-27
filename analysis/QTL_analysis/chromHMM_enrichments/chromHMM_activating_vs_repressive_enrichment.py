import argparse
import numpy as np
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

    # read in QTLs
    caqtls = pd.read_csv(f"{args.caqtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_snp.all_significant.annot.txt", sep="\t")
    eqtls = pd.read_csv(f"{args.eqtl_dir}/{args.cell_type}/{args.cell_type}.cis_qtl_snp.all_significant.annot.txt.cis_qtl_snp.all_significant.annot.txt", sep="\t")
    merged_df = caqtls.merge(eqtls, on=["chr", "pos", "ref", "alt"], how="inner", suffixes=("_caqtl", "_eqtl"))
    merged_df["Shared QTL type"] = merged_df["slope_caqtl"].apply(lambda x: np.sign(x)) == merged_df["slope_eqtl"].apply(lambda x: np.sign(x))
    merged_df["Shared QTL type"] = merged_df["Shared QTL type"].map({True: "Activating", False: "Repressive"})

    # add chromHMM annotations to QTLs
    qtls = {"Activating": merged_df[merged_df["Shared QTL type"] == "Activating"], 
            "Repressive": merged_df[merged_df["Shared QTL type"] == "Repressive"]}
    for k in qtls.keys():
        qtls[k].rename(columns={"chr": "Chromosome", "pos": "End"}, inplace=True)
        qtls[k]["Start"] = qtls[k]["End"] - 1

        # Add chromHMM annotations
        qtls[k] = pr.PyRanges(qtls[k]).join(chromhmm_states, how="left", suffix="_chromHMM").df

    # compute enrichment
    prop_activating_chromhmm = (qtls["Activating"]["chromHMM State"].value_counts()/len(qtls["Activating"]))
    prop_repressive_chromhmm = (qtls["Repressive"]["chromHMM State"].value_counts()/len(qtls["Repressive"]))
    enrichment_df = pd.DataFrame({"Activating QTL enrichment relative to repressive QTLs":prop_activating_chromhmm/prop_repressive_chromhmm, 
                                "Repressive QTL enrichment relative to activating QTLs":prop_repressive_chromhmm/prop_activating_chromhmm})
    enrichment_df["Enrichment"] = enrichment_df.max(axis=1) - 1
    enrichment_df.loc[enrichment_df["Repressive QTL enrichment relative to activating QTLs"] > enrichment_df["Activating QTL enrichment relative to repressive QTLs"], "Enrichment"] *= -1
    enrichment_df = enrichment_df[enrichment_df.index != "-1"]
    enrichment_df = enrichment_df[enrichment_df.index != "Bivalent Enhancer"]
    enrichment_df.to_csv(f"{args.outfile}.tsv", sep="\t", header=True, index=True) 

    # plot
    fig, ax = plt.subplots(figsize=(4,4))
    enrichment_df["Enrichment"].sort_values(ascending=False).plot.bar()
    locs, labels = plt.yticks()
    locs, labels = plt.yticks(locs, [round(l - 1, 1) if l < 0 else round(l + 1, 1) for l in locs])
    plt.plot([-1,20], [0,0], c="black")
    plt.ylabel("Enrichment", fontsize=12)
    plt.xlabel("Enriched in activating QTLs --> Enriched in repressive QTLs", fontsize=12)
    plt.title(f"chromHMM State Enrichment in activating vs. repressive QTLs - {cell_type}")
    plt.savefig(f"{args.outfile}.png", dpi=300)
    plt.savefig(f"{args.outfile}.pdf", dpi=300)


if __name__ == "__main__":
    main()
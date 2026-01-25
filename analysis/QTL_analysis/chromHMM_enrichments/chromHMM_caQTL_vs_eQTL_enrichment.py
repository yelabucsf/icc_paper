import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns

# chromHMM states in PBMCs
# downloaded from https://www.encodeproject.org/annotations/ENCSR333EIL/
chromhmm_bed = "/data/References/chromHMM/ENCFF863IVL.bed.gz"

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


sns.set_context('talk',font_scale=1.4)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--eqtldir",
                      action = "store", type=str, required=False,
                      dest = "eqtl_dir", default = "/data/FINAL/QTL/rna/L2_manual",
                      help = "eQTL dir") 
    parser.add_argument("-c", "--caqtldir",
                      action = "store", type=str, required=False,
                      dest = "caqtl_dir", default = "/data/FINAL/QTL/atac/permissive/L2_manual",
                      help = "caQTL dir") 
    parser.add_argument("--celltype",
                      action = "store", type=str, required=False,
                      dest = "cell_type", default = "B_Memory",
                      help = "cell type") 
    parser.add_argument("-o", "--outfile_prefix",
                      action = "store", type=str, required=False,
                      dest = "outfile",
                      help = "outfile") 
    args = parser.parse_args()

    if not os.path.exists(f"{args.outfile}_odds_ratio.tsv"):
        # read in chromHMM states
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
            qtls[k].drop_duplicates(subset=["chrom", "pos"], inplace=True)
            qtls[k].rename(columns={"chrom": "Chromosome", "pos": "End"}, inplace=True)
            qtls[k]["Chromosome"] = qtls[k]["Chromosome"].apply(lambda x: f"chr{x}")
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
        enrichment_df.to_csv(f"{args.outfile}_fold_enrichment.tsv", sep="\t", header=True, index=True) 

        # compute odds ratio and p-value using Fisher's exact test
        overlapping_snps = np.intersect1d(qtls["caQTLs"]["variant_id"].unique(), qtls["eQTLs"]["variant_id"].unique())

        odds_ratio_df = []
        for state in chromhmm_states.df["chromHMM State"].unique():
            ncaqtls = len(qtls["caQTLs"][(~qtls["caQTLs"]["variant_id"].isin(overlapping_snps))])
            neqtls = len(qtls["eQTLs"][(~qtls["eQTLs"]["variant_id"].isin(overlapping_snps))])
            n_caqtls_state = len(qtls["caQTLs"][(qtls["caQTLs"]["chromHMM State"] == state) & (~qtls["caQTLs"]["variant_id"].isin(overlapping_snps))])
            n_eqtls_state = len(qtls["eQTLs"][(qtls["eQTLs"]["chromHMM State"] == state) & (~qtls["eQTLs"]["variant_id"].isin(overlapping_snps))])
            
            contingency_table = np.array([[n_caqtls_state, ncaqtls - n_caqtls_state],
                                        [n_eqtls_state, neqtls - n_eqtls_state]])
            caqtl_stats = stats.fisher_exact(contingency_table)

            contingency_table = np.array([[n_eqtls_state, neqtls - n_eqtls_state],
                                          [n_caqtls_state, ncaqtls - n_caqtls_state]])
            eqtl_stats = stats.fisher_exact(contingency_table)
            odds_ratio_df.append({"chromHMM state": state, "Odds ratio (caQTL)": caqtl_stats.statistic, 
                                  "Odds ratio (eQTL)": eqtl_stats.statistic, "p-value": caqtl_stats.pvalue})
        odds_ratio_df = pd.DataFrame(odds_ratio_df)
        odds_ratio_df["Odds ratio"] = odds_ratio_df[["Odds ratio (caQTL)", "Odds ratio (eQTL)"]].max(axis=1)
        _, pvals, _, _ = multipletests(odds_ratio_df["p-value"].values,
                                        alpha=0.01, method="fdr_bh")
        odds_ratio_df["p-value (corrected)"] = pvals
        odds_ratio_df.to_csv(f"{args.outfile}_odds_ratio.tsv", sep="\t", header=True, index=False)

    # plot fold enrichment results
    enrichment_df = pd.read_csv(f"{args.outfile}_fold_enrichment.tsv", sep="\t", index_col=0)

    fig, ax = plt.subplots(figsize=(8,7))
    enrichment_df["Enrichment"].sort_values(ascending=False).plot.barh(color="black")
    plt.plot([0,0], [-1,20], c="black", linewidth=0.5)

    xmax = max(abs(plt.xlim()[0]), abs(plt.xlim()[1]))
    xmax_rounded = np.ceil(xmax*2)/2  # ceiling to nearest 0.5
    plt.xlim(-xmax_rounded, xmax_rounded)
    plt.xticks(np.linspace(-xmax_rounded, xmax_rounded, 5))
    locs, labels = plt.xticks()
    locs, labels = plt.xticks(locs, [abs(l - 1) if l < 0 else l + 1 for l in locs])
    
    xlim = plt.xlim()
    ylim = plt.ylim()
    plt.fill_betweenx(ylim, [xlim[0], xlim[0]], [0, 0],  facecolor='tab:red', alpha=0.1)
    plt.fill_betweenx(ylim, [xlim[1], xlim[1]], [0, 0],  facecolor='tab:blue', alpha=0.1)

    plt.xlabel("Enrichment")
    plt.text(xlim[0] + 0.04, 0.1, "Enriched\nin eQTLs\nrelative\nto caQTLs", multialignment='center', fontsize=21, fontweight="bold", color="tab:red")
    plt.text(0.05, 8, "Enriched\nin caQTLs\nrelative\nto eQTLs", multialignment='center', fontsize=21, fontweight="bold", color="tab:blue")
    plt.title(args.cell_type.replace("_", " "))
    plt.tight_layout()
    plt.savefig(f"{args.outfile}_fold_enrichment.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{args.outfile}_fold_enrichment.pdf", dpi=300, bbox_inches="tight")

    # plot Odds ratio results
    odds_ratio_df = pd.read_csv(f"{args.outfile}_odds_ratio.tsv", sep="\t", index_col=0)
    odds_ratio_df["Odds ratio (plot vals)"] = odds_ratio_df.apply(lambda x: (x["Odds ratio"] - 1)*-1 if x["Odds ratio (eQTL)"] > 1 else x["Odds ratio"] - 1, axis=1)
    odds_ratio_df["Significance label"] = odds_ratio_df["p-value (corrected)"].apply(lambda x: "*" if x < 0.01 else " ")
    odds_ratio_df = odds_ratio_df.sort_values(by="Odds ratio (plot vals)", ascending=False)

    fig, ax = plt.subplots(figsize=(8,7))
    rects = odds_ratio_df["Odds ratio (plot vals)"].sort_values(ascending=False).plot.barh(color="black")
    ax.bar_label(rects.containers[0], odds_ratio_df["Significance label"].values, padding=3)
    plt.plot([0,0], [-1, 20], c="black", linewidth=0.5)

    xmax = max(abs(plt.xlim()[0]), abs(plt.xlim()[1]))
    xmax_rounded = np.ceil(xmax*2)/2  # ceiling to nearest 0.5
    plt.xlim(-xmax_rounded, xmax_rounded)
    plt.xticks(np.linspace(-xmax_rounded, xmax_rounded, 5))
    locs, labels = plt.xticks()
    locs, labels = plt.xticks(locs, [abs(l - 1) if l < 0 else l + 1 for l in locs])
    
    xlim = plt.xlim()
    ylim = plt.ylim()
    plt.fill_betweenx(ylim, [xlim[0], xlim[0]], [0, 0],  facecolor='tab:red', alpha=0.1)
    plt.fill_betweenx(ylim, [xlim[1], xlim[1]], [0, 0],  facecolor='tab:blue', alpha=0.1)

    plt.xlabel("Odds ratio")
    plt.text(xlim[0] + 0.04, 0.1, "Enriched\nin eQTLs\nrelative\nto caQTLs", multialignment='center', fontsize=21, fontweight="bold", color="tab:red")
    plt.text(0.05, 8, "Enriched\nin caQTLs\nrelative\nto eQTLs", multialignment='center', fontsize=21, fontweight="bold", color="tab:blue")
    plt.title(args.cell_type.replace("_", " "))
    plt.tight_layout()
    plt.savefig(f"{args.outfile}_odds_ratio.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{args.outfile}_odds_ratio.pdf", dpi=300, bbox_inches="tight")


if __name__ == '__main__':
    main()
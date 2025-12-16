#!/usr/bin/env python3

import pandas as pd


def parse_attributes(attr_str):
    """Parse the GFF3 attributes column into a dictionary."""
    attrs = {}
    for attr in attr_str.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key] = value
    return attrs


def pull_gene_id(attr_str):
    """Parse the GFF3 attributes column into a dictionary."""
    attr_l = attr_str.split(' ')
    gene_id_ix = attr_l.index("gene_id") + 1
    gene_id = attr_l[gene_id_ix][1:-2].split(".")[0]
    return gene_id


# get gene locations
gtf_url = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
gtf = pd.read_csv(gtf_url, sep='\t', comment="#", header=None)

genes_gtf = gtf.loc[gtf[2] == "gene"].copy()
gene_id = genes_gtf[8].apply(pull_gene_id)
genes_gtf.loc[:, "gene_id"] = gene_id
genes_gtf.loc[:, "tss"] = genes_gtf[3]
neg_k = genes_gtf[6] == '-'
genes_gtf.loc[neg_k, "tss"] = genes_gtf.loc[neg_k][4]
genes_gtf["tss0"] = genes_gtf["tss"] - 1

# get peak locations
peak_s3_url = "s3://gracie-icc/production_data/analysis_outputs/FINAL/references/peak_locs_L2_permissive.txt"
peaks_df = pd.read_csv(peak_s3_url, header=0, sep='\t')

# create BED files
bed_col_names = ["#chrom", "chromStart", "chromEnd", "name", "score", "strand"]

gene_bed = genes_gtf[[0, "tss0", "tss", "gene_id"]].copy()
gene_bed["score"] = 1000
gene_bed["strand"] = genes_gtf[6]
gene_bed.columns = bed_col_names

peak_bed = peaks_df[["chr", "s1", "s2", "gene"]].copy()
peak_bed["score"] = 1000
peak_bed["strand"] = "."
peak_bed.columns = bed_col_names
peak_bed.loc[:, "chromStart"] = peak_bed["chromStart"] - 1

# save
gene_bed_fn = "./dat/processed/cm_il1b/reqtl/gene_tss_loc.bed"
gene_bed.to_csv(gene_bed_fn, header=True, index=False, sep='\t')
peak_bed_fn = "./dat/processed/cm_il1b/reqtl/peak_loc.bed"
peak_bed.to_csv(peak_bed_fn, header=True, index=False, sep='\t')

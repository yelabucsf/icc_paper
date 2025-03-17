### Basic count matrix
### Take in H5 for rna or atac data and create pseudobulks per celltype
### change the cell by gene/peak matrix into a donor by gene/peak matrix
### This saves just the raw counts and some summary information about the runs 
### that are used for DEseq2

import argparse
import pandas as pd
import numpy as np
import anndata as adata
import scanpy as sc
from scipy.io import mmread
import os
#import pickle as pkl

#CMD: python make_basic_counts.py -R /data/march2025/MultiVI_posttraining_RNA_Celltypist_VaryingClusteringResolution_manual_annots.h5ad -O /data/processed/pseudobulk/rna/ -S best_sample L1_manual 

parser = argparse.ArgumentParser(description="This pseudobulk script accepts ATAC-seq, RNA-seq and/or protein sequencing data")
parser.add_argument("-A", "--atac", type=str, required=False, help='ATAC-seq ArchR project folder')
parser.add_argument("-R", "--rna", type=str, required=False, help='RNA-seq scanpy H5 file')
parser.add_argument("-O", "--pseudobulk_outdir", type=str, required=True, help='File output directory for running pseudobulk script')
parser.add_argument("-C", "--covars_atac", nargs='*', type=str, required=False, help='ATAC Column names by individuals and cell calls (ex. -C SNG.BEST.GUESS.ATAC leiden)')
parser.add_argument("-S", "--covars_rna", nargs='*', type=str, required=False, help='RNA Column names by individuals and cell calls (ex. -C SNG.BEST.GUESS.ATAC leiden)')
parser.add_argument("-G", "--rna_cell", type=str, required=False, help='metadata column name for cell type rna')

args = parser.parse_args()
try:
    os.stat(args.pseudobulk_outdir)
except:
    os.mkdir(args.pseudobulk_outdir)
print(str("out directory: "+str(args.pseudobulk_outdir)))


##############
##RNA ANALYSES
##############
if args.rna != None:
    print(str("rna file: "+str(args.rna)))
    #print(str("rna donor column: "+str(args.donor_col)))
    print(str("columns rna: "+str(args.covars_rna)))
    covars=args.covars_rna
    exp_all=sc.read_h5ad(args.rna)
    print(exp_all.obs.columns)
    print(exp_all.raw)
    exp=exp_all.raw.to_adata()
    exp_all=0
    get_tupes_exp = exp.obs[covars][~exp.obs[covars].duplicated()].itertuples
    df_exp = pd.DataFrame(np.nan,index=pd.MultiIndex.from_tuples([i[1:] for i in get_tupes_exp()], names=covars),columns=exp.var_names)
    df_cellcount = pd.DataFrame(columns=list(set(exp.obs[covars[1]])),index=list(set(exp.obs[covars[0]])))
    #df_cell_counts=pd.DataFrame(np.nan,index=pd.MultiIndex.from_tuples([i[1:] for i in get_tupes_exp()], names=covars),columns=exp.var_names)
    print('RNA matrix set up complete')

    for bc, indv, cluster in list(get_tupes_exp()):
        #print(str(bc+" "+indv+" "+cluster))
        subset_bool_exp = ((exp.obs[covars[0]] == indv) & (exp.obs[covars[1]] == cluster))
        view_exp = exp[subset_bool_exp]
        num_cells_exp = view_exp.shape[0]
        arr_exp = view_exp.X.sum(0).A.flatten()
        df_exp.loc[(indv, cluster)] = arr_exp
        df_cellcount.loc[(indv,cluster)]=num_cells_exp
        # test code
        print(df_exp)
        print(view_exp)
        print(num_cells_exp)
        #df_cell_counts.loc[(indv, cluster)] =num_cells_exp

        #df_exp.loc[(indv, cluster)] = np.log1p(arr_exp/num_cells_exp)
    print("done creating RNA pseudobulk matrix")
    print(df_exp)

    #save dataframes of metrics
    donor_counts=df_exp.sum(axis=1).unstack(level=0).T
    donor_cell_counts=df_exp.count(axis=1).unstack(level=0).T
    print(donor_cell_counts.head())
    #"cell_group"
    gene_counts=(df_exp.groupby(by=[args.rna_cell]).sum(0)).T
    gene_counts.index.name = 'genes'

    df_cellcount=df_cellcount.fillna(0)


    #save outputs
    donor_counts.to_csv(str(args.pseudobulk_outdir)+"/donor_counts_exp.txt",sep="\t")
    df_cellcount.to_csv(str(args.pseudobulk_outdir)+"/donor_cell_counts_exp.txt",sep="\t")
    gene_counts.to_csv(str(args.pseudobulk_outdir)+"/gene_counts_exp.txt",sep="\t")

    cell_types = np.unique(df_exp.index.get_level_values(1))
    for ct in cell_types:
        print(ct)
        exp_mat = df_exp.loc[(slice(None), ct), :]
        exp_mat.index = exp_mat.index.droplevel(1)
        exp_mat = exp_mat.T
        indv_list = exp_mat.columns.dropna()
        exp_mat = exp_mat[indv_list]
        print(exp_mat.shape)
        temp_expname = str(args.pseudobulk_outdir)+"/"+str(ct)+".exp"
        exp_mat.to_csv(temp_expname, sep="\t")  # write expression file for RNA-seq data
    print("done writing RNA matrices")


##################
##ATAC ANALYSIS ##
##################

if args.atac != None:
    print(str("atac file: "+str(args.atac)))
    print(str("columns atac: "+str(args.covars_atac)))
    covars=args.covars_atac
    atac=sc.read_h5ad(args.atac)
    get_tupes_atac = atac.obs[covars][~atac.obs[covars].duplicated()].itertuples
    df_atac = pd.DataFrame(np.nan,index=pd.MultiIndex.from_tuples([i[1:] for i in get_tupes_atac()], names=covars),columns=atac.var_names)
    print('ATAC matrix set up complete')

    for bc, indv, cluster in list(get_tupes_atac()):
        #print(str(bc+" "+indv+" "+cluster))
        subset_bool_atac = ((atac.obs[covars[0]] == indv) & (atac.obs[covars[1]] == cluster))
        view_atac = atac[subset_bool_atac]
        num_cells_atac = view_atac.shape[0]
        arr_atac = view_atac.X.sum(0).A.flatten()
        df_atac.loc[(indv, cluster)] = arr_atac
        #df_atac.loc[(indv, cluster)] = np.log1p(arr_atac/num_cells_atac)
    print("done creating ATAC pseudobulk matrix")
    print(df_atac)

    #save dataframes of metrics
    donor_counts=df_atac.sum(axis=1).unstack(level=0).T
    gene_counts=(df_atac.groupby(by=[args.rna_cell]).sum(0)).T
    gene_counts.index.name = 'genes'


    #save outputs
    donor_counts.to_csv(str(args.pseudobulk_outdir)+"/donor_counts_atac.txt",sep="\t")
    gene_counts.to_csv(str(args.pseudobulk_outdir)+"/gene_counts_atac.txt",sep="\t")

    #save outputs
    cell_types = np.unique(df_atac.index.get_level_values(1))
    for ct in cell_types:
        print(ct)
        atac_mat = df_atac.loc[(slice(None), ct), :]
        atac_mat.index = atac_mat.index.droplevel(1)
        atac_mat = atac_mat.T
        indv_list = atac_mat.columns.dropna()
        atac_mat = atac_mat[indv_list]
        print(atac_mat.shape)
        temp_atacname = str(args.pseudobulk_outdir)+"/"+str(ct)+".atac"
        atac_mat.to_csv(temp_atacname, sep="\t")  # write expression file for RNA-seq data
    print("done writing ATAC matrices")



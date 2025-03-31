
#####CONCATENATE PEAK FILES INTO h5

###LOAD LIBS #####

import argparse
import anndata
import scanpy as sc
import sys,os
import pandas as pd
import numpy as np

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100) 
import random
random.seed(5)


########################
######FUCNTIONS#########
########################

def make_adata_from_mtx(mtx_file,var_file,obs_file):
    x=sc.read_mtx(mtx_file).T
    #x=mtx_file.T
    var=pd.DataFrame(pd.read_csv(var_file,sep='\t',header=None))
    #print(var.head())
    # chr1_804712_805212
    #peaks=([ str(x+":"+str(y)+"-"+str(z)) for x,y,z in zip(var['seqnames'], var['start'], var['end']) ])
    chr=[loc.split('_')[0] for loc in var[0] ]
    start=[loc.split('_')[1] for loc in var[0] ]
    end=[loc.split('_')[2] for loc in var[0] ]
    peaks=[ str(x+":"+str(y)+"-"+str(z)) for x,y,z in zip(chr, start, end) ]
    peaks[:5]
    var.index=peaks
    #var['peaks']=peaks

    obs=pd.DataFrame(pd.read_csv(obs_file,sep='\t',header=None))
    print('obs')
    print(obs)
    obs.index=obs[0].tolist()
    #obs['cells']=obs.index
    #print(var.shape)
    #print(obs.shape)
    #print(x.shape)
    x.var=var
    x.obs=obs
    x.var['peaks']=x.var.index
    x.obs.columns=['cell_barcode']
    return(x)
    
def concat_adatas(adatas):
    #adata = adatas[0].concatenate(adatas[1:])
    adata=anndata.concat(adatas,join='outer')
    return(adata)


########################
###read inputs #########
########################
##CMD: python combine_pools_peak_mtx.py -D /data/processed/aggregate_phase2_ICC_L1_peaks/ -O /data/processed/pseudobulk/atac/ -A aggregate_phase2_ICC_L1_peaks_counts_mtx  -P  
print('start')
parser = argparse.ArgumentParser(description="This pseudobulk script accepts ATAC-seq, RNA-seq and/or protein sequencing data")
parser.add_argument("-D", "--dir", type=str, required=False, help='directory conatining mtx files, files should be saved as ${pool_prefix}_obs.tsv,${pool_prefix}_var.tsv, ${pool_prefix}.mtx')
parser.add_argument("-O", "--outdir", type=str, required=True, help='File output directory')
parser.add_argument("-A", "--adataName", type=str, required=True, help='name to save adata object')
parser.add_argument("-P", "--poolPrefixList", nargs='*', type=str, required=False, help='list of pools to concat')
args = parser.parse_args()

try:
    os.stat(args.outdir)
except:
    os.mkdir(args.outdir)

print(str("input directory: "+str(args.dir)))
print(str("out directory: "+str(args.outdir)))
print(str("new file name: "+str(args.adataName)))
print(str("List of pools to merge: "+str(args.poolPrefixList)))

########################
### process data #########
########################
# create adata from each pool mtx files and add to list
adatas=[]
for p in args.poolPrefixList:
    print(str(str(args.dir)+"/"+p+"/peak_clusters.tsv"))
	#adatas.append(make_adata_from_mtx(str(str(args.dir)+"/"+p+".mtx"),str(str(args.dir)+"/"+p+"_var.tsv"),str(str(args.dir)+"/"+p+"_obs.tsv")))
    adatas.append(make_adata_from_mtx(str(str(args.dir)+"/"+p+"/peak_coords.tsv"),str(str(args.dir)+"/"+p+"/peak_clusters.tsv"),str(str(args.dir)+"/"+p+"/peak_barcodes.tsv")))

#peak_barcodes.tsv  peak_coords.tsv XHLT1-POOL-XX1-SC1-Doublet-Summary.rds  peak_clusters.tsv
print('adata list')
print(adatas[1].obs.head())
print(adatas[1].var.head())
print(adatas)
# concatenate all adatas
adata=concat_adatas(adatas)
#ensure type casts are correct
adata.obs.columns = adata.obs.columns.astype(str)
adata.var.columns = adata.var.columns.astype(str)

print('out file')
print(adata)
print(adata.obs.head())
print(adata.var.head())

# save file
adata.write(str(args.outdir+args.adataName+".h5"))

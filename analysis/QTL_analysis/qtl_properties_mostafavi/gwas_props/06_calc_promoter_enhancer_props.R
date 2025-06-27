
## Script to compute promoter-enhancer enrichments of GWAS SNPs.

# This script, as provided, is run on GWAS SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.
set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

args <- commandArgs(TRUE)
gwasfile <- args[1]
outfile <- args[2]

data_dir="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi/"
gwasfile=paste0(data_dir, gwasfile)
infofile=paste0(data_dir, "gwas_eqtl/snp_annotations/promoter_enhancer.counts")
outfile=paste0(data_dir, outfile)

d_info=fread(infofile)

d_gwas=fread(gwasfile)
d_gwas_annots=left_join(d_gwas,d_info,by="SNP")

if ("Pval" %in% colnames(d_gwas_annots)) {
    p1=sapply((d_gwas_annots %>% select(-pheno,-SNP,-Pval)), function(x) length(x[x>0])/length(x)) 
} else {
    p1=sapply((d_gwas_annots %>% select(-SNP)), function(x) length(x[x>0])/length(x)) 
}
p0=sapply((d_info %>% select(-SNP)), function(x) length(x[x>0])/length(x) ) 

enrichments=data.frame(t(p1/p0))

write.table(enrichments,file=outfile,quote=F,sep=",",row.names=F)

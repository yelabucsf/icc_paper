
## Script to compute promoter-enhancer enrichments of eQTL SNPs.

# This script, as provided, is run on eQTL SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

args <- commandArgs(TRUE)
eqtlfile <- args[1]
outfile <- args[2]

data_dir="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi/"
eqtlfile=paste0(data_dir, eqtlfile)
infofile=paste0(data_dir, "gwas_eqtl/snp_annotations/promoter_enhancer.counts")
outfile=paste0(data_dir, outfile)

d_info=fread(infofile)

d_eqtl=fread(eqtlfile)
if ("Pval" %in% colnames(d_eqtl)) {
    d_eqtl=d_eqtl[d_eqtl$Pval<(5e-8),] #restrict to top eQTLs using a p-value cutoff matching GWAS
}

d_eqtl_annots=left_join(d_eqtl,d_info,by=c("SNP"))

if ("Pval" %in% colnames(d_eqtl_annots)) {
    p1=sapply((d_eqtl_annots %>% select(-Tissue,-SNP,-Pval,-GeneSymbol)), function(x) length(x[x>0])/length(x) ) 
} else {
    p1=sapply((d_eqtl_annots %>% select(-SNP)), function(x) length(x[x>0])/length(x) ) 
}
# p1=sapply((d_eqtl_annots %>% select(-Tissue,-SNP,-Pval,-GeneSymbol)), function(x) length(x[x>0])/length(x) ) 
p0=sapply((d_info %>% select(-SNP)), function(x) length(x[x>0])/length(x) ) 

enrichments=data.frame(t(p1/p0))

write.table(enrichments,file=outfile,quote=F,sep=",",row.names=F)

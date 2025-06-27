## Script to filter caQTLs -- only difference from eQTLs is that we don't filter by protein-coding genes
args <- commandArgs(TRUE)
level <- args[1]

#file concatenating the output of codes in "eqtl_process" across all tissue-eGene pairs
data_dir="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi/"
eqtlfiles = list.files(path = paste0(data_dir, "atac/", level), 
                       recursive = TRUE, pattern = "\\.clumped\\.hg19\\.snps$", full.names = TRUE)
snpfile=paste0(data_dir, "gwas_eqtl/snp_annotations/filter_snps.txt")
outfile1=paste0(data_dir, "atac/", level, "/filter_qtls.hg19.assoc")

library(tidyverse)
library(data.table)
set.seed(0)

d_snp=fread(snpfile)

# read in and concatenate all eQTL files
d_eqtl=data.frame()
for (f in eqtlfiles) {
    d_eqtl=rbind(d_eqtl,fread(f,header=F))
}

colnames(d_eqtl)=c("Tissue","Gene","SNP","Pval","rank")
d_eqtl$SNP=as.character(d_eqtl$SNP)
d_eqtl <- (d_eqtl) %>% separate(Gene,c("GeneSymbol","nn"),sep="\\.",remove=F) %>% select(-nn)
d_eqtl <- d_eqtl %>% select(Tissue,GeneSymbol,SNP,Pval)

df=d_eqtl[(d_eqtl$SNP %in% d_snp$SNP),]

#######

write.table(df,file=outfile1,quote=F,sep="\t",row.names=F)

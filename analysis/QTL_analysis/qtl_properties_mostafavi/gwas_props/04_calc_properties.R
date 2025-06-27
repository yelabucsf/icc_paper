
## Script to compute basic properties of GWAS SNPs.

# This script, as provided, is run on GWAS SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.

# setwd("/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi/gwas_eqtl")
set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

args <- commandArgs(TRUE)
gwasfile <- args[1]
outfile <- args[2]

data_dir="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi/"
gwasfile=paste0(data_dir, gwasfile)
infofile=paste0(data_dir, "gwas_eqtl/snp_annotations/filter_snps.txt")
genefile=paste0(data_dir, "gwas_eqtl/gene_annotations/all_annots_pc_genes.txt")
outfile=paste0(data_dir, outfile)

d_info=fread(infofile) #SNP features
d_gene=fread(genefile) #genic features
d_gene$gene=d_gene$hgnc_id

#bin continuous genic features
d_gene$LOEUF_cat=ntile(d_gene$LOEUF,5)
d_gene$Road_length_cat=ntile(d_gene$Roadmap_length_per_type,5)
d_gene$Road_count_cat=ntile(d_gene$Roadmap_count,5)
d_gene$ABC_length_cat=ntile(d_gene$ABC_length_per_type,5)
d_gene$ABC_count_cat=ntile(d_gene$ABC_count,5)
d_gene$TSS_cat=ntile(d_gene$promoter_count,5)

d_all=left_join(d_info,(d_gene %>% select(-TSSD)),by="gene")
d_all$d_tss_signed=(d_all$d_tss)
d_all$d_tss=abs(d_all$d_tss)

#bin distance to TSS, MAF, LD score, and gene density for regression models
d_all$d_tss_cat=ntile(d_all$d_tss,n = 20)
d_all$MAF_cat=ntile(d_all$MAF,n = 20)
d_all$L2_cat=ntile(d_all$L2,n = 20)
d_all$TSSD_cat=ntile(d_all$TSSD,n = 20)

d_gwas=fread(gwasfile)
d_gwas_annots=left_join(d_gwas,d_all,by="SNP")

########

#compute props
d_test=d_gwas_annots 

prop_HI=mean(d_test$HI,na.rm=T)
prop_TF=mean(d_test$TF,na.rm=T)
prop_connect_top_decile=nrow(d_test[d_test$connect_decile==10,])/nrow(d_test)
prop_connect_top_half=nrow(d_test[d_test$connect_decile>=5,])/nrow(d_test)
mean_promoter_count=mean(d_test$promoter_count,na.rm=T)
mean_connectedness=mean(d_test$connectedness,na.rm=T)

d_x1=d_test[d_test$connectedness==1,] %>% group_by(connect_quantile) %>% count()
prop_connect_qs=t(d_x1$n/nrow(d_test))
colnames(prop_connect_qs)=paste0("prop_connect_q",as.character(1:5))

d_x1=d_test[d_test$PPI_degree_cat==1,] %>% group_by(PPI_degree_quantile) %>% count()
prop_PPI_qs=t(d_x1$n/nrow(d_test))
colnames(prop_PPI_qs)=paste0("prop_PPI_q",as.character(1:5))

d_x1=d_test[!is.na(d_test$LOEUF_cat),] %>% group_by(LOEUF_cat) %>% count()
prop_LOEUF_qs=t(d_x1$n/nrow(d_test))
colnames(prop_LOEUF_qs)=paste0("prop_LOEUF_q",as.character(1:5))

d_x1=d_test %>% group_by(ABC_length_cat) %>% count()
prop_ABC_length_qs=t(d_x1$n/nrow(d_test))
colnames(prop_ABC_length_qs)=paste0("prop_ABC_length_qs",as.character(1:5))

d_x1=d_test %>% group_by(ABC_count_cat) %>% count()
prop_ABC_count_qs=t(d_x1$n/nrow(d_test))
colnames(prop_ABC_count_qs)=paste0("prop_ABC_count_qs",as.character(1:5))

d_x1=d_test %>% group_by(Road_length_cat) %>% count()
prop_Roadmap_length_qs=t(d_x1$n/nrow(d_test))
colnames(prop_Roadmap_length_qs)=paste0("prop_Roadmap_length_qs",as.character(1:5))

d_x1=d_test %>% group_by(Road_count_cat) %>% count()
prop_Roadmap_count_qs=t(d_x1$n/nrow(d_test))
colnames(prop_Roadmap_count_qs)=paste0("prop_Roadmap_count_qs",as.character(1:5))

d_out1=data.frame(prop_HI,prop_TF,prop_connect_top_decile,prop_connect_top_half,mean_promoter_count,mean_connectedness,prop_connect_qs,prop_PPI_qs,prop_LOEUF_qs,prop_Roadmap_length_qs,prop_Roadmap_count_qs,prop_ABC_length_qs,prop_ABC_count_qs)
write.table(d_out1,file=outfile,quote=F,sep=",",row.names=F)


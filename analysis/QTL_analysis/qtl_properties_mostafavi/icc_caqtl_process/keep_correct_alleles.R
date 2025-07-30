args <- commandArgs(TRUE)
maffile <- args[1]
hg38_snps <- args[2]
infile1 <- args[3]
infile2 <- args[4]
infile3 <- args[5]

library(tidyverse)
library(data.table)

d_var=read.csv(infile1,sep="",header=F)
d_pos1=read.csv(infile2,sep="",header=F)
d_pos2=read.csv(infile3,sep="",header=F)
print("done reading first 3 input files")

colnames(d_var)=c("var","gene","P","snp")
d_var <- separate(d_var, var, c("chr","pos_38","A1","A2","build"), sep = "_")
print("done separating d_var")

d_pos1=d_pos1[c("V1","V3","V4")]; colnames(d_pos1)=c("chr","pos_38","snp")
d_pos2=d_pos2[c("V1","V3","V4")]; colnames(d_pos2)=c("chr","pos_19","snp")
d_pos=merge(d_pos1,d_pos2,by=c("chr","snp"))
print("done merging d_pos1 and d_pos2")
rm(d_pos1) 
rm(d_pos2) 

df=merge(d_var,d_pos,by=c("chr","pos_38","snp"))
rm(d_var)
rm(d_pos)
print("done merging d_var and d_pos")

df$chr_num=gsub("[a-z]","", df$chr)
df$pos_38_num=as.numeric(as.character(df$pos_38))
df$pos_19_num=as.numeric(as.character(df$pos_19))
df$chr_num=as.numeric(as.character(df$chr_num))
df <- df[(df$chr_num>=1 & df$chr_num<=22),]
print("done filtering to autosomal chromosomes")

cols <- c('chr_num','pos_38_num','A1','A2')
df <- unite(df,"allele1_hg38", cols, sep=":", na.rm = TRUE, remove = FALSE)
cols <- c('chr_num','pos_38_num','A2','A1')
df <- unite(df,"allele2_hg38", cols, sep=":", na.rm = TRUE, remove = FALSE)

d_hg38_snps=fread(hg38_snps,header = F)
print("done reading hg38 snps file")

dg1=df[df$allele1_hg38 %in% d_hg38_snps$V1,]
dg2=df[df$allele2_hg38 %in% d_hg38_snps$V1,]
rm(d_hg38_snps)
rm(df)

dg1$allele_hg38 = dg1$allele1_hg38
dg2$allele_hg38 = dg2$allele2_hg38
dg=rbind(dg1,dg2)
rm(dg1)
rm(dg2)
print("done merging dg1 and dg2")

cols <- c('chr_num','pos_19_num','A1','A2')
dg <- unite(dg,"allele1_hg19", cols, sep=":", na.rm = TRUE, remove = FALSE)
cols <- c('chr_num','pos_19_num','A2','A1')
dg <- unite(dg,"allele2_hg19", cols, sep=":", na.rm = TRUE, remove = FALSE)
print("done merging allele1 and allele2")

d_maf=fread(maffile,header = F)
print("done reading maf file")

dh1=dg[dg$allele1_hg19 %in% d_maf$V1,]
dh2=dg[dg$allele2_hg19 %in% d_maf$V1,]
rm(d_maf)
rm(dg)
print("done merging dh1 and dh2")

dh1=dh1[c("chr_num","gene","allele_hg38","P")]; colnames(dh1)=c("CHR","GENE","SNP","P")
dh2=dh2[c("chr_num","gene","allele_hg38","P")]; colnames(dh2)=c("CHR","GENE","SNP","P")
dh=rbind(dh1,dh2)
dh=dh[order(dh$GENE),]
rm(dh1)
rm(dh2)
print("done ordering dh")

write.table(dh,file=infile1,quote = FALSE,sep=" ",row.names = F, col.names=F)

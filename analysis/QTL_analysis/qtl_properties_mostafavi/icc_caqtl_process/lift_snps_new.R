args <- commandArgs(TRUE)
maffile <- args[1]
infile1 <- args[2]
infile2 <- args[3]
infile3 <- args[4]
outfile <- args[5]

library(tidyverse)
library(data.table)

d_var=read.csv(infile1,sep="",header=F)
d_pos1=read.csv(infile2,sep="",header=F)
d_pos2=read.csv(infile3,sep="",header=F)
d_maf=fread(maffile,header = F)

colnames(d_var)=c("cell_type","gene","snp","P", "num")
d_var <- separate(d_var, snp, c("chr","pos_38","A1","A2"), sep = ":")

d_pos1=d_pos1[c("V1","V3","V4")]; colnames(d_pos1)=c("chr","pos_38","snp")
d_pos2=d_pos2[c("V1","V3","V4")]; colnames(d_pos2)=c("chr","pos_19","snp")
d_pos=merge(d_pos1,d_pos2,by=c("chr","snp"))
d_pos$chr=gsub("[a-z]","", d_pos$chr)

df=merge(d_var,d_pos,by=c("chr","pos_38"))
df$chr_num=gsub("[a-z]","", df$chr)
df$pos_19_num=as.numeric(as.character(df$pos_19))
df$chr_num=as.numeric(as.character(df$chr_num))
df <- df[(df$chr_num>=1 & df$chr_num<=22),]

cols <- c('chr_num','pos_19_num','A1','A2')
df <- unite(df,"allele1", cols, sep=":", na.rm = TRUE, remove = FALSE)

cols <- c('chr_num','pos_19_num','A2','A1')
df <- unite(df,"allele2", cols, sep=":", na.rm = TRUE, remove = FALSE)

dg1=df[df$allele1 %in% d_maf$V1,]
dg2=df[df$allele2 %in% d_maf$V1,]

dg1=dg1[c("cell_type","gene","allele1","P", "num")]
colnames(dg1)=c("cell_type","gene","snp","P", "num")
dg2=dg2[c("cell_type","gene","allele2","P", "num")]
colnames(dg2)=c("cell_type","gene","snp","P", "num")
dg=rbind(dg1,dg2)
dg=unique(dg[order(dg$snp),])

write.table(dg,file=outfile,quote = FALSE,sep=" ",row.names = F, col.names=F)

library(data.table)
library(stringr)
library(ggplot2)
library(optparse)

#######################################################
### Author Gracie 2021
### 
### This script takes in the .evec output of make_genosPC (smartpca)
### expression/peak PCs the output of make_pheno_peakPC
### and the output of extract_clin_covs.R to merge into one covariates matrix
#################################################

## CMD: Rscript combine_covariates.R -G genopc.evec -P phenpc.txt -M clincovs.txt -C ${cell_type} -P prefix -O output
###############
### READ IN ARGS
################

##accept commandline args
option_list = list(
make_option(c("-G", "--genoPC"), action="store", default=NULL, type='character',
              help="evec file output of smartpca"),
make_option(c("-P", "--phenoPC"), action="store", default=NULL, type='character',
              help="PCs for phenotype, ex: expression, peaks"),
make_option(c("-M", "--metaData"), action="store", default=NULL, type='character',
              help="donor clinical data"),
make_option(c("-C", "--cell"), action="store", default=NULL, type='character',
              help="cell type"),
make_option(c("-P", "--prefix"), action="store", default=NULL, type='character',
              help="prefix to use"),
##make_option(c("-A", "--ngenoPC"), action="store", default=NULL, type='character',
##              help="number of genotype pcs to use"),
##make_option(c("-B", "--nphenoPC"), action="store", default=NULL, type='character',
##              help="number of phenotype to use"),		   
make_option(c("-O", "--out_dir"), action="store", default=NULL, type='character',
              help="outdirectory"))

opt = parse_args(OptionParser(option_list=option_list))

gPC_f=opt$genoPC
pPC_f=opt$phenoPC
meta_f=opt$metaData
phen_f=opt
cell=opt$cell
pref=opt$prefix
##ngenoPC=opt$ngenoPC
##nphenoPC=opt$nphenoPC
outdir=opt$out_dir


##############
## FUNCTIONS #
##############


clean_genoPC<-function(file,out){
    genoPC=fread(file)
    cols=t(genoPC)[1,]
    #rows=c('geno_pc1', 'geno_pc2', 'geno_pc3', 'geno_pc4', 'geno_pc5', 'geno_pc6', 'geno_pc7', 'geno_pc8', 'geno_pc9', 'geno_pc10')

    genoPC2=head(tail(as.data.frame(t(genoPC)),-1),-1)
    npc=length(rownames(genoPC2))
    rows=paste0('geno_pc',seq(1,npc))
    colnames(genoPC2)=cols
    rownames(genoPC2)=rows
    print(head(genoPC2))
    write.table(genoPC2,file=out,sep='\t',quote=FALSE,row.names = TRUE)
    
}

merge_covariates<-function(meta_f,geno_f,expr_f,cell,outdir,prefix){
    meta=as.data.frame(fread(meta_f))
    geno=as.data.frame(fread(geno_f))
    expr=as.data.frame(fread(expr_f))
    g_r=(geno$V1)
    e_r=(expr$V1)
    g_inds=colnames(geno)[-1]
    e_inds=colnames(expr)[-1]
    geno=as.data.frame(t(geno))[-1,]
    colnames(geno)=g_r
    expr=as.data.frame(t(expr))[-1,]
    colnames(expr)=e_r
    #print(geno)
    geno$geno_ID=g_inds
    expr$geno_ID=e_inds
    
    df=merge(meta,geno,by=c("geno_ID"))
    df=merge(df,expr,by=c("geno_ID"))

    inds=df$geno_ID
    df=(t(df))[-1,]
    colnames(df)=inds
    df=df[,e_inds]
    print(df)
    write.table(df,paste0(outdir,"/",prefix,"_",cell,".covs"),quote=FALSE, col.names=NA,sep='\t')
    
    
}


##############
## MAIN ######
##############

clean_genoPC(gPC_f,paste0(outdir,"/",pref,"_geno.cov"))
merge_covariates(m,gPC_f,pPC_f,cell,outdir,pref)

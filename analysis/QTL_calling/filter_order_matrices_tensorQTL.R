suppressPackageStartupMessages(require(optparse))
library(data.table)

##CMD: Rscript filter_order_matrices_tensorQTL.R -X expresion.expr -V cov_file.txt  -P pref  -O ${out_dir} 
 

## Settings

##accept commandline args
option_list = list(
##make_option(c("-G", "--genos"), action="store", default=NULL, type='character',
##              help="genotype file"),
make_option(c("-X", "--gex"), action="store", default=NULL, type='character',
              help="gene expression file"),
make_option(c("-V", "--covariate"), action="store", default=NULL, type='character',
              help="covariate file"),
make_option(c("-P", "--prefix"), action="store", default=NULL, type='character',
              help="prefix to use for outputs"),
##make_option(c("-C", "--chromosome"), action="store", default=NULL, type='character',
##              help="chromosome number"),
make_option(c("-I", "--inds"), action="store", default="none", type='character',
              help="file contianing list of individuals"),
make_option(c("-F", "--percent_inds"), action="store", default="0.8", type='double',
              help="percent inds with expression"),
make_option(c("-L", "--gene_locs"), action="store", default="/data/references/gene_locs.txt", type='character',
              help="percent inds with expression"),
make_option(c("-O", "--outdir"), action="store", default=NULL, type='character',
              help="outfile prefix"))


opt = parse_args(OptionParser(option_list=option_list))

##read in files
## These are formated for V1 being the first column and not a donor name
#phenos=as.data.frame(as.matrix(fread(opt$gex),rownames=1))
phenos=as.data.frame(as.matrix(fread(opt$gex)))
#covs=as.data.frame(as.matrix(fread(opt$covariate,header=TRUE),rownames=1))
covs=as.data.frame(as.matrix(fread(opt$covariate,header=TRUE)))
rownames(covs)=covs$PersonID

#get set of donors to keep in the analyses

print(head(colnames(phenos)))
print(head(colnames(covs)))
donors=intersect(colnames(phenos[-1]),colnames(covs[-1]))

print(length(donors))
print(setdiff(colnames(phenos),colnames(covs)))
print(setdiff(colnames(covs),colnames(phenos)))


if(opt$inds!="none"){
	print("filter")
	inds=as.data.frame(fread(opt$inds, header=FALSE),header=FALSE)$V1
    print(length(inds))
	donors=sort(intersect(donors,inds))

}
print(paste("number donors retained: ",length(donors)))

print("order individuals")
#print(head(phenos))
#print(donors)
phenos2=phenos[c("V1",donors)]
min_inds=(dim(phenos2)[2]-1)*opt$percent_inds
phenos3=phenos2[rowSums(phenos2 == 0) <= min_inds, ]
covs2=covs[c("V1",donors)]

# unname V1
#colnames(phenos3)=unname(phenos3[1,])
#colnames(covs2)=unname(covs2[1,])
colnames(phenos3)=c("",colnames(phenos3)[-1])
colnames(covs2)=c("",colnames(covs2)[-1])


print(dim(phenos3))
print(dim(covs2))

##print(head(phenos3))
##print(head(covs2))

##save outputs
dir.create(file.path(dirname(opt$outdir), basename(opt$outdir)), showWarnings = FALSE)

write.table(covs2, file=paste0(opt$outdir,"/",opt$prefix,".cov"),quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")
write.table(phenos3, file=paste0(opt$outdir,"/",opt$prefix,".phenos"),quote=FALSE,row.names=FALSE,sep="\t")

## generate pheno bedfile
print(opt$gene_locs)
gene_locs=as.data.frame(fread(opt$gene_locs))
colnames(gene_locs)=c('gene_id','#chr','start','end')
gene_locs=gene_locs[,c('#chr','start','end','gene_id')]
print(head(gene_locs))
colnames(phenos3)=c("gene_id",colnames(phenos3)[-1])
bed_phenos=merge(gene_locs,phenos3, by='gene_id',how='right')
bed_phenos[['#chr']]=gsub('chr','',bed_phenos[['#chr']])
bed_phenos=bed_phenos[,c('#chr','start','end','gene_id',donors)]

print(head(bed_phenos))
write.table(bed_phenos, file=paste0(opt$outdir,"/",opt$prefix,"_phenos_tensorqtl.bed"),quote=FALSE, col.names=TRUE,row.names=FALSE,sep='\t')



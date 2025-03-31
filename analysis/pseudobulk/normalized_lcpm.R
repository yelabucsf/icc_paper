library(edgeR)
library("optparse")
library(data.table)
library(stringr)


##CMD Rscript normalized_lcpm.R -E -G -D -C B -O test
option_list = list(
  make_option(c("-E", "--expr"), action="store", default=NA, type='character',
              help="pseudobulk count matrix"),
  make_option(c("-G", "--gene_count"), action="store", default=NA, type='character',
              help="count matrix of genes by cell types"),
  make_option(c("-D", "--donor_count"), action="store", default=NA, type='character',
              help="count matrix of donor by cell types"),
  make_option(c("-C", "--celltype"), action="store", default=NA, type='character',
              help="cell type to run (must match columns of donor/gene count matrices"),		   
  make_option(c("-O", "--outdir"), action="store", default=NA, type='character',
              help="output directory"),
  make_option(c("-M", "--mode"), action="store", default='exp', type='character',
              help="define modality used, atac or exp, expression is default, this just changes suffix"),
  make_option(c("-N", "--ncell"), action="store", default='0', type='integer',
              help="minimum number of cells per donor per group to include"),
  make_option(c("-F", "--countfile"), action="store", default=NA, type='character',
              help="file of cell counts per donor")
)
opt = parse_args(OptionParser(option_list=option_list))

print(paste("expression matrix:",opt$expr))
print(paste("gene matrix:",opt$gene_count))
print(paste("donor matrix:",opt$donor_count))
print(paste("cell type:",opt$celltype))
print(paste("output directory:",opt$outdir))
print(paste("mode:",opt$mode))
print(paste("min cells:", opt$ncell))
print(paste("file of cell counts:", opt$countfile))

##read in data
exprMat=as.matrix(fread(opt$expr),rownames=1)
geneMat=as.data.frame(fread(opt$gene_count))
rownames(geneMat)=geneMat$genes
donorMat=as.data.frame(fread(opt$donor_count))
rownames(donorMat)=donorMat$best_sample
cellcountMat=as.data.frame(fread(opt$countfile))
colnames(cellcountMat)=c('donor',colnames(cellcountMat)[-1])

#replace NAs with 0
exprMat[is.na(exprMat)] <- 0
geneMat[is.na(geneMat)] <- 0
donorMat[is.na(donorMat)] <- 0
cellcountMat[is.na(cellcountMat)] <- 0

print(head(geneMat))
print(head(donorMat))
print((exprMat[1:6,1:6]))

#filter donors with too few cells
opt$celltype
keep_donor=cellcountMat[which(cellcountMat[opt$celltype]>=opt$ncell),]$donor

##verify ordering of individuals 
all_donors=intersect(rownames(donorMat),colnames(exprMat))
all_donors=intersect(all_donors,keep_donor)
donorMat <- donorMat[which(rownames(donorMat) %in% all_donors),]
donorMat <- donorMat[order(all_donors),]

exprMat <- exprMat[,all_donors ]
exprMat <- exprMat[,order(all_donors) ]


##filter for cell type
sub("_", " ", opt$celltype)
geneMat=geneMat[c(sub("_", " ", opt$celltype))]
donorMat=donorMat[c(sub("_", " ", opt$celltype))]


print("lib")
print(head(donorMat[,1]))
print("samp")
print(head(geneMat))

x=DGEList(counts=exprMat,lib.size=as.numeric(donorMat[,1]),samples=rownames(donorMat),genes=geneMat)

##filter genes without expression
table(rowSums(x$counts==0)==length(rownames(donorMat)))
#keep.exprs <- filterByExpr(x)

###MINIMUM UMI COUNTS HARDCODED HERE
keep.exprs <- rowSums(x$counts) >= 15
x <- x[keep.exprs, keep.lib.sizes=FALSE]
x <- calcNormFactors(x, method = "TMM")

print(x$lib.size)
print(x$norm.factors)

lcpm <- edgeR::cpm(x, log=TRUE)
s_lcpm <- t(scale(t(lcpm)))
n <- ncol(lcpm)
zvalues <- qnorm(ppoints(n))
z <- lcpm
for (i in 1:nrow(z)) z[i,] <- zvalues[order(order(z[i,]))]

if(opt$mode=='exp'){
    write.table(lcpm,file=paste0(opt$outdir,"/",opt$celltype,"_expr.lcpm"),quote=FALSE)
    write.table(s_lcpm,file=paste0(opt$outdir,"/",opt$celltype,"_scaled_expr.lcpm"),quote=FALSE)
    write.table(z,file=paste0(opt$outdir,"/",opt$celltype,"_expr.lcpm.qnorm"),quote=FALSE)
}else{
    write.table(lcpm,file=paste0(opt$outdir,"/",opt$celltype,"_atac.lcpm"),quote=FALSE)
    write.table(s_lcpm,file=paste0(opt$outdir,"/",opt$celltype,"_scaled_atac.lcpm"),quote=FALSE)
    write.table(z,file=paste0(opt$outdir,"/",opt$celltype,"_atac.lcpm.qnorm"),quote=FALSE)

}



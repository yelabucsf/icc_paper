#############################
#### ICC ARCHR BIG MERGES ###
#### Gracie 2022
####
####
#############################

## LOAD PACKAGES
suppressPackageStartupMessages(library(ArchR))
addArchRGenome("hg38")
set.seed(5)
library("optparse")
library(data.table)
library(ggplot2)
library(stringr)

####################################################
##DEFINE POOLS and WELLS
##TODO: CHANGE THIS TO FEED IN LISTS WITH CMD ARGS
####################################################

outdir <- "/data/processed/icc_phase2/"
#dir_free="/data/ambimux_per_well/"
dir_free='/data/'

pools=c('XHLT1-POOL-XX1','XHLT1-POOL-XX2','XHLT1-POOL-XX3','XHLT1-POOL-XX4','XHLT1-POOL-XX5',
'XHLT1-POOL-XX6','XHLT1-POOL-XX7','XHLT1-POOL-XX8','XHLT1-POOL-XX9','XHLT1-POOL-XX10',
'XHLT1-POOL-XX11','XHLT1-POOL-XX12','XHLT1-POOL-XX13','XHLT1-POOL-XX14','XHLT1-POOL-XX15',
'XHLT1-POOL-XX16','XHLT1-POOL-XX19','XHLT1-POOL-XX20','XHLT1-POOL-XX21','XHLT1-POOL-XX22')

samples=c('XHLT1-POOL-XX1-SC1',
'XHLT1-POOL-XX1-SC2',
'XHLT1-POOL-XX1-SC3',
'XHLT1-POOL-XX1-SC4',
'XHLT1-POOL-XX10-SC1',
'XHLT1-POOL-XX10-SC2',
'XHLT1-POOL-XX10-SC3',
'XHLT1-POOL-XX10-SC4',
'XHLT1-POOL-XX11-SC1',
'XHLT1-POOL-XX11-SC2',
'XHLT1-POOL-XX11-SC3',
'XHLT1-POOL-XX11-SC4',
'XHLT1-POOL-XX12-SC1',
'XHLT1-POOL-XX12-SC2',
'XHLT1-POOL-XX12-SC3',
'XHLT1-POOL-XX12-SC4',
'XHLT1-POOL-XX13-SC1',
'XHLT1-POOL-XX13-SC2',
'XHLT1-POOL-XX13-SC3',
'XHLT1-POOL-XX14-SC1',
'XHLT1-POOL-XX14-SC3',
'XHLT1-POOL-XX14-SC4',
'XHLT1-POOL-XX15-SC1',
'XHLT1-POOL-XX15-SC2',
'XHLT1-POOL-XX15-SC3',
'XHLT1-POOL-XX15-SC4',
'XHLT1-POOL-XX16-SC1',
'XHLT1-POOL-XX16-SC2',
'XHLT1-POOL-XX16-SC3',
'XHLT1-POOL-XX16-SC4',
'XHLT1-POOL-XX19-SC1',
'XHLT1-POOL-XX19-SC2',
'XHLT1-POOL-XX19-SC3',
'XHLT1-POOL-XX19-SC4',
'XHLT1-POOL-XX2-SC1',
'XHLT1-POOL-XX2-SC2',
'XHLT1-POOL-XX2-SC3',
'XHLT1-POOL-XX2-SC4',
'XHLT1-POOL-XX20-SC1',
'XHLT1-POOL-XX20-SC2',
'XHLT1-POOL-XX20-SC3',
'XHLT1-POOL-XX21-SC1',
'XHLT1-POOL-XX21-SC2',
'XHLT1-POOL-XX21-SC3',
'XHLT1-POOL-XX21-SC4',
'XHLT1-POOL-XX22-SC1',
'XHLT1-POOL-XX22-SC2',
'XHLT1-POOL-XX22-SC3',
'XHLT1-POOL-XX22-SC4',
'XHLT1-POOL-XX3-SC1',
'XHLT1-POOL-XX3-SC2',
'XHLT1-POOL-XX3-SC3',
'XHLT1-POOL-XX3-SC4',
'XHLT1-POOL-XX4-SC1',
'XHLT1-POOL-XX4-SC2',
'XHLT1-POOL-XX4-SC3',
'XHLT1-POOL-XX4-SC4',
'XHLT1-POOL-XX5-SC1',
'XHLT1-POOL-XX5-SC2',
'XHLT1-POOL-XX5-SC3',
'XHLT1-POOL-XX5-SC4',
'XHLT1-POOL-XX6-SC1',
'XHLT1-POOL-XX6-SC2',
'XHLT1-POOL-XX6-SC3',
'XHLT1-POOL-XX6-SC4',
'XHLT1-POOL-XX7-SC1',
'XHLT1-POOL-XX7-SC2',
'XHLT1-POOL-XX7-SC3',
'XHLT1-POOL-XX7-SC4',
'XHLT1-POOL-XX8-SC1',
'XHLT1-POOL-XX8-SC2',
'XHLT1-POOL-XX8-SC3',
'XHLT1-POOL-XX8-SC4',
'XHLT1-POOL-XX9-SC2',
'XHLT1-POOL-XX9-SC4')

#'XHLT1-POOL-XX8',


###ArrowFiles=c()
###rna_files=c()
###for(x in samples){
###    ArrowFiles=c(ArrowFiles,paste0('/data/ArrowFiles/atac_only/ArrowFiles/',x,'.arrow'))             
###    rna_files=c(rna_files,paste0('/data/cellranger_outs/',x,'/filtered_feature_bc_matrix.h5'))
###}
###print(ArrowFiles)
markerGenes  <- c(
    #"IL7R","CD3D",  #CD4T
    "CD8A",#"CD3D", #CD8T
    "CD14", #"LYZ", #CD14 Monocytes
    "FCGR3A", "MS4A7", #CD16 Monocytes
    "MS4A1", #B cells
    "GNLY", "NKG7", #NK
    "FCER1A", #"CST3", #DC
    "IRF4","PRDM1", #PB
    "CPNE3",#"CAMK2D","DNASE1L3","TACSTD2","CLEC9A","XCR1",#CDC1
    "ALOX5AP",#"SAMSN1","TYROBP","SORL1","CCND2",
    "PRDM1",#,"PTGS2","TAGAP","EREG","INHBA","CD1C","ZEB2"#CDC2
    "IRF7", "IRF8" #DC
    #"PPBP" #MK
    )


#######################################
#### DATA INTAKE, MERGING, AND     ####
####  ORGANIZATION                 ####
#######################################

###READ IN ARROW FILES AND H5
proj <- ArchRProject( ArrowFiles=ArrowFiles, outputDirectory=outdir, copyArrows = FALSE)
#

##rename barcodes to add -1 as the well number so freemux files can be merged

tmp=as.data.frame(t(as.data.frame(strsplit((rownames(proj@cellColData)),"#"),header=FALSE,stringsAsFactors=FALSE)),stringsAsFactors=FALSE)
#num=str_sub(tmp$V1,-1)
#bc=str_sub(tmp$V2,1,17)
#new=paste0(bc,num)
new=tmp$V2
new_idx=paste0(tmp$V1,"#",new)
proj@cellColData$full_name=new_idx
print(head(proj@cellColData))
#read merged freemux files and add a donor and pool column to the files
frees=as.data.frame(fread(paste0(dir_free,'ambimux.summary.txt')))
frees$full_name=paste0(frees$well_id,"#",frees$Barcode)
frees$pool=sub("\\-SC.*", "", frees$well_id)

# add freemux annotations to metadata
proj2=proj
print(proj2)
proj2@cellColData$rows=rownames(proj2@cellColData)
print(head(frees))
print(head(proj2@cellColData))
test=merge(proj2@cellColData,frees,by="full_name",all.x = TRUE)
rownames(test)=test$full_name
proj2@cellColData=test
print(proj2)
proj2 <- proj2[proj2$TSSEnrichment > 6 & proj2$nFrags > 2500 ]
proj2 <- addDoubletScores(proj2)

print(proj2)
print(unique(proj2$best_type))
#Filter for only Freemux called Singlets
idxPass <- which(proj2$best_type == "Singlet")
cellsPass <- proj2$cellNames[idxPass]
proj2=proj2[cellsPass, ]
print(proj2)
#proj2 <- filterDoublets(proj2)

#proj2 <- filterDoublets(proj2)


###################################
#### DIMENSIONALITY REDUCTION  ####
#### AND CLUSTERING            ####
###################################

#LSI-ATAC
proj2 <- addIterativeLSI(
    ArchRProj = proj2,
    projectCellsPre = TRUE,
    sampleCellsFinal=400000,
    clusterParams = list(
      resolution = 0.2,
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix",
    depthCol = "nFrags",
    name = "LSI_ATAC"
)


#Combined Dims
#proj2 <- addCombinedDims(proj2, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
#Create UMAPs for ATAC/RNA and COMBINED
proj2 <- addUMAP(proj2, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
## CLUSTER BASED ON COMBINED DIMENSIONS
proj2 <- addClusters(proj2, reducedDims = "LSI_ATAC", name = "Clusters", resolution = 0.4, force = TRUE)

#Plot and save Embedding no BATCH correction
umap_atac_nh <- plotEmbedding(proj2, name = "Clusters", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)

plotPDF(umap_atac_nh,name="ATAC_umap_noHarmony.pdf", ArchRProj = proj2, addDOC = FALSE)

sample_umap_comb_nh <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_ATAC",force=TRUE)
cluster_umap_comb_nh <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_ATAC",force=TRUE)
donor_umap_comb_nh <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "best_sample", embedding = "UMAP_ATAC",force=TRUE)
pool_umap_comb_nh <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "pool", embedding = "UMAP_ATAC",force=TRUE)

plotPDF(sample_umap_comb_nh,name="sample_comb_umap_noHarmony.pdf", ArchRProj = proj2, addDOC = FALSE)
plotPDF(cluster_umap_comb_nh,name="cluster_comb_umap_noHarmony.pdf", ArchRProj = proj2, addDOC = FALSE)
plotPDF(donor_umap_comb_nh,name="donor_comb_umap_noHarmony.pdf", ArchRProj = proj2, addDOC = FALSE)
plotPDF(pool_umap_comb_nh,name="pool_comb_umap_noHarmony.pdf", ArchRProj = proj2, addDOC = FALSE)

###VISUALIZE MARKERS###

genes_umap_comb_nh<- plotEmbedding(
    ArchRProj = proj2,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP_ATAC"
)

plotPDF(plotList=genes_umap_comb_nh,name="genes_comb_umap_noHarmony.pdf", ArchRProj = proj2, addDOC = FALSE)

saveArchRProject(ArchRProj = proj2, outputDirectory = "/data/processed/aggregate_phase2_ICC_noHarmony", load = FALSE)

proj2=loadArchRProject("/data/processed/aggregate_phase2_ICC_noHarmony")
#####################################################
### RUN HARMONY ACROSS POOLS, Batch Correction
#####################################################

proj2 <- addHarmony(
    ArchRProj = proj2,
    reducedDims = "LSI_ATAC",
    name = "Harmony",
    groupBy = "pool"
)
proj2 <- addUMAP(
    ArchRProj = proj2,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)

#umap_comb_h <- plotEmbedding(proj2, name = "Clusters", embedding = "UMAPHarmony", size = 1.5, labelAsFactors=F, labelMeans=F)
sample_umap_comb_h <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony",force=TRUE)
cluster_umap_comb_h <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony",force=TRUE)
donor_umap_comb_h <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "best_sample", embedding = "UMAPHarmony",force=TRUE)
pool_umap_comb_h <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "pool", embedding = "UMAPHarmony",force=TRUE)

plotPDF(sample_umap_comb_h,name="sample_comb_umap_Harmony.pdf", ArchRProj = proj2, addDOC = FALSE)
plotPDF(cluster_umap_comb_h,name="cluster_comb_umap_Harmony.pdf", ArchRProj = proj2, addDOC = FALSE)
plotPDF(donor_umap_comb_h,name="donor_comb_umap_Harmony.pdf", ArchRProj = proj2, addDOC = FALSE)
plotPDF(pool_umap_comb_h,name="pool_comb_umap_Harmony.pdf", ArchRProj = proj2, addDOC = FALSE)

genes_umap_comb_h <- plotEmbedding(
    ArchRProj = proj2,
    colorBy = "GeneScoreMatrix", 
    name = markerGenes,
    embedding = "UMAPHarmony",
    quantCut = c(0.1, 0.9),
    imputeWeights = NULL
)
plotPDF(plotList=genes_umap_comb_h,name="genes_comb_umap_Harmony.pdf", ArchRProj = proj2, addDOC = FALSE)


### SAVE ARCHR PROJECT
saveArchRProject(ArchRProj = proj2, outputDirectory = "/data/processed/aggregate_phase2_ICC", load = FALSE)






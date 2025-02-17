


call_atac_peaks<-function(proj, outdir, group_name){
    ## take in project, output directory name, and name of cellColData column to use for peak calling	
    pathToMacs2 <- findMacs2()
    print('add group cov')
    proj <- addGroupCoverages(ArchRProj = proj, groupBy = group_name)
    print('add peak set')
    proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = group_name, pathToMacs2 = pathToMacs2)
    print('add peak matrix')
    proj <- addPeakMatrix(proj)
    print('save proj')
    saveArchRProject(ArchRProj = proj, outputDirectory = outdir, load = FALSE)

    ## save each sample separately, this code breaks if the object is too big..
    for( sample in unique(proj@cellColData$Sample)){
        print(sample)
        tmp=proj[which(proj@cellColData$Sample==sample),]
        peaks_mat <- getMatrixFromProject(tmp,useMatrix = "PeakMatrix")
        peaks = as.data.frame(peaks_mat@rowRanges)
	peaks$chr = as.character(peaks$seqnames)
        peaks$peak = paste0(peaks$chr, "_", peaks$start, "_", peaks$end)
        rownames(peaks_mat) = peaks$peak
        
        dir.create(file.path(outdir, sample))
        write(colnames(assay(peaks_mat)), file =  paste0(outdir,'/',sample, '/peak_barcodes.tsv'))
        write(rownames(assay(peaks_mat)), file = paste0(outdir, '/',sample,'/peak_clusters.tsv'))
        writeMM(assay(peaks_mat), file = paste0(outdir, '/',sample,'/peak_coords.tsv'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

    }
    
    ##print('extract peak matrix')
    ##peaks_mat <- getMatrixFromProject(proj,useMatrix = "PeakMatrix")
    ##print(peaks_mat)
    ##print('organize files')

    ###saveArchRProject(proj,outputDirectory = outdir)
    #### peak files
    ##peaks = as.data.frame(peaks_mat@rowRanges)
    ##peaks$chr = as.character(peaks$seqnames)
    ##peaks$peak = paste0(peaks$chr, "_", peaks$start, "_", peaks$end)
    ##rownames(peaks_mat) = peaks$peak
    ##print('write')
    ##write(colnames(assay(peaks_mat)), file =  paste0(outdir, '/peak_barcodes.tsv'))
    ##write(rownames(assay(peaks_mat)), file = paste0(outdir, '/peak_clusters.tsv'))
    ##writeMM(assay(peaks_mat), file = paste0(outdir, '/peak_coords.tsv'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

    ####print('write matrix files')
    ##write(rownames(colData(peaks.mat)), paste0(outdir, '/peak_barcodes.tsv'))
    ##write(colData(peaks.mat)$Clusters, paste0(outdir, '/peak_clusters.tsv'))
    ##write.table(as.data.frame(rowRanges(peaks.mat)), paste0(outdir, '/peak_coords.tsv'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    ##writeMM(assay(peaks.mat), paste0(outdir, '/peaks.mtx'))

    return(proj)
}

#test = getMatrixFromProject(proj, useMatrix = "PeakMatrix")
#peaks = as.data.frame(test@rowRanges)
#peaks$chr = as.character(peaks$seqnames)
#peaks$peak = paste0(peaks$chr, "_", peaks$start, "_", peaks$end)
#rownames(test) = peaks$peak
#write(colnames(assay(test)), file = "MatrixData/ATACPeak_cells_col.txt")
#write(rownames(assay(test)), file = "MatrixData/ATACPeak_peaks_row.txt")
#writeMM(assay(test), file = "MatrixData/ATACPeak_Matrix.txt.gz")

#########################################
## FUNCTIONs TO LOAD PACKAGES AND DATA ##
#########################################

###Load Libraries###
## adjusted to run peaks for all PBMCs together

suppressPackageStartupMessages(library(ArchR))
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(5)
addArchRThreads(1)
library("optparse")
library(data.table)

proj_file='/data/processed/aggregate_phase2_ICC'
outdir='/data/processed/aggregate_phase2_ICC_PBMC_peaks_all'
cellannotfile='/data/rna_files/celltypist/rna.celltypist_allen.csv'
#group_name='AIFI_L1_prediction'
group_name="Tissue"
print(paste('proj file: ',proj_file))
print(paste('outdir: ', outdir))
print(paste('cell annotations: ',cellannotfile))
print(paste('group column: ',group_name))


proj=loadArchRProject(proj_file)
proj

#add cell annotations
df=as.data.frame(fread(cellannotfile))
#df=as.data.frame(fread(cellannotfile))
df$V1=sub(';','#',df$V1)


tmp=merge(proj@cellColData, df, by.x=0, by.y='V1', how='left')
rownames(tmp)=tmp$Row.names
proj@cellColData=tmp
proj@cellColData$Tissue='PBMC'
print(head(proj@cellColData))

#print('save test')
#saveArchRProject(ArchRProj = proj, outputDirectory = outdir, load = FALSE)

#saveArchRProject(proj,outputDirectory = outdir)


proj=call_atac_peaks(proj, outdir, group_name)
#saveArchRProject(proj,outputDirectory = outdir)





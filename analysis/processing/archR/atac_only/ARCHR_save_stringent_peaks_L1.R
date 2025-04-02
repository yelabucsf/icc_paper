# conda activate archr

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(ArchR))
addArchRGenome("hg38")
addArchRThreads(16)

# load archR project
proj=loadArchRProject("/data/processed/aggregate_phase2_ICC_L1_peaks")

# load filtered peak set
filtered_peaks <- readRDS("/data/processed/aggregate_phase2_ICC_L1_peaks/PeakCalls/stringent_peaks/stringent_L1_peaks_spm5e250_rulen12/All_Samples.fwp.filter.non_overlapping.rds")

# save peak matrix per pool
outdir="/data/processed/aggregate_phase2_ICC_L1_peaks_archr_export/stringent_L1_peaks"
dir.create(outdir)
for( sample in unique(proj@cellColData$Sample)){
	print(sample)
	if (!dir.exists(file.path(outdir, sample))) {
	    tmp=proj[which(proj@cellColData$Sample==sample),]
	    tmp <- saveArchRProject(ArchRProj = tmp, outputDirectory = "/data/processed/aggregate_phase2_ICC_L1_peaks_archr_export/stringent_L1_peaks/tmp")
	    print(getOutputDirectory(tmp))

	    tmp <- addPeakSet(tmp, peakSet=filtered_peaks, force=TRUE)
	    tmp <- addPeakMatrix(tmp, threads=48)
	    peaks_mat <- getMatrixFromProject(tmp,useMatrix = "PeakMatrix")
	    peaks = as.data.frame(peaks_mat@rowRanges)
	    print(dim(peaks))
		peaks$chr = as.character(peaks$seqnames)
	    peaks$peak = paste0(peaks$chr, "_", peaks$start, "_", peaks$end)
	    rownames(peaks_mat) = peaks$peak
	        
	    dir.create(file.path(outdir, sample))
	    write(colnames(assay(peaks_mat)), file =  paste0(outdir,'/',sample, '/peak_barcodes.tsv'))
	    write(rownames(assay(peaks_mat)), file = paste0(outdir, '/',sample,'/peak_clusters.tsv'))
	    writeMM(assay(peaks_mat), file = paste0(outdir, '/',sample,'/peak_coords.tsv'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

	    file.remove(paste0("/data/processed/aggregate_phase2_ICC_L1_peaks_archr_export/stringent_L1_peaks/tmp/ArrowFiles/", sample, ".arrow"))
	}
}


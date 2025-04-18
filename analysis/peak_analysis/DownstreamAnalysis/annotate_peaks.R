suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))

option_list <- list(
  make_option(c("--peak_bed"), action="store", type="character", default=NULL),
  make_option(c("--out_bed"), action="store", type="character", default=NULL)
)
opt = parse_args(OptionParser(option_list=option_list))
peak_bed <- opt$peak_bed
out_bed <- opt$out_bed
print(peak_bed)
print(out_bed)

genome <- "hg38"
geneAnno <- createGeneAnnotation(genome)

peaks =  import(peak_bed)
peaksAnno <- ArchR:::.fastAnnoPeaks(peaks, genome, geneAnno)
df <- as.data.frame(peaksAnno)
write.table(df, out_bed, sep="\t", quote=FALSE, row.names=FALSE)




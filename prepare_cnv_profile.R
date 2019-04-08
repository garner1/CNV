library(QDNAseq)
library(DNAcopy)

#bins <- readRDS("bin_annotations_1Mbp.rds")
bins <- getBinAnnotations(binSize=5, genome="hg19")
#saveRDS(bins,file = "bin_annotations_10Mbp.rds")

files <- list.files(path="/media/garner1/hdd2/cutseq/XZ174BICRO158/outdata", pattern="\\.bam$", full.names=T, recursive=FALSE)

for (x in files){
  readCounts <- binReadCounts(bins,bamfiles=x)  
  readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt",smoothBy = 1L)
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  copyNumbersCalled <- callBins(copyNumbersSegmented,method = "cutoff")
  filename<- paste(x, ".pdf", sep="")
  pdf(filename)
  plot(copyNumbersCalled)  
  dev.off()
  filename <- paste(x,".bed", sep="")
  exportBins(copyNumbersSegmented, filename, format = "bed")  
  exportBins(copyNumbersCalled, format="vcf")
  filename <- paste(x,".tsv",sep="")
  exportBins(copyNumbersCalled, filename, type = "segments", format="tsv")
}



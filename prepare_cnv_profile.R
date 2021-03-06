#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(QDNAseq)
library(DNAcopy)

#bins <- readRDS("bin_annotations_1Mbp.rds")
#saveRDS(bins,file = "bin_annotations_10Mbp.rds")

bins <- getBinAnnotations(binSize=as.numeric(args[2]), genome="hg19")
files <- list.files(path=args[1], pattern="\\.bam$", full.names=T, recursive=FALSE)

print(files)
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


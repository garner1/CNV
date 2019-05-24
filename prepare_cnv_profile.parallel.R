#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(QDNAseq)
library(DNAcopy)

bins <- getBinAnnotations(binSize=as.numeric(args[2]), genome="hg19")
readCounts <- binReadCounts(bins,bamfiles=args[1])

readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt",smoothBy = 1L)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
copyNumbersCalled <- callBins(copyNumbersSegmented,method = "cutoff")

filename<- paste(args[1], ".pdf", sep="")
pdf(filename)
plot(copyNumbersCalled)  
dev.off()

#filename <- paste(args[1],".bed", sep="")
#exportBins(copyNumbersSegmented, filename, format = "bed")  
#exportBins(copyNumbersCalled, format="vcf")

filename <- paste(args[1],".tsv",sep="")
exportBins(copyNumbersCalled, filename, type = "segments", format="tsv")


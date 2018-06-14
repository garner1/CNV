library(QDNAseq)
library(DNAcopy)
library(GenVisR)
require(data.table)

bins <- readRDS("bin_annotations_100kbp.rds")

files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/cellLines/HindIII/bamfiles/selection-RC100K", pattern = "\\.bam$", full.names=T, recursive=FALSE)
files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/cellLines/NlaIII/bamfiles/selection-RC100K", pattern = "\\.bam$", full.names=T, recursive=FALSE)

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
  plot(copyNumbersCalled)  
  exportBins(copyNumbersCalled, type = "calls", format="seg")
}
##############
# rm data.seg
# cat *.seg | grep -v SAMPLE > data.seg
##############
data <- read.table("data.seg", header = FALSE)
names(data) <- c("sample","chromosome","start","end","bins","segmean") 

cnSpec(data, genome = "hg19",CNscale = "relative")

plot(copyNumbersCalled)  



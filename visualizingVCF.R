library(QDNAseq)
library(DNAcopy)
library(GenVisR)
require(data.table)

bins <- readRDS("bin_annotations_100kbp.rds")

files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/cellLines/HindIII/bamfiles", pattern = "\\.bam$", full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/cellLines/NlaIII/bamfiles", pattern = "\\.bam$", full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/BreastCancer/HindIII/bamfiles", pattern = "\\.bam$", full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/BreastCancer/NlaIII/bamfiles", pattern = "\\.bam$", full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/TurinSamples/HindIII/bamfiles", pattern = "\\.bam$", full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/TurinSamples/NlaIII/bamfiles", pattern = "\\.bam$", full.names=T, recursive=FALSE)

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
  filename<- paste(x, ".tsv", sep="")
  # exportBins(copyNumbersCalled, type = "segments", format="seg")
  exportBins(copyNumbersCalled, filename, type = "segments", format="tsv")
}
##############

# rm data.seg
# cat *.seg | grep -v SAMPLE > data.seg
##############
data <- read.table("data.tsv", header = FALSE)
names(data) <- c("bins","chromosome","start","end","segmean","sample") 

cnSpec(data, genome = "hg19",CNscale = "relative")

plot(copyNumbersCalled)  



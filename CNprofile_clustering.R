library(QDNAseq)
# bins <- getBinAnnotations(binSize=100, genome="hg19")
# saveRDS(bins,file = "bin_annotations_100kbp.rds")
bins <- readRDS("bin_annotations_100kbp.rds")

# files <- list.files(path="/home/garner1/Work/dataset/CNV/bamfiles/patient_1/selection", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/CNV/bamfiles/patient_2/selection", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/CNV/bamfiles/patient_3", pattern=NULL, full.names=T, recursive=FALSE)

# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/bamfiles_2_cluster_wo33-37/selection_GE100k/selection_GE300k", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/bamfiles_2_cluster_xz33/selection_GE500k", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/bamfiles_2_cluster_xz37/selection_GE500k", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/selection/all_bamfile/selection_GE1M", pattern=NULL, full.names=T, recursive=FALSE)

files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/BreastCancer/NlaIII/bedfiles/", pattern=NULL, full.names=T, recursive=FALSE)

for (x in files){
  readCounts <- binReadCounts(bins,bamfiles=x)  
  readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt",smoothBy = 1L)
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  filename <- paste(x,".bed", sep="")
  exportBins(copyNumbersSegmented, filename, format = "bed")  
}
##########################
#GO TO SELECTION DIRECTORY AND PARSE THE BED FILE TO GET A DATA FILE:
# parallel "tail -n+2 {} |cut -f5 > {.}.signal" ::: *.bed
# paste *.signal | datamash transpose > data.tsv
#TO GET THE LABELS FROM THE SAMPLES
# ls *signal|tr '_.' '\t\t'|cut -f1|datamash transpose |sed 's/\t/","/g'
#########################
temp = read.csv("/home/garner1/Work/dataset/reduced_sequencing/BreastCancer/NlaIII/bedfiles/group_3/data.tsv", header = FALSE, sep="\t")
d <- dist(as.matrix(temp), method = "euclidean")
hc1 <- hclust(d, method = "complete" )
par(cex=0.6, mar=c(5, 8, 4, 1))
plot(hc1,hang=0.1,
     label=c("M10B","M11","M12A","M13A","M13B","M14A","M14B","M15A","M15B","M1A","M1B","M2A","M2B","M3A","M3B","M4A","M5A","M9A","M9B","T10A","T10B","T11","T12","T13","T14","T15","T16","T1","T2","T3","T4A","T4B","T5A","T5B"),
     main="Patients clustering")
#########################
temp = read.csv("/home/garner1/Work/dataset/reduced_sequencing/bamfiles_2_cluster_xz37/selection_GE500k/data.tsv", header = FALSE, sep="\t")
d <- dist(as.matrix(temp), method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1,hang=0.1, 
     label=c("M12A","M13B","M14A","M15A","M1A","M1B","M2A","M2B","M4A","M5A","M5B","M8A","M9B","T12","T1","T2","T4A","T4B","T5A","T5B","T8"),
     main="Patients clustering from XZ37")
#########################
temp = read.csv("/home/garner1/Work/dataset/reduced_sequencing/bamfiles_2_cluster_xz33/selection_GE500k/selection_GE1M/data.tsv", header = FALSE, sep="\t")
d <- dist(as.matrix(temp), method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1,hang=0.1, 
     label=c("M10A","M11","M12A","M13A","M14A","M1A","M1B","M2A","M2B","M5A","M5B","M8B","M9B","T10B","T11","T12","T14","T16","T1","T2","T4A","T4B","T5A","T5B","T8"),
     main="Patients clustering from XZ33")
########################################################
temp = read.csv("/home/garner1/Work/dataset/reduced_sequencing/bamfiles_2_cluster_wo33-37/selection_GE100k/selection_GE300k/selection_GE500k/selection_GE1M/data.tsv", header = FALSE, sep="\t")
d <- dist(as.matrix(temp), method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1,hang=0.1, 
     label=c("M1A-1","M1A-2","M1A-3","M1A-4","M1A-5","M1B-5","M2A-1","M2A-3","M2B-1","M2B-2","M2B-3","M3A-1","T1-2","T1-3","T1-4","T1-5","T2-2","T2-3","T2-4","T2-5","T2-6","T5B-7","T9-4"),
     main="Patients clustering")
########################################################
# temp = read.csv("/home/garner1/Work/dataset/CNV/bamfiles/patient_1/selection/data.tsv", header = FALSE, sep="\t")
# d <- dist(as.matrix(temp), method = "euclidean")
# hc1 <- hclust(d, method = "complete" )
# plot(hc1,hang=0.1,
#      labels=c("M1A-1","M1A-2","M1A-3","M1A-4","M1A-5","M1A","M1B-5","M1B-6","M1B","T1-2","T1-3","T1-5","T1"),
#      main="Patient 1")
# ############################
# temp = read.csv("/home/garner1/Work/dataset/CNV/bamfiles/patient_2/selection/data.tsv", header = FALSE, sep="\t")
# d <- dist(as.matrix(temp), method = "euclidean")
# hc1 <- hclust(d, method = "complete" )
# plot(hc1,hang=0.1,
#      labels=c("M2A-1","M2A-3","M2A","M2B-1","M2B-2","M2B-3","M2B","T2-2","T2-3","T2-4","T2-5","T2-6","T2"),
#      main="Patient 2")
# ###################################
# temp = read.csv("/home/garner1/Work/dataset/CNV/100kbp_30_segmented/PrimaryMetastatBC/XZ33/data.tsv", header = FALSE, sep="\t")
# d <- dist(as.matrix(temp), method = "euclidean")
# hc1 <- hclust(d, method = "complete" )
# plot(hc1,hang=0.1,
#      labels=c("M1B","T5B","T4B","M11","T8","M2A","T1","M8B","M9B","M12A","T11","M2B","M1A","T5A","T4A","M10A","M13A","T12","M5B","M14A","T2","M5A"),
#      main="XZ33 hierarchical clustering ")


library(QDNAseq)
# bins <- getBinAnnotations(binSize=30, genome="hg19")
# saveRDS(bins,file = "bin_annotations_30kbp.rds")
# bins <- readRDS("bin_annotations_100kbp.rds")

# files <- list.files(path="/home/garner1/Work/dataset/CNV/bamfiles/patient_1/selection", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/CNV/bamfiles/patient_2/selection", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/CNV/bamfiles/patient_3", pattern=NULL, full.names=T, recursive=FALSE)

# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/bamfiles_2_cluster_wo33-37/selection_GE100k/selection_GE300k", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/bamfiles_2_cluster_xz33/selection_GE500k", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/bamfiles_2_cluster_xz37/selection_GE500k", pattern=NULL, full.names=T, recursive=FALSE)
# files <- list.files(path="/home/garner1/Work/dataset/reduced_sequencing/selection/all_bamfile/selection_GE1M", pattern=NULL, full.names=T, recursive=FALSE)

# files <- list.files(path="/media/garner1/hdd/cutseq/BreastCancer/HindIII/tsvfiles/selection", pattern=NULL, full.names=T, recursive=FALSE)

#for (x in files){
#  readCounts <- binReadCounts(bins,bamfiles=x)  
#  readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
#  readCountsFiltered <- estimateCorrection(readCountsFiltered)
#  copyNumbers <- correctBins(readCountsFiltered)
#  copyNumbersNormalized <- normalizeBins(copyNumbers)
#  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
#  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt",smoothBy = 1L)
#  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
#  filename <- paste(x,".bed", sep="")
#  exportBins(copyNumbersSegmented, filename, format = "bed")  
#}
##########################
#GO TO SELECTION DIRECTORY AND PARSE THE BED FILE TO GET A DATA FILE:
# parallel "tail -n+2 {} |cut -f5 > {.}.signal" ::: *.bed
# paste *.signal | datamash transpose > data.tsv
#TO GET THE LABELS FROM THE SAMPLES
# ls *signal|tr '_.' '\t\t'|cut -f1|datamash transpose |sed 's/\t/","/g'
#########################
#temp = read.csv("/media/garner1/hdd2/cutseq/clustering/patient_1/data.tsv", header = FALSE, sep="\t")
#d <- dist(as.matrix(temp), method = "euclidean")
#hc1 <- hclust(d, method = "complete" )
#par(cex=0.6, mar=c(5, 8, 4, 1))
#plot(hc1,hang=0.1,
#     label=c("Ma_FS","Ma_LR1","Ma_LR2","Ma_SR2","Ma_SR3","Ma_SR4","Ma_SR5","Ma_SR6","Mb_FS","Mb_LR1","Mb_LR2","Mb_SR2","Mb_SR3","Mb_SR4","#Mb_SR5","T_FS","T_LR1","T_LR2","T_SR2","T_SR3","T_SR4","T_SR5","T_SR6"),
#     main="Patient #1")
#########################
temp = read.csv("/home/garner1/Work/dataset/tsvfiles/data.tsv", header = FALSE, sep="\t")
d <- dist(as.matrix(temp), method = "euclidean")
hc1 <- hclust(d, method = "complete" )
par(cex=0.4, mar=c(5, 8, 4, 1))
plot(hc1,hang=0.1,
     label=c("M10a_FS","M10a_LR1","M10a_LR2","M10b_FS","M10b_LR1","M10b_LR2","M11_FS","M11_LR1","M11_LR2","M12_FS","M12_LR1","M12_LR2","M13a_FS","M13a_LR1","M13a_LR2","M13b_FS","M13b_LR1","M13b_LR2","M14a_FS","M14a_LR1","M14a_LR2","M14b_FS","M14b_LR1","M14b_LR2","M1a_FS","M1a_LR1","M1a_LR2","M1b_FS","M1b_LR1","M1b_LR2","M2a_FS","M2a_LR1","M2a_LR2","M2b_FS","M2b_LR1","M2b_LR2","M3a_FS","M3a_LR1","M3a_LR2","M3b_FS","M3b_LR1","M3b_LR2","M4a_FS","M4a_LR1","M4a_LR2","M4b_FS","M4b_LR1","M4b_LR2","M5_FS","M5_LR1","M5_LR2","M6_FS","M6_LR1","M6_LR2","M7_FS","M7_LR1","M7_LR2","M9_FS","M9_LR1","M9_LR2","T10_FS","T10_LR1","T10_LR2","T11_FS","T11_LR1","T11_LR2","T12_FS","T12_LR1","T12_LR2","T13_FS","T13_LR1","T13_LR2","T14_FS","T14_LR1","T14_LR2","T2a_FS","T2a_LR1","T2a_LR2","T2b_FS","T2b_LR1","T2b_LR2","T3a_FS","T3a_LR1","T3a_LR2","T3b_FS","T3b_LR1","T3b_LR2","T4_FS","T4_LR1","T4_LR2","T5_FS","T5_LR1","T5_LR2","T6_FS","T6_LR1","T6_LR2","T7_FS","T7_LR1","T7_LR2","T8_FS","T8_LR1","T8_LR2","T9_FS","T9_LR1","T9_LR2
"),
     main="Clustering")


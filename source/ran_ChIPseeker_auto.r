
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")

ran.konglab <- function(file){
  library(GenomicFeatures)
  library(ChIPseeker)
  txdb <- makeTxDbFromGFF('genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf', format='gtf', organism='Mus musculus')
  peak <- readPeakFile(file)
  peakAnno <- annotatePeak(peak, TxDb=txdb)
  # add annotations
  peak2 <- read.table(file, header=T)
  peak2 <- peak2[peak2$chr!="chrM",] #This is essential for some cases.
  peak2$feature <- peakAnno@anno$annotation
  peak2$distance <- peakAnno@anno$distanceToTSS
  peak2$gene <- peakAnno@anno$geneId
  file_dir <- dirname(file)
  file_base <- basename(file)
  
  write.table(peak2, paste(file, "_anno.csv", sep = ""), sep='\t', quote=F, row.names=F)
  return()
}

ran.konglab("~/Documents/Anne_Methyl/bam/cov/wk25_tumor_all_results.csv")




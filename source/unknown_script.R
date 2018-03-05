# source("https://bioconductor.org/biocLite.R")
# biocLite("DSS")
# # R
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF('/home/administrator/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf',
                        format='gtf', organism='Mus musculus')
library(ChIPseeker)
peak <- readPeakFile('/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_fr1c1x2_WJ.csv')
peakAnno <- annotatePeak(peak, TxDb=txdb)

# add annotations to comb_def.csv
peak2 <- read.table('/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_fr1c1x2_WJ.csv', header=T)
peak2 <- peak2[peak2$chr!="chrM",]
peak2$feature <- peakAnno@anno$annotation
peak2$distance <- peakAnno@anno$distanceToTSS
peak2$gene <- peakAnno@anno$geneId
write.table(peak2, '/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_fr1c1x2_WJ_anno.csv.csv', sep='\t', quote=F, row.names=F)
#2
peak <- readPeakFile('/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_r1c1x2_WJ.csv')
peakAnno <- annotatePeak(peak, TxDb=txdb)

# add annotations to comb_def.csv
peak2 <- read.table('/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_r1c1x2_WJ.csv', header=T)
peak2 <- peak2[peak2$chr!="chrM",]
peak2$feature <- peakAnno@anno$annotation
peak2$distance <- peakAnno@anno$distanceToTSS
peak2$gene <- peakAnno@anno$geneId
write.table(peak2, '/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_r1c1x2_WJ_anno.csv.csv', sep='\t', quote=F, row.names=F)

##############
##

peak <- readPeakFile('/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_fr1c1x2m1_WJ.csv')
peakAnno <- annotatePeak(peak, TxDb=txdb)

# add annotations to comb_def.csv
peak2 <- read.table('/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_fr1c1x2m1_WJ.csv', header=T)
peak2 <- peak2[peak2$chr!="chrM",]
peak2$feature <- peakAnno@anno$annotation
peak2$distance <- peakAnno@anno$distanceToTSS
peak2$gene <- peakAnno@anno$geneId
write.table(peak2, '/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_fr1c1x2m1_WJ_anno.csv.csv', sep='\t', quote=F, row.names=F)
#2
peak <- readPeakFile('/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_r1c1x2m1_WJ.csv')
peakAnno <- annotatePeak(peak, TxDb=txdb)

# add annotations to comb_def.csv
peak2 <- read.table('/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_r1c1x2m1_WJ.csv', header=T)
peak2 <- peak2[peak2$chr!="chrM",]
peak2$feature <- peakAnno@anno$annotation
peak2$distance <- peakAnno@anno$distanceToTSS
peak2$gene <- peakAnno@anno$geneId
write.table(peak2, '/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/COV_Files/combined_r1c1x2m1_WJ_anno.csv.csv', sep='\t', quote=F, row.names=F)



#or 
write.table(peak2, 'results_anno_T.csv', sep='\t', quote=T, row.names=F)
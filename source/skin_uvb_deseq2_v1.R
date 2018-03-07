# |----------------------------------------------------------------------------------|
# | Project: Skin UVB SKH1 mouse model treated with UA/SFN                           |
# | Script: RNA-seq data analysis and visualization using DESeq2                     |
# | Coordinator: Ran Yin, Renyi Wu                                                   |
# | Author: Davit Sargsyan                                                           |
# | Created: 03/05/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_skin_uvb_deseq2_v1.R")
date()

# Workflow: https://www.bioconductor.org/help/workflows/rnaseqGene/
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

require(data.table)
require(DESeq2)
require(readxl)
require(BiocParallel)
require(ggplot2)

# Load data----
dt1 <- fread("data/renyi_dedup_rnaseq_data/featurescounts_uvb-skin_dedup_renyi_2-9-2018.csv",
             skip = 1)
dt1

# Remove unused columns----
dt1 <- dt1[, c(1, 7:ncol(dt1)), with = FALSE]
dt1

cnames <- colnames(dt1)[-1]
cnames <- gsub(x = cnames,
               pattern = ".dedup.bam",
               replacement = "")
cnames
colnames(dt1)[-1] <- cnames
dt1

# Read-in legend
dt2 <- readxl::read_xlsx("docs/legend.xlsx")
dt2

# Specify treatment groups
mat <- data.frame(sample = cnames,
                  treatment = dt2$Description,
                  trt = rep(rep(c("ng", "pc", "ua", "sfn"), 
                            each = 2),
                            4),
                  week = dt2$Time,
                  time = factor(rep(1:4, each = 8)),
                  repl = factor(rep(1:16, each = 2)))
mat

dtm <- as.matrix(dt1[, -1, with = FALSE])
rownames(dtm) <- dt1$Geneid

# Part I: build a DESeq2 data----
dds <- DESeqDataSetFromMatrix(countData = dtm, 
                              colData = mat,
                              ~ trt)
dds

# Remove genes with low counts (<100 in all samples combined)
dds <- dds[rowSums(counts(dds)) >= 100, ]
nrow(counts(dds))
# 15,385 genes left, down from 24,421
# If all samples contain zeros, geometric means cannot be
# estimated. Change default 'type = "ratio"' to 'type = "poscounts"'.
# Type '?DESeq2::estimateSizeFactors' for more details.
dds <- estimateSizeFactors(dds,
                           type = "poscounts")
dds

# Clustering and PCA----
rld <- rlog(dds, 
            blind = FALSE)
head(assay(rld))

# Clustering----
sampleDists <- dist(t(assay(rld)))
sampleDists

tiff(filename = "tmp/skin_uvb_samples_cluster.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
heatmap(as.matrix(sampleDists),
        symm = TRUE)
graphics.off()

# PCA----
m1 <- prcomp(t(assay(rld)),
             center = TRUE,
             scale. = TRUE)
summary(m1)
plot(m1)

# Biplot while keep only the most important variables (Javier)----
# Select PC-s to pliot (PC1 & PC2)
choices <- 1:2

# Scores, i.e. points (df.u)
dt.scr <- data.table(m1$x[, choices])
# Add grouping variable
# dt.scr$grp <- colnames(assay(rld))
dt.scr$grp <- mat$treatment
dt.scr$sample <- mat$sample
dt.scr

# Loadings, i.e. arrows (df.v)
dt.rot <- as.data.frame(m1$rotation[, choices])
dt.rot$feat <- rownames(dt.rot)
dt.rot <- data.table(dt.rot)
dt.rot

# Axis labels
u.axis.labs <- paste(colnames(dt.rot)[1:2], 
                     sprintf('(%0.1f%% explained var.)', 
                             100*m1$sdev[choices]^2/sum(m1$sdev^2)))
u.axis.labs

# Based on Figure p0, keep only a few variables with high loadings in PC1 and PC2----
# var.keep.ndx <- which(dt.rot$feat %in% c("V1",
#                                          "V4",
#                                          "V6",
#                                          "V7",
#                                          "V8"))
gl1 <-dt.rot$feat[order(abs(dt.rot$PC1),
                        decreasing = TRUE)][1:5]
gl2 <-dt.rot$feat[order(abs(dt.rot$PC2),
                        decreasing = TRUE)][1:5]
var.keep.ndx <- which(dt.rot$feat %in% unique(c(gl1, gl2)))
# Or select all
# var.keep.ndx <- 1:ncol(assay(rld))

p1 <- ggplot(data = dt.rot[var.keep.ndx,],
             aes(x = PC1,
                 y = PC2)) +
  # coord_equal() +
  geom_point(data = dt.scr,
             aes(fill = grp),
             shape = 21,
             size = 3,
             alpha = 0.5) +
  # geom_segment(aes(x = 0,
  #                  y = 0,
  #                  xend = 10000*PC1,
  #                  yend = 10000*PC2),
  #              arrow = arrow(length = unit(1/2, 'picas')),
  #              color = "black",
  #              size = 0.5) +
  # geom_text(aes(x = 11000*PC1,
  #               y = 11000*PC2,
  #               label = dt.rot$feat[var.keep.ndx]),
  #           size = 2,
  #           hjust = 0.2) +
geom_text(data = dt.scr,
          aes(x = PC1 + 20,
              y = PC2,
              label = dt.scr$sample),
          size = 2,
          hjust = 0.5) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  scale_fill_manual(name = "Treatment",
                    values = c("white", "red", "blue", "green")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20))
p1
tiff(filename = "tmp/pca_biplot.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# # Alternative PCA (form the package)----
# tiff(filename = "tmp/skin_uvb_samples_pca.tiff",
#      height = 8,
#      width = 8,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# plotPCA(rld, 
#         intgroup = c("trt"))
# graphics.off()

# Part II: build a DESeq2 model with treatment alone----
# Set cores for parallel processing of DESeq----
snowparam <- SnowParam(workers = snowWorkers(), 
                       type = "SOCK")
register(snowparam, 
         default = TRUE)

# Run DESeq----
dds <- DESeq(dds,
             fitType = "local",
             parallel = TRUE)
results(dds)
resultsNames(dds)

# Contrasts----
# a. Positive (UVB) vs. negative controls----
pc_nc <- results(dds,
                 name = "trt_pc_vs_ng")

t1 <- data.table(do.call("cbind", 
                         pc_nc@listData))
write.csv(t1,
          file = "tmp/pc_vs_nc.csv",
          row.names = FALSE)

t1$clr <- "black"
t1$clr[t1$pvalue < 0.05 & abs(t1$log2FoldChange) >= 2] <- "red"
t1$pch <- 46
t1$pch[t1$pvalue < 0.05 & abs(t1$log2FoldChange) >= 2] <- 3
t1
summary(t1)
# 'baseMean' is the average of the normalized count values, 
# divided by the size factors, taken over all samples in the DESeqDataSet
tiff(filename = "tmp/skin_uvb_pc.vs.nc.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(t1$log2FoldChange ~ log(t1$baseMean + 1),
     col = t1$clr,
     pch = t1$pch,
     xlab = "log(baseMean+1)",
     ylab = "log2 Fold Change",
     main = "Positive vs. Negative Control")
abline(h = c(-2, 2),
       lty = 2)
graphics.off()

# b. Ursolic Acid (UA) vs. Positive (UVB)----
ua_pc <- results(dds,
                 contrast = c("trt", "ua", "pc"))
t2 <- data.table(do.call("cbind", 
                         ua_pc@listData))
write.csv(t2,
          file = "tmp/ua_vs_pc.csv",
          row.names = FALSE)

t2$clr <- "black"
t2$clr[t2$pvalue < 0.05 & abs(t2$log2FoldChange) >= 1] <- "red"
t2$pch <- 46
t2$pch[t2$pvalue < 0.05 & abs(t2$log2FoldChange) >= 1] <- 3
t2
summary(t2)

tiff(filename = "tmp/skin_uvb_ua.vs.pc.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(t2$log2FoldChange ~ log(t2$baseMean + 1),
     col = t2$clr,
     pch = t2$pch,
     xlab = "log(baseMean+1)",
     ylab = "log2 Fold Change",
     main = "Ursolic Acid vs. Positive Control")
abline(h = c(-1, 1),
       lty = 2)
graphics.off()

# c. SFN vs. Positive (UVB)----
sfn_pc <- results(dds,
                 contrast = c("trt", "sfn", "pc"))
t3 <- data.table(do.call("cbind", 
                         sfn_pc@listData))
write.csv(t3,
          file = "tmp/sfn_vs_pc.csv",
          row.names = FALSE)

t3$clr <- "black"
t3$clr[t3$pvalue < 0.05 & abs(t3$log2FoldChange) >= 1] <- "red"
t3$pch <- 46
t3$pch[t3$pvalue < 0.05 & abs(t3$log2FoldChange) >= 1] <- 3
t3
summary(t3)

tiff(filename = "tmp/skin_uvb_sfn.vs.pc.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(t3$log2FoldChange ~ log(t3$baseMean + 1),
     col = t3$clr,
     pch = t3$pch,
     xlab = "log(baseMean+1)",
     ylab = "log2 Fold Change",
     main = "SFN vs. Positive Control")
abline(h = c(-1, 1),
       lty = 2)
graphics.off()

# Part III: build a DESeq2 model with treatment*time interaction----
dds <- DESeqDataSetFromMatrix(countData = as.matrix(dt1[, -1, with = FALSE]),
                              colData = mat,
                              ~ trt + time + trt:time)
dds

# Remove genes with low counts (<100 in all samples combined)
dds <- dds[rowSums(counts(dds)) >= 100, ]
nrow(counts(dds))
# 15,385 genes left, down from 24,421

# If all samples contain zeros, geometric means cannot be
# estimated. Change default 'type = "ratio"' to 'type = "poscounts"'.
# Type '?DESeq2::estimateSizeFactors' for more details.
dds <- estimateSizeFactors(dds,
                           type = "poscounts")
dds

# Set cores for parallel processing of DESeq----
snowparam <- SnowParam(workers = snowWorkers(), 
                       type = "SOCK")
register(snowparam, 
         default = TRUE)

# Run DESeq----
dds <- DESeq(dds,
             fitType = "local",
             parallel = TRUE)
results(dds)
resultsNames(dds)

# Contrasts----
# a. Positive (UVB) vs. negative controls at week 2----
pc_nc_w2 <- results(dds,
                    name = "trt_pc_vs_ng")
t1w2 <- data.table(do.call("cbind", 
                           pc_nc@listData))
write.csv(t1w2,
          file = "tmp/pc_vs_nc_w2.csv",
          row.names = FALSE)

t1w2$clr <- "black"
t1w2$clr[t1w2$pvalue < 0.05 & abs(t1w2$log2FoldChange) >= 2] <- "red"
t1w2$pch <- 46
t1w2$pch[t1w2$pvalue < 0.05 & abs(t1w2$log2FoldChange) >= 2] <- 3
t1w2
summary(t1w2)

tiff(filename = "tmp/skin_uvb_pc.vs.nc.w2.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(t1w2$log2FoldChange ~ log(t1w2$baseMean + 1),
     col = t1w2$clr,
     pch = t1w2$pch,
     xlab = "log(baseMean+1)",
     ylab = "log2 Fold Change",
     main = "Positive vs. Negative Control")
abline(h = c(-2, 2),
       lty = 2)
graphics.off()

# b. Positive (UVB) vs. negative controls at week 15----
pc_nc_w15 <- results(dds,
                     contrast = list(c("trt_pc_vs_ng",
                                       "trtpc.time2")))
t1w15 <- data.table(do.call("cbind", 
                            pc_nc_w15@listData))
write.csv(t1w15,
          file = "tmp/pc_vs_nc_w15.csv",
          row.names = FALSE)

t1w15$clr <- "black"
t1w15$clr[t1w15$pvalue < 0.05 & abs(t1w15$log2FoldChange) >= 2] <- "red"
t1w15$pch <- 46
t1w15$pch[t1w15$pvalue < 0.05 & abs(t1w15$log2FoldChange) >= 2] <- 3
t1w15
summary(t1w15)

tiff(filename = "tmp/skin_uvb_pc.vs.nc.w15.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(t1w15$log2FoldChange ~ log(t1w15$baseMean + 1),
     col = t1w15$clr,
     pch = t1w15$pch,
     xlab = "log(baseMean+1)",
     ylab = "log2 Fold Change",
     main = "Positive vs. Negative Control")
abline(h = c(-2, 2),
       lty = 2)
graphics.off()

# sink()
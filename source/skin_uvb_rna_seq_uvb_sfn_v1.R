# |----------------------------------------------------------------------------------|
# | Project: Skin UVB SKH1 mouse model treated with UA/SFN                           |
# | Script: RNA-seq data analysis and visualization using DESeq2                     |
# | Coordinator: Yuqing (SFN) Yang, Ran Yin, Renyi Wu                                |
# | Author: Davit Sargsyan                                                           |
# | Created: 03/05/2018                                                              |
# | Modified: 03/10/2018, DS: added trt*time model and pairwise comparisons          |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_skin_uvb_rna_seq_uvb_sfn_v1.txt")
date()

# Workflow: https://www.bioconductor.org/help/workflows/rnaseqGene/
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

require(data.table)
require(DESeq2)
require(readxl)
require(BiocParallel)
require(ggplot2)
require(knitr)

# NOTE on DESeq2 Output: 'baseMean' is the average of the normalized count values, 
# divided by the size factors, taken over all samples in the DESeqDataSet

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
dt2 <- readxl::read_xlsx("docs/legend_v2.xlsx")
dt2

# Specify treatment groups
mat <- data.frame(sample = cnames,
                  treatment = dt2$Description,
                  trt = rep(rep(c("nc", "sfn", "ua", "pc"), 
                            each = 2),
                            4),
                  week = dt2$Time,
                  time = factor(rep(c(1, 2, 4, 3), each = 8)),
                  repl = factor(rep(1:16, each = 2)))
mat

dtm <- as.matrix(dt1[, -1, with = FALSE])
rownames(dtm) <- dt1$Geneid

# Keep only controls and UA treatment----
mat <- droplevels(mat[mat$trt != "sfn", ])
mat
dtm <- dtm[, as.character(mat$sample)]
dtm

# Part I: build a DESeq2 data----
dds <- DESeqDataSetFromMatrix(countData = dtm, 
                              colData = mat,
                              ~ trt + time)
dds

# Remove genes with low counts (<100 in all samples combined)
dds <- dds[rowSums(counts(dds)) >= 100, ]
nrow(counts(dds))
# 15,055 genes left, down from 24,421
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

# PCA, no tumor----
# Keep only controls and UA treatment----
mat <- droplevels(mat[mat$time != 4, ])
mat
dtm <- dtm[, as.character(mat$sample)]
dtm

# Part I: build a DESeq2 data----
dds <- DESeqDataSetFromMatrix(countData = dtm, 
                              colData = mat,
                              ~ trt + time + trt:time)
dds

# Remove genes with low counts (<100 in all samples combined)
dds <- dds[rowSums(counts(dds)) >= 100, ]
nrow(counts(dds))
# 13,372 genes left, down from ??

rld <- rlog(dds, 
            blind = FALSE)
head(assay(rld))

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

p2 <- ggplot(data = dt.rot[var.keep.ndx,],
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
p2
tiff(filename = "tmp/pca_biplot_no_tumor.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
# gridExtra::grid.arrange(p1, p2)
plot(p2)
graphics.off()

# Part II: build a DESeq2 model with treatment*time interaction----
dds <- DESeqDataSetFromMatrix(countData = dtm,
                              colData = mat,
                              ~ trt + time + trt:time)
dds

# Remove genes with low counts (<100 in all samples combined)
dds <- dds[rowSums(counts(dds)) >= 100, ]
nrow(counts(dds))
# 13,372 genes left, down from 24,421

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
colData(dds)

# Contrasts----
# i. Week 2----
# a. Positive (UVB) vs. negative controls at week 2----
pc_nc_w2 <- results(dds,
                    name = "trt_pc_vs_nc")
t1w2 <- data.table(gene = rownames(pc_nc_w2),
                   do.call("cbind", 
                           pc_nc_w2@listData))
t1w2
write.csv(t1w2,
          file = "tmp/pc_vs_nc_w2.csv",
          row.names = FALSE)

t1w2$clr <- "black"
t1w2$clr[t1w2$pvalue < 0.05] <- "purple"
t1w2$clr[t1w2$pvalue < 0.05 & abs(t1w2$log2FoldChange) >= 2] <- "red"
t1w2$pch <- 46
t1w2$pch[t1w2$pvalue < 0.05] <- 3
t1w2$pch[t1w2$pvalue < 0.05 & abs(t1w2$log2FoldChange) >= 2] <- 4
t1w2$pch <- factor(t1w2$pch,
                 levels = c(46, 3, 4))
t1w2
kable(table(t1w2$clr,
            t1w2$pch))
  # |       |   46|    3|   4|
  # |:------|----:|----:|---:|
  # |black  | 9078|    0|   0|
  # |purple |    0| 3471|   0|
  # |red    |    0|    0| 823|

p1 <- ggplot(t1w2,
             aes(x = log(baseMean + 1),
                 y = log2FoldChange,
                 colour = clr,
                 shape = pch)) +
  scale_shape_manual(name = "Significance",
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2"),
                     values = c(46, 3, 4)) +
  scale_color_manual(name = "Significance",
                     #guide = "legend",
                     values = c("black",
                                "purple",
                                "red"),
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2")) +
  geom_hline(yintercept = c(-2, 2),
             lty = 2) +
  ggtitle("UVB vs. Negative Control, Week 2")  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_pc.vs.nc.w2.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# b. UA vs. Positive (UVB) week 2----
ua_pc_w2 <- results(dds,
                    contrast = list(c("trt_ua_vs_nc",
                                      "trt_pc_vs_nc")))
t2w2 <- data.table(gene = rownames(ua_pc_w2),
                   do.call("cbind", 
                           ua_pc_w2@listData))
t2w2
write.csv(t2w2,
          file = "tmp/ua_vs_pc_w2.csv",
          row.names = FALSE)

t2w2$clr <- "black"
t2w2$clr[t2w2$pvalue < 0.05] <- "purple"
t2w2$clr[t2w2$pvalue < 0.05 & abs(t2w2$log2FoldChange) >= 2] <- "red"
t2w2$pch <- 46
t2w2$pch[t2w2$pvalue < 0.05] <- 3
t2w2$pch[t2w2$pvalue < 0.05 & abs(t2w2$log2FoldChange) >= 2] <- 4
t2w2$pch <- factor(t2w2$pch,
                   levels = c(46, 3, 4))
t2w2
kable(table(t2w2$clr,
            t2w2$pch))
  # |       |    46|    3|    4|
  # |:------|-----:|----:|----:|
  # |black  | 11440|    0|    0|
  # |purple |     0| 2263|    0|
  # |red    |     0|    0| 1682|

p1 <- ggplot(t2w2,
             aes(x = log(baseMean + 1),
                 y = log2FoldChange,
                 colour = clr,
                 shape = pch)) +
  scale_shape_manual(name = "Significance",
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2"),
                     values = c(46, 3, 4)) +
  scale_color_manual(name = "Significance",
                     values = c("black",
                                "purple",
                                "red"),
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2")) +
  geom_hline(yintercept = c(-2, 2),
             lty = 2) +
  ggtitle("UA vs. UVB, Week 2")  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_ua.vs.pc.w2.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# ii. Week 15----
resultsNames(dds)

# a. Positive (UVB) vs. negative controls at week 2----
pc_nc_w15 <- results(dds,
                    name = "trtpc.time2")
t1w15 <- data.table(gene = rownames(pc_nc_w15),
                   do.call("cbind", 
                           pc_nc_w15@listData))
t1w15
write.csv(t1w15,
          file = "tmp/pc_vs_nc_w15.csv",
          row.names = FALSE)

t1w15$clr <- "black"
t1w15$clr[t1w15$pvalue < 0.05] <- "purple"
t1w15$clr[t1w15$pvalue < 0.05 & abs(t1w15$log2FoldChange) >= 2] <- "red"
t1w15$pch <- 46
t1w15$pch[t1w15$pvalue < 0.05] <- 3
t1w15$pch[t1w15$pvalue < 0.05 & abs(t1w15$log2FoldChange) >= 2] <- 4
t1w15$pch <- factor(t1w15$pch,
                   levels = c(46, 3, 4))
t1w15
kable(table(t1w15$clr,
            t1w15$pch))
  # |       |    46|   3|   4|
  # |:------|-----:|---:|---:|
  # |black  | 14348|   0|   0|
  # |purple |     0| 639|   0|
  # |red    |     0|   0| 398|

p1 <- ggplot(t1w15,
             aes(x = log(baseMean + 1),
                 y = log2FoldChange,
                 colour = clr,
                 shape = pch)) +
  scale_shape_manual(name = "Significance",
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2"),
                     values = c(46, 3, 4)) +
  scale_color_manual(name = "Significance",
                     #guide = "legend",
                     values = c("black",
                                "purple",
                                "red"),
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2")) +
  geom_hline(yintercept = c(-2, 2),
             lty = 2) +
  ggtitle("UVB vs. Negative Control, Week 15")  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_pc.vs.nc.w15.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# b. UA vs. Positive (UVB) week 15----
ua_pc_w15 <- results(dds,
                    contrast = list(c("trtua.time2",
                                      "trtpc.time2")))
t2w15 <- data.table(gene = rownames(ua_pc_w15),
                   do.call("cbind", 
                           ua_pc_w15@listData))
t2w15
write.csv(t2w15,
          file = "tmp/ua_vs_pc_w15.csv",
          row.names = FALSE)

t2w15$clr <- "black"
t2w15$clr[t2w15$pvalue < 0.05] <- "purple"
t2w15$clr[t2w15$pvalue < 0.05 & abs(t2w15$log2FoldChange) >= 2] <- "red"
t2w15$pch <- 46
t2w15$pch[t2w15$pvalue < 0.05] <- 3
t2w15$pch[t2w15$pvalue < 0.05 & abs(t2w15$log2FoldChange) >= 2] <- 4
t2w15$pch <- factor(t2w15$pch,
                   levels = c(46, 3, 4))
t2w15
kable(table(t2w15$clr,
            t2w15$pch))
# |       |    46|   3|   4|
# |:------|-----:|---:|---:|
# |black  | 14073|   0|   0|
# |purple |     0| 534|   0|
# |red    |     0|   0| 778|

p1 <- ggplot(t2w15,
             aes(x = log(baseMean + 1),
                 y = log2FoldChange,
                 colour = clr,
                 shape = pch)) +
  scale_shape_manual(name = "Significance",
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2"),
                     values = c(46, 3, 4)) +
  scale_color_manual(name = "Significance",
                     #guide = "legend",
                     values = c("black",
                                "purple",
                                "red"),
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2")) +
  geom_hline(yintercept = c(-2, 2),
             lty = 2) +
  ggtitle("UA vs. UVB, Week 15")  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_ua.vs.pc.w15.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# iii. Week 15----
resultsNames(dds)

# a. Positive (UVB) vs. negative controls at week 2----
pc_nc_w25 <- results(dds,
                     name = "trtpc.time3")
t1w25 <- data.table(gene = rownames(pc_nc_w25),
                    do.call("cbind", 
                            pc_nc_w25@listData))
t1w25
write.csv(t1w25,
          file = "tmp/pc_vs_nc_w25.csv",
          row.names = FALSE)

t1w25$clr <- "black"
t1w25$clr[t1w25$pvalue < 0.05] <- "purple"
t1w25$clr[t1w25$pvalue < 0.05 & abs(t1w25$log2FoldChange) >= 2] <- "red"
t1w25$pch <- 46
t1w25$pch[t1w25$pvalue < 0.05] <- 3
t1w25$pch[t1w25$pvalue < 0.05 & abs(t1w25$log2FoldChange) >= 2] <- 4
t1w25$pch <- factor(t1w25$pch,
                    levels = c(46, 3, 4))
t1w25
kable(table(t1w25$clr,
            t1w25$pch))
# |       |    46|   3|   4|
# |:------|-----:|---:|---:|
# |black  | 14348|   0|   0|
# |purple |     0| 639|   0|
# |red    |     0|   0| 398|

p1 <- ggplot(t1w25,
             aes(x = log(baseMean + 1),
                 y = log2FoldChange,
                 colour = clr,
                 shape = pch)) +
  scale_shape_manual(name = "Significance",
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2"),
                     values = c(46, 3, 4)) +
  scale_color_manual(name = "Significance",
                     #guide = "legend",
                     values = c("black",
                                "purple",
                                "red"),
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2")) +
  geom_hline(yintercept = c(-2, 2),
             lty = 2) +
  ggtitle("UVB vs. Negative Control, Week 25")  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_pc.vs.nc.w25.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# b. UA vs. Positive (UVB) week 25----
ua_pc_w25 <- results(dds,
                     contrast = list(c("trtua.time2",
                                       "trtpc.time2")))
t2w25 <- data.table(gene = rownames(ua_pc_w25),
                    do.call("cbind", 
                            ua_pc_w25@listData))
t2w25
write.csv(t2w25,
          file = "tmp/ua_vs_pc_w25.csv",
          row.names = FALSE)

t2w25$clr <- "black"
t2w25$clr[t2w25$pvalue < 0.05] <- "purple"
t2w25$clr[t2w25$pvalue < 0.05 & abs(t2w25$log2FoldChange) >= 2] <- "red"
t2w25$pch <- 46
t2w25$pch[t2w25$pvalue < 0.05] <- 3
t2w25$pch[t2w25$pvalue < 0.05 & abs(t2w25$log2FoldChange) >= 2] <- 4
t2w25$pch <- factor(t2w25$pch,
                    levels = c(46, 3, 4))
t2w25
kable(table(t2w25$clr,
            t2w25$pch))
# |       |    46|   3|   4|
# |:------|-----:|---:|---:|
# |black  | 14073|   0|   0|
# |purple |     0| 534|   0|
# |red    |     0|   0| 778|

p1 <- ggplot(t2w25,
             aes(x = log(baseMean + 1),
                 y = log2FoldChange,
                 colour = clr,
                 shape = pch)) +
  scale_shape_manual(name = "Significance",
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2"),
                     values = c(46, 3, 4)) +
  scale_color_manual(name = "Significance",
                     #guide = "legend",
                     values = c("black",
                                "purple",
                                "red"),
                     labels = c("No significance",
                                "p-Value < 0.05",
                                "p-Value < 0.05 & abs(log2) >= 2")) +
  geom_hline(yintercept = c(-2, 2),
             lty = 2) +
  ggtitle("UA vs. UVB, Week 25")  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_ua.vs.pc.w25.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# sessionInfo()
# sink()
# |----------------------------------------------------------------------------------|
# | Project: Skin UVB SKH1 mouse model treated with UA/SFN                           |
# | Script: RNA-seq data analysis and visualization using DESeq2                     |
# | Coordinator: Ran Yin, Renyi Wu                                                   |
# | Author: Davit Sargsyan                                                           |
# | Created: 03/05/2018                                                              |
# | Modified: 03/10/2018, DS: added trt*time model and pairwise comparisons          |
# |----------------------------------------------------------------------------------|
sink(file = "tmp/log_skin_uvb_deseq2_v1.txt")
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
dt2 <- readxl::read_xlsx("docs/legend.xlsx")
dt2

# Specify treatment groups
mat <- data.frame(sample = cnames,
                  treatment = dt2$Description,
                  trt = rep(rep(c("nc", "pc", "ua", "sfn"), 
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
                 name = "trt_pc_vs_nc")

t1 <- data.table(gene = rownames(pc_nc),
                 do.call("cbind", 
                         pc_nc@listData))
t1
write.csv(t1,
          file = "tmp/pc_vs_nc.csv",
          row.names = FALSE)

t1$clr <- "black"
t1$clr[t1$pvalue < 0.05] <- "purple"
t1$clr[t1$pvalue < 0.05 & abs(t1$log2FoldChange) >= 2] <- "red"
t1$pch <- 46
t1$pch[t1$pvalue < 0.05] <- 3
t1$pch[t1$pvalue < 0.05 & abs(t1$log2FoldChange) >= 2] <- 4
t1$pch <- factor(t1$pch,
                 levels = c(46, 3, 4))
t1
kable(table(t1$clr,
            t1$pch))
  # |       |    46|    3|   4|
  # |:------|-----:|----:|---:|
  # |black  | 12114|    0|   0|
  # |purple |     0| 2898|   0|
  # |red    |     0|    0| 373|

p1 <- ggplot(t1,
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
  ggtitle("UVB vs. Negative Control, All Timepoints") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_pc.vs.nc.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# b. Ursolic Acid (UA) vs. Positive (UVB)----
ua_pc <- results(dds,
                 contrast = c("trt", "ua", "pc"))
t2 <- data.table(gene = rownames(ua_pc),
                 do.call("cbind", 
                         ua_pc@listData))
t2
write.csv(t2,
          file = "tmp/ua_vs_pc.csv",
          row.names = FALSE)

t2$clr <- "black"
t2$clr[t2$pvalue < 0.05] <- "purple"
t2$clr[t2$pvalue < 0.05 & abs(t2$log2FoldChange) >= 2] <- "red"
t2$pch <- 46
t2$pch[t2$pvalue < 0.05] <- 3
t2$pch[t2$pvalue < 0.05 & abs(t2$log2FoldChange) >= 2] <- 4
t2$pch <- factor(t2$pch,
                 levels = c(46, 3, 4))
t2
kable(table(t2$clr,
            t2$pch))
  # |       |    46|   3|  4|
  # |:------|-----:|---:|--:|
  # |black  | 15160|   0|  0|
  # |purple |     0| 210|  0|
  # |red    |     0|   0| 15|

p1 <- ggplot(t2,
             aes(x = log(baseMean + 1),
                 y = log2FoldChange,
                 colour = clr,
                 shape = pch)) +
  scale_shape_manual(name = "Significance",
                     labels = c("No significance",
                                "p-value < 0.05",
                                "p-value < 0.05 & abs(log2) >= 2"),
                     values = c(46, 3, 4)) +
  scale_color_manual(name = "Significance",
                     #guide = "legend",
                     values = c("black",
                                "purple",
                                "red"),
                     labels = c("No significance",
                                "p-value < 0.05",
                                "p-value < 0.05 & abs(log2) >= 2")) +
  geom_hline(yintercept = c(-2, 2),
             lty = 2) +
  ggtitle("UA vs. UVB") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_ua.vs.pc.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# c. SFN vs. Positive (UVB)----
sfn_pc <- results(dds,
                  contrast = c("trt", "sfn", "pc"))
t3 <- data.table(gene = rownames(sfn_pc),
                 do.call("cbind", 
                         sfn_pc@listData))
t3
write.csv(t3,
          file = "tmp/sfn_vs_pc.csv",
          row.names = FALSE)

t3$clr <- "black"
t3$clr[t3$pvalue < 0.05] <- "purple"
t3$clr[t3$pvalue < 0.05 & abs(t3$log2FoldChange) >= 2] <- "red"
t3$pch <- 46
t3$pch[t3$pvalue < 0.05] <- 3
t3$pch[t3$pvalue < 0.05 & abs(t3$log2FoldChange) >= 2] <- 4
t3$pch <- factor(t3$pch,
                 levels = c(46, 3, 4))
t3
kable(table(t3$clr,
            t3$pch))
  # |       |    46|   3|  4|
  # |:------|-----:|---:|--:|
  # |black  | 15184|   0|  0|
  # |purple |     0| 180|  0|
  # |red    |     0|   0| 21|

p1 <- ggplot(t3,
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
  ggtitle("SFN vs. UVB, All Timepoints") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_sfn.vs.pc.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Part III: build a DESeq2 model with treatment*time interaction----
dds <- DESeqDataSetFromMatrix(countData = dtm,
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
  # |       |    46|    3|   4|
  # |:------|-----:|----:|---:|
  # |black  | 11294|    0|   0|
  # |purple |     0| 3303|   0|
  # |red    |     0|    0| 788|

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

# c. SFN vs. Positive (UVB) week 2----
sfn_pc_w2 <- results(dds,
                    contrast = list(c("trt_sfn_vs_nc",
                                      "trt_pc_vs_nc")))
t3w2 <- data.table(gene = rownames(sfn_pc_w2),
                   do.call("cbind", 
                           sfn_pc_w2@listData))
t3w2
write.csv(t3w2,
          file = "tmp/sfn_vs_pc_w2.csv",
          row.names = FALSE)

t3w2$clr <- "black"
t3w2$clr[t3w2$pvalue < 0.05] <- "purple"
t3w2$clr[t3w2$pvalue < 0.05 & abs(t3w2$log2FoldChange) >= 2] <- "red"
t3w2$pch <- 46
t3w2$pch[t3w2$pvalue < 0.05] <- 3
t3w2$pch[t3w2$pvalue < 0.05 & abs(t3w2$log2FoldChange) >= 2] <- 4
t3w2$pch <- factor(t3w2$pch,
                   levels = c(46, 3, 4))
t3w2
kable(table(t3w2$clr,
            t3w2$pch))
  # |       |    46|    3|    4|
  # |:------|-----:|----:|----:|
  # |black  | 10023|    0|    0|
  # |purple |     0| 3013|    0|
  # |red    |     0|    0| 2349|

p1 <- ggplot(t3w2,
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
  ggtitle("SFN vs. UVB, Week 2")  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_sfn.vs.pc.w2.tiff",
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

# c. SFN vs. Positive (UVB) week 15----
sfn_pc_w15 <- results(dds,
                     contrast = list(c("trtsfn.time2",
                                       "trtpc.time2")))
t3w15 <- data.table(gene = rownames(sfn_pc_w15),
                   do.call("cbind", 
                           sfn_pc_w15@listData))
t3w15
write.csv(t3w15,
          file = "tmp/sfn_vs_pc_w15.csv",
          row.names = FALSE)

t3w15$clr <- "black"
t3w15$clr[t3w15$pvalue < 0.05] <- "purple"
t3w15$clr[t3w15$pvalue < 0.05 & abs(t3w15$log2FoldChange) >= 2] <- "red"
t3w15$pch <- 46
t3w15$pch[t3w15$pvalue < 0.05] <- 3
t3w15$pch[t3w15$pvalue < 0.05 & abs(t3w15$log2FoldChange) >= 2] <- 4
t3w15$pch <- factor(t3w15$pch,
                   levels = c(46, 3, 4))
t3w15
kable(table(t3w15$clr,
            t3w15$pch))
# |       |    46|   3|   4|
# |:------|-----:|---:|---:|
# |black  | 13934|   0|   0|
# |purple |     0| 630|   0|
# |red    |     0|   0| 821|

p1 <- ggplot(t3w15,
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
  ggtitle("SFN vs. UVB, Week 2")  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") + 
  geom_point()
p1

tiff(filename = "tmp/skin_uvb_sfn.vs.pc.w15.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

sink()
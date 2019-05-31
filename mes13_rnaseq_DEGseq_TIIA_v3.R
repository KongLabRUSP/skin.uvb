# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: RNA-seq data analysis and visualization using edgeR, TIIA only (Wenji)   |
# | Author: Davit Sargsyan                                                           |
# | Created: 01/29/2018                                                              |
# | Modified: 05/21/2018(DS): added dendrogram (phylogenic style)                    |
# |----------------------------------------------------------------------------------|
sink(file = "tmp/log_mes13_rnaseq_DEGseq_TIIA_v3.R")
# Source: 
# https://bioconductor.org/packages/release/bioc/html/DEGseq.html

# if (!requireNamespace("BiocManager",
#                       quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("DEGseq")

# Header----
require(data.table)
require(ggplot2)
require(DESeq2)
require(DEGseq)
require(knitr)
require(ggdendro)

# MES13 data----
# Treatment legend----
trt.names <- c("LG",
               "HG",
               "MIC_1.5uM",
               "TIIA_5uM",
               "FX_1uM",
               "Gen_10uM",
               "Ber_6uM")

# Load data----
# Question to Renyi: how was the data processed and annotated?
# dt1 <- fread("data/Renyi_RNAseq_12292017/mes13_fpkm_Dec2017_david.csv")
dt1 <- fread("data/Renyi_RNAseq_12292017/mes13_featurecounts_Dec2017_david.csv")
# dt1 <- fread("data/rna_seq/Renyi_12292017/mes13_featurecounts_Dec2017_david.csv")
dt1

# CHECK
dt1[Geneid == "Tnfrsf25",]
dt1[Geneid == "Fcer1g",]

# Keep only gene IDs and counts
dt1 <- dt1[, c("Geneid",
               "WJ1.dedup.bam",
               "WJ2.dedup.bam",
               "WJ3.dedup.bam",
               "WJ4.dedup.bam",
               "WJ5.dedup.bam",
               "WJ6.dedup.bam",
               "WJ7.dedup.bam")]
dt1
colnames(dt1) <- c("GeneNames",
                   trt.names)

dt1

# Remove genes with low counts----
summary(dt1[, -1])
tmp <- rowSums(dt1[, -1])
# Remove if total across 3 samples is no more than 10
dt1 <- droplevels(subset(dt1,
                         tmp > 20))
dt1
# 13,954 genes left, down from 24,421 genes

# Leave the 3 treatments only----
dt1 <- subset(dt1,
              select = c("GeneNames",
                         "LG",
                         "HG",
                         "TIIA_5uM"))
dt1

# Normalize data to fragments per million (FPM)----
mat <- data.frame(sample = colnames(dt1)[-1],
                  trt = colnames(dt1)[-1],
                  repl = factor(rep(1, ncol(dt1) - 1)))
mat

dtm <- as.matrix(dt1[, -1, with = FALSE])
rownames(dtm) <- dt1$GeneNames
head(dtm)

dds <- DESeqDataSetFromMatrix(countData = dtm, 
                              colData = mat,
                              ~ trt)
dds <- estimateSizeFactors(dds)
dds

# Fragments per million (FPM) normalization----
dt1.fpm <- data.table(GeneNames = dt1$GeneNames,
                      fpm(dds,
                          robust = FALSE))
colnames(dt1.fpm)[-1] <- paste(colnames(dt1.fpm)[-1],
                               "fpm",
                               sep = "_")
dt1.fpm

# DEGseq----
# a. (HG - LG)----
DEGexp(geneExpMatrix1 = dt1,
       geneCol1 = which(colnames(dt1) == "GeneNames"), 
       expCol1 = which(colnames(dt1) == "HG"), 
       groupLabel1 = "HG",
       
       geneExpMatrix2 = dt1,
       geneCol2 = which(colnames(dt1) == "GeneNames"), 
       expCol2 = which(colnames(dt1) == "LG"),
       groupLabel2 = "LG",
       
       foldChange = 2^0.3,
       qValue = 0.5,
       thresholdKind = 5, 
       rawCount = TRUE,
       normalMethod = "none",
       method = "MARS",
       outputDir = "tmp/hg_lg")

hg_lg <- fread("tmp/hg_lg/output_score.txt")
hg_lg
# CHECK:
hg_lg[hg_lg$GeneNames == "Nmu", ]
hg_lg[hg_lg$`Signature(q-value(Storey et al. 2003) < 0.5)`,]

# Write as CSV----
write.csv(hg_lg,
          file = "tmp/hg_lg/MES13_RNA_DEGseq_HG-LG.csv",
          row.names = FALSE)

# MA Plot----
hg_lg <- merge(hg_lg,
               dt1.fpm,
               by = "GeneNames")
hg_lg[, mu := (log2(HG_fpm + 1) + log2(LG_fpm + 1))/2]
hg_lg[, diff := log2(HG_fpm + 1) - log2(LG_fpm + 1)]

tiff(filename = "tmp/hg_lg/MES13_RNA_DEGseq_HG-LG_maplot.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(hg_lg$diff ~ hg_lg$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "MES13 FPM-Normalized Gene Expressions\nHG-LG, FDR < 0.5")
points(hg_lg$diff[hg_lg$`q-value(Storey et al. 2003)` < 0.5 & hg_lg$diff > 0] ~ hg_lg$mu[hg_lg$`q-value(Storey et al. 2003)` < 0.5 & hg_lg$diff > 0] ,
       pch = "x",
       col = "green")
points(hg_lg$diff[hg_lg$`q-value(Storey et al. 2003)` < 0.5 & hg_lg$diff < 0] ~ hg_lg$mu[hg_lg$`q-value(Storey et al. 2003)` < 0.5 & hg_lg$diff < 0] ,
       pch = "x",
       col = "red")
abline(h = c(-0.3, 0.3),
       lty = 2)
graphics.off()

# b. (TIIA - HG)----
DEGexp(geneExpMatrix1 = dt1,
       geneCol1 = which(colnames(dt1) == "GeneNames"), 
       expCol1 = which(colnames(dt1) == "TIIA_5uM"), 
       groupLabel1 = "TIIA_5uM",
       
       geneExpMatrix2 = dt1,
       geneCol2 = which(colnames(dt1) == "GeneNames"), 
       expCol2 = which(colnames(dt1) == "HG"),
       groupLabel2 = "HG",
       
       foldChange = 2^0.3,
       qValue = 0.5,
       thresholdKind = 5, 
       rawCount = TRUE,
       normalMethod = "none",
       method = "MARS",
       outputDir = "tmp/tiia_hg")

tiia_hg <- fread("tmp/tiia_hg/output_score.txt")
tiia_hg
#CHECK:
tiia_hg[tiia_hg$GeneNames == "Nmu"]
tiia_hg[tiia_hg$`Signature(q-value(Storey et al. 2003) < 0.5)`, ]

# Write as CSV----
write.csv(tiia_hg,
          file = "tmp/tiia_hg/MES13_RNA_DEGseq_TIIA-HG.csv",
          row.names = FALSE)

# MA Plot----
tiia_hg <- merge(tiia_hg,
                 dt1.fpm,
                 by = "GeneNames")
tiia_hg[, mu := (log2(TIIA_5uM_fpm + 1) + log2(HG_fpm + 1))/2]
tiia_hg[, diff := log2(TIIA_5uM_fpm + 1) - log2(HG_fpm + 1)]

tiff(filename = "tmp/tiia_hg/MES13_RNA_DEGseq_TIIA-HG_maplot.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(tiia_hg$diff ~ tiia_hg$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "MES13 FPM-Normalized Gene Expressions\nTIIA-HG, FDR < 0.5")
points(tiia_hg$diff[tiia_hg$`q-value(Storey et al. 2003)` < 0.5 & tiia_hg$diff > 0] ~ tiia_hg$mu[tiia_hg$`q-value(Storey et al. 2003)` < 0.5 & tiia_hg$diff > 0] ,
       pch = "x",
       col = "green")
points(tiia_hg$diff[tiia_hg$`q-value(Storey et al. 2003)` < 0.5 & tiia_hg$diff < 0] ~ tiia_hg$mu[tiia_hg$`q-value(Storey et al. 2003)` < 0.5 & tiia_hg$diff < 0] ,
       pch = "x",
       col = "red")
abline(h = c(-0.3, 0.3),
       lty = 2)
graphics.off()

# Venn diagram----
g1 <- hg_lg[`q-value(Storey et al. 2003)` < 0.5 & 
              `log2(Fold_change) normalized` > 0.3,]$GeneNames
# 263 genes
g2 <- hg_lg[`q-value(Storey et al. 2003)` < 0.5 & 
              `log2(Fold_change) normalized` < -0.3,]$GeneNames
# 207 genes

g3 <- tiia_hg[`q-value(Storey et al. 2003)` < 0.5 & 
                `log2(Fold_change) normalized` > 0.3,]$GeneNames
# 1,120 genes
g4 <- tiia_hg[`q-value(Storey et al. 2003)` < 0.5 & 
                `log2(Fold_change) normalized` < -0.3,]$GeneNames
# 1,393 genes

up.dn <- g1[g1 %in% g4]
# 124 genes
dn.up <- g2[g2 %in% g3]
# 89 genes

# Heatmap----
ll <- unique(c(up.dn,
               dn.up))

t1 <- merge(hg_lg[hg_lg$GeneNames %in% ll, 
                  c("GeneNames",
                    "log2(Fold_change) normalized")],
            tiia_hg[tiia_hg$GeneNames %in% ll, 
                    c("GeneNames",
                      "log2(Fold_change) normalized")],
            by = "GeneNames")
colnames(t1) <- c("Gene",
                  "HG-LG",
                  "TIIA-HG")
t1 <- t1[order(t1$`HG-LG`,
               decreasing = TRUE), ]
t1
t1[t1$Gene == "Nmu", ]
# Gene    HG-LG      TIIA-HG
# Nmu     4.128832   -0.8490297
write.csv(t1,
          file = "tmp/mes13_tiia_rnaseq_degseq_genes_q-0.5_log2-0.3.csv",
          row.names = FALSE)

ll <- melt.data.table(data = t1,
                      id.vars = 1,
                      measure.vars = 2:3,
                      variable.name = "Comparison",
                      value.name = "Gene Expression Diff")
ll$Comparison <- factor(ll$Comparison,
                        levels = c("TIIA-HG", 
                                   "HG-LG"))
lvls <- ll[ll$Comparison == "HG-LG", ]
ll$Gene <- factor(ll$Gene,
                  levels = lvls$Gene[order(lvls$`Gene Expression Diff`)])

# Add dendrogram----
# Sources: 
# https://stackoverflow.com/questions/43794870/plotting-a-clustered-heatmap-with-dendrograms-using-rs-plotly
# https://stats.stackexchange.com/questions/4062/how-to-plot-a-fan-polar-dendrogram-in-r
# https://cran.r-project.org/web/packages/ggdendro/vignettes/ggdendro.html
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
dt.dndr <- data.frame(t1[Gene %in% levels(ll$Gene), ])
rownames(dt.dndr) <- dt.dndr$Gene
dt.dndr <- dt.dndr[, -1]

# Compute distances between genes----
sampleDists <- dist(dt.dndr)
as.matrix(sampleDists)

# Make dendrogram data----
dhc <- as.dendrogram(hclust(d = sampleDists),
                     horiz = TRUE)
ddata <- dendro_data(dhc, 
                     type = "rectangle")
# Gene orger----
ddata$labels

# Segment data----
dtp1 <- segment(ddata)

# Hitmap data----
dtp2 <- ll
dtp2$Gene <- factor(dtp2$Gene,
                    levels = ddata$labels$label)
offset.size <- 10



p1 <- ggplot(data = dtp2) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_tile(aes(x =  as.numeric(Comparison),
                y = Gene, 
                fill = `Gene Expression Diff`),
            color = "white") +
  # geom_text(data = dtp2[Comparison == "HG-LG", ],
  #           aes(x = rep(2.25,
  #                       nlevels(Gene)),
  #               y = Gene,
  #               angle = 90 + seq(from = 0,
  #                                to = 360,
  #                                length.out = nlevels(Gene))[as.numeric(Gene)] +
  #                 offset.size,
  #               label = 1:nlevels(Gene)),
  #           hjust = 0,
  #           size = 6) +
  geom_text(data = dtp2[Comparison == "HG-LG", ],
            aes(x = rep(1.75,
                        nlevels(Gene)),
                y = Gene,
                angle = 90 + seq(from = 0,
                                 to = 360,
                                 length.out = nlevels(Gene))[as.numeric(Gene)] + 
                  offset.size,
                label = unique(Gene)),
            hjust = 0) +
  geom_text(data = dtp2[Gene == levels(dtp2$Gene)[1], ],
            aes(x = 1:nlevels(Comparison),
                y = rep(-offset.size,
                        nlevels(Comparison)),
                angle = 0,
                label = levels(Comparison)),
            hjust = 1,
            size = 8) +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "grey", 
                       midpoint = 0, 
                       name = "") +
  scale_y_discrete("",
                   expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 25),
        legend.direction = "horizontal",
        legend.key.width = unit(1.5, "in"),
        legend.key.height = unit(0.5, "in")) +
  geom_segment(data = dtp1,
               aes(x = -sqrt(y) + 0.5,
                   y = x, 
                   xend = -sqrt(yend) + 0.5,
                   yend = xend),
               size = 1) 
p1

tiff(filename = "tmp/MES13_RNA_DEGseq_TIIA-HG-LG_hitmap_with_phylo.tiff",
     height = 12,
     width = 12,
     units = 'in',
     res = 1200,
     compression = "lzw+p")
plot(p1)
graphics.off()

# CHECK:
dt1.fpm[GeneNames == "Iigp1"]
dtp2[Gene == "Iigp1"]
# Correct

sessionInfo()
sink()

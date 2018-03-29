# |----------------------------------------------------------------------------------|
# | Project: Skin UVB SKH1 mouse model treated with UA/SFN                           |
# | Script: Methyl-seq data analysis and visualization using DSS                     |
# | Coordinator: Ran Yin, Renyi Wu                                                   |
# | Author: Davit Sargsyan                                                           |
# | Created: 03/17/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_skin_uvb_dna_v1.txt")
date()

# NOTE: several packages, e.g. Rcpp, MASS, etc., might be deleted manually and reinstalled
# Workflow: https://www.bioconductor.org/help/workflows/rnaseqGene/
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
# biocLite("ChIPseeker")
# biocLite("DO.db")
# biocLite("GenomicRanges")
# biocLite("org.Mm.eg.db")
# biocLite("DSS")
# biocLite("bsseq")

require(data.table)
require(ggplot2)
require(knitr)

require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(DSS)

# Load and view raw counts (no annoation)----
dt01 <- fread("data/renyi_methylseq_02092018/combined_dmr_default__uvb-skin_renyi_02092018.csv")
dt01

# NOTE: there are 14 rows representing mitochondrial DNA
unique(dt01$chr)
dt01[dt01$chr == "chrM",]

# Load and view percent methylated (no annotation)----
dt02 <- fread("data/renyi_methylseq_02092018/combined_dmr_fraction_s32_uvb-skin_renyi_02092018.csv")
dt02

# Annotate----
peakAnno1 <- annotatePeak(peak = "data/renyi_methylseq_02092018/combined_dmr_default__uvb-skin_renyi_02092018.csv", 
                          tssRegion = c(-3000, 3000), 
                          TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                          annoDb = "org.Mm.eg.db")

head(peakAnno1@detailGenomicAnnotation)
t1 <- peakAnno1@annoStat
t1$Feature <- factor(t1$Feature,
                     levels = as.character(t1$Feature[order(t1$Frequency,
                                                            decreasing = TRUE)]))
t1
p1 <- ggplot(t1,
             aes(x = "",
                 y = Frequency,
                 fill = Feature)) +
  geom_bar(width = 1, 
           stat = "identity",
           color = "black") +
  coord_polar("y",
              start = 0,
              direction = -1) +
  scale_x_discrete("") +
  ggtitle("Annotation by Region (%)")
p1 

tiff(filename = "tmp/skin_uvb_anno_by_reg.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# # NOTE: mitochondrial DNA (chrM) not mapped.
# #       Also, chrX = 20 and chrY = 21. Compare below:
# unique(dt1$geneChr)
# dt1[dt1$geneChr == 20,]
# dt01[dt01$chr == "chrX",]
# 
# dt1[dt1$geneChr == 21,]
# dt01[dt01$chr == "chrY",]

# Make data table----
dt1 <- data.table(start = peakAnno1@anno@ranges@start,
                  as.data.frame(peakAnno1@anno@elementMetadata@listData))

# Remove unmapped regions
dt1 <- dt1[!is.na(dt1$SYMBOL == "NA"), ]
# Removed 12 rows

# Subset data; remove Ursolic Acid samles----
dt1 <- data.table(gene = dt1$SYMBOL,
                  anno = dt1$annotation,
                  geneId = dt1$geneId,
                  chr = dt1$geneChr,
                  pos = dt1$start,
                  reg = NA,
                  dt1[, CpG:X02w_SFN_1.X],
                  dt1[, X02w_UVB_0.N:X15w_SFN_1.X], 
                  dt1[, X15w_UVB_0.N:X25t_SFN_1.X],
                  dt1[, X25t_UVB_0.N:X25w_SFN_1.X],
                  dt1[, X25w_UVB_0.N:X25w_UVB_1.X],
                  geneName = dt1$GENENAME)
dt1

snames <- c("W2Ctrl1",
            "W2Ctrl2",
            "W2UVB1",
            "W2UVB2",
            "W2SFN1",
            "W2SFN2",
            
            "W15Ctrl1",
            "W15Ctrl2",
            "W15UVB1",
            "W15UVB2",
            "W15SFN1",
            "W15SFN2",
            
            "W25Ctrl1",
            "W25Ctrl2",
            "W25UVB1",
            "W25UVB2",
            "W25SFN1",
            "W25SFN2")
# "W25tCtrl1",
# "W25tCtrl2",
# "W25tUVB1",
# "W25tUVB2",
# "W25tSFN1",
# "W25tSFN2")

# Regions----
kable(data.table(table(substr(dt1$anno, 1, 9))))
  # |V1        |     N|
  # |:---------|-----:|
  # |3' UTR    |  4869|
  # |5' UTR    |   839|
  # |Distal In | 69042|
  # |Downstrea |  2939|
  # |Exon (uc0 | 12726|
  # |Intron (u | 60617|
  # |Promoter  | 86826|

# Separate Promoter, Body and Downstream----
dt1$reg <- as.character(dt1$anno)

# a. Promoter: up to 3kb upstream
dt1$reg[substr(dt1$anno, 
               1,
               8) == "Promoter"] <- "Promoter"

# b. Body: exons and introns
dt1$reg[substr(dt1$anno, 
               1, 
               4) %in% c("Exon",
                         "Intr")] <- "Body"

# c. Downstream: Distal Intergenic and  Downstream
dt1$reg[substr(dt1$anno, 
               1, 
               4) %in% c("Dist",
                         "Down")] <- "Downstream"

dt1$reg <- factor(dt1$reg,
                  levels = c("5' UTR",
                             "Promoter",
                             "Body",
                             "Downstream",
                             "3' UTR"))
kable(data.table(table(dt1$reg)))
  # |V1         |     N|
  # |:----------|-----:|
  # |5' UTR     |   839|
  # |Promoter   | 86826|
  # |Body       | 73343|
  # |Downstream | 71981|
  # |3' UTR     |  4869|

# CpG distribution and coverage----
summary(dt1$CpG)
p2 <- ggplot(dt1,
             aes(x = CpG)) +
  geom_histogram(color = "black",
                 fill = "grey",
                 binwidth = 5) +
  scale_x_continuous(name = "Number of CpG",
                     breaks = c(3, 
                                seq(from = 5,
                                    to = 60,
                                    by = 5))) +
  coord_cartesian(xlim=c(3, 60)) +
  scale_y_continuous(name = "Counts (x10,000)",
                     breaks = seq(from = 0, 
                                  to = 14*10^4,
                                  by = 10^4),
                     labels = seq(from = 0, 
                                  to = 14,
                                  by = 1)) +
  ggtitle("Distribution of DMR by Number of CpG")
p2

tiff(filename = "tmp/skin_uvb_CpG_hist.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

p3 <- ggplot(dt1,
             aes(x = CpG)) +
  facet_wrap(~ reg,
             scale = "free_y") +
  geom_histogram(color = "black",
                 fill = "grey",
                 binwidth = 5) +
  scale_x_continuous(name = "Number of CpG",
                     breaks = c(3, 
                                seq(from = 5,
                                    to = 60,
                                    by = 5))) +
  coord_cartesian(xlim=c(3, 60)) +
  scale_y_continuous(name = "Counts") +
  ggtitle("Distribution of DMR by Number of CpG and Region")
p3

tiff(filename = "tmp/skin_uvb_CpG_by_reg_hist.tiff",
     height = 8,
     width = 12,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p3)
graphics.off()

# Percent methylation----
tmp <- as.matrix(dt1[, X02w_CON_0.N:X25w_UVB_1.X])
# tmp <- apply(tmp,
#              MARGIN = 2,
#              function(a) {
#                a[is.na(a)] <- 0
#                return(a)
#              })
head(tmp)

dtN <- tmp[, seq(1,
                 ncol(tmp) - 1, 
                 by = 2)]
head(dtN)

dtX <- tmp[, seq(2,
                 ncol(tmp), 
                 by = 2)]
head(dtX)

pct <- dtX/dtN
colnames(pct) <- substr(colnames(pct),
                        1,
                        nchar(colnames(pct)) - 2)
head(pct)

# pct <- apply(pct,
#              MARGIN = 2,
#              function(a) {
#                a[is.nan(a)] <- 0
#                return(a)
#              })
head(pct)

# Remove rows with all zeros----
dim(pct[rowSums(pct) == 0, ])
dim(pct[is.na(rowSums(pct)), ])
dim(pct)
# 2,548 out of 237,858 such rows. Remove.

ndx.keep <- rowSums(pct) != 0
pct <- pct[ndx.keep, ]
dt1 <- dt1[ndx.keep, ]
dtN <- dtN[ndx.keep, ]
dtX <- dtX[ndx.keep, ]

# Hits per CpG average (i.e. vertical coverage)----
t1 <- apply(dtN,
            2,
            function(a) {
              return(round(a/dt1$CpG,
                           1))
            })

mu <- list()
for (i in 1:(ncol(t1)/2)) {
  x1 <- aggregate(x = t1[, 2*i - 1],
                       FUN = mean,
                       by = list(dt1$reg))
  x2 <- aggregate(x = t1[, 2*i],
                  FUN = mean,
                  by = list(dt1$reg))
  x3 <- merge(x1, x2, by = "Group.1")
  mu[[i]] <- data.table(reg = x3[, 1],
                        mu = (x3[, 2] + x3[, 3])/2)
}
names(mu) <- unique(substr(colnames(t1),
                           1,
                           8))
mu

# Average methylation per region per treatment/time
mumth <- list()
for (i in 1:(ncol(pct)/2)) {
  x1 <- aggregate(x = c(pct[, 2*i - 1],
                        pct[, 2*i]),
                  FUN = mean,
                  by = list(rep(dt1$reg, 2)))
  
  x2 <- aggregate(x = c(pct[, 2*i - 1],
                        pct[, 2*i]),
                  FUN = sd,
                  by = list(rep(dt1$reg, 2)))
  mumth[[i]] <- data.table(rep(colnames(pct)[[2*i]],
                               5),
                           merge(x1, 
                                 x2,
                                 by = "Group.1"))

  colnames(mumth[[i]]) <- c("trt",
                            "reg",
                            "mu",
                            "std")
}
mumth <- do.call("rbind",
                 mumth)
mumth
mumth$Time <- substr(mumth$trt, 1, 4)
mumth$Time <- mumth$Time <- factor(mumth$Time,
                                   levels = c("X02w",
                                              "X15w",
                                              "X25w",
                                              "X25t"),
                                   labels = c("Week 2",
                                              "Week 15",
                                              "Week 25",
                                              "Week 25 Tumor"))

mumth$Treatment <- substr(mumth$trt, 6, 8)
mumth$Treatment <- factor(mumth$Treatment,
                          levels = c("CON",
                                     "UVB",
                                     "SFN"),
                          labels = c("Control",
                                     "UVB",
                                     "SFN + UVB"))

mumth$`Methylation (%)` <- 100*mumth$mu

p1 <- ggplot(mumth,
             aes(x = reg,
                 y = `Methylation (%)`,
                 group = Treatment,
                 fill = Treatment)) +
  facet_wrap(~ Time) +
  geom_bar(position = position_dodge(),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous(limits = c(0, 60)) +
  ggtitle("Percent of Methylated CpG by Region") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p1
tiff(filename = "tmp/skin_uvb_avg_methyl_by_reg.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

p2 <- ggplot(mumth,
             aes(x = reg,
                 y = mu,
                 group = Treatment,
                 fill = Treatment)) +
  facet_wrap(~ Time) +
  geom_errorbar(aes(ymax = mu + std,
                    ymin = mu),
                width = 0.5,
                position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9),
             size = 1) +
  geom_bar(position = position_dodge(0.9),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous() +
  ggtitle("Proportion of Methylated CpG by Region") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p2
tiff(filename = "tmp/skin_uvb_avg_sd_methyl_by_reg.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Clustering----
dst <- dist(t(pct))
sampleDists <- as.matrix(dst)
sampleDists

# Make zeros less influential
diag(sampleDists) <- min(dst) - 5

tiff(filename = "tmp/skin_uvb_samples_cluster.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
heatmap(sampleDists,
        symm = TRUE,
        col = heat.colors(50))
graphics.off()

# PCA----
m1 <- prcomp(t(pct),
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
dt.scr$grp <- substr(colnames(pct),
                     6,
                     8)
dt.scr$grp <- factor(dt.scr$grp,
                     levels = c("CON",
                                "UVB",
                                "SFN"),
                     labels = c("Control",
                                "UVB",
                                "SFN + UVB"))
  
dt.scr$sample <- colnames(pct)
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

p1 <- ggplot(data = dt.rot,
             aes(x = PC1,
                 y = PC2)) +
  geom_point(data = dt.scr,
             aes(fill = grp),
             shape = 21,
             size = 3,
             alpha = 0.5) +
  geom_text(data = dt.scr,
            aes(x = PC1 + 40,
                y = PC2,
                label = dt.scr$sample),
            size = 2,
            hjust = 0.5) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  scale_fill_manual(name = "Treatment",
                    values = c("white", "red", "blue", "green")) +
  ggtitle("PCA of Percent Methylation") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10))
p1
tiff(filename = "tmp/skin_uvb_pca_plot.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Clustering, no tumor----
ndx.keep <- which(!(colnames(pct) %in% c("X25t_CON_0",
                                         "X25t_CON_1",
                                         "X25t_SFN_0",
                                         "X25t_SFN_1",
                                         "X25t_UVB_0",
                                         "X25t_UVB_1")))

dst <- dist(t(pct[, ndx.keep]))
sampleDists <- as.matrix(dst)
sampleDists

# Make zeros less influential
diag(sampleDists) <- min(dst) - 5

tiff(filename = "tmp/skin_uvb_samples_cluster_no_tumor.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
heatmap(sampleDists,
        symm = TRUE,
        col = heat.colors(50))
graphics.off()

# PCA, no tumor----
xx <- apply(t(pct[, ndx.keep]), 
            2,
            scale)
head(t(xx))

row.keep <- apply(t(xx),
                  1,
                  function(a) {
                    sum(is.na(a)) == 0
                  })

m1 <- prcomp(t(pct[row.keep, ndx.keep]),
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
dt.scr$grp <- substr(colnames(pct[row.keep, ndx.keep]),
                     6,
                     8)
dt.scr$grp <- factor(dt.scr$grp,
                     levels = c("CON",
                                "UVB",
                                "SFN"),
                     labels = c("Control",
                                "UVB",
                                "SFN + UVB"))

dt.scr$sample <- colnames(pct[row.keep, ndx.keep])
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

p1 <- ggplot(data = dt.rot,
             aes(x = PC1,
                 y = PC2)) +
  geom_point(data = dt.scr,
             aes(fill = grp),
             shape = 21,
             size = 3,
             alpha = 0.5) +
  geom_text(data = dt.scr,
            aes(x = PC1 + 40,
                y = PC2,
                label = dt.scr$sample),
            size = 2,
            hjust = 0.5) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  scale_fill_manual(name = "Treatment",
                    values = c("white", "red", "blue", "green")) +
  ggtitle("PCA of Percent Methylation") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10))
p1
tiff(filename = "tmp/skin_uvb_pca_plot_no_tumor.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Dispersion Shrinkage for Sequencing data (DSS)----
# This is based on Wald test for beta-binomial distribution.
# Source: https://www.bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.pdf
# The DM detection procedure implemented in DSS is based on a rigorous Wald test for betabinomial
# distributions. The test statistics depend on the biological variations (characterized
# by dispersion parameter) as well as the sequencing depth. An important part of the algorithm
# is the estimation of dispersion parameter, which is achieved through a shrinkage estimator
# based on a Bayesian hierarchical model [1]. An advantage of DSS is that the test can be
# performed even when there is no biological replicates. That’s because by smoothing, the
# neighboring CpG sites can be viewed as “pseudo-replicates", and the dispersion can still be
# estimated with reasonable precision.

# # One-way analysis (our example)----
# dtl <- list(data.table(dt1[, c("chr", "pos")],
#                        N = dt1$X02w_CON_0.N,
#                        X = dt1$X02w_CON_0.X),
#             data.table(dt1[, c("chr", "pos")],
#                        N = dt1$X02w_CON_1.N,
#                        X = dt1$X02w_CON_1.X),
#             data.table(dt1[, c("chr", "pos")],
#                        N = dt1$X02w_UVB_0.N,
#                        X = dt1$X02w_UVB_0.X),
#             data.table(dt1[, c("chr", "pos")],
#                        N = dt1$X02w_UVB_1.N,
#                        X = dt1$X02w_UVB_1.X))
# dtl
# renyi_methylseq_02092018
# BSobj <- makeBSseqData(dat = dtl,
#                        sampleNames = c("W2Ctrl1",
#                                        "W2Ctrl2",
#                                        "W2UVB1",
#                                        "W2UVB2"))
# BSobj
# 
# dmlTest <- DMLtest(BSobj,
#                    group1 = c("W2Ctrl1", 
#                               "W2Ctrl2"), 
#                    group2 = c("W2UVB1", 
#                               "W2UVB2"))

# Multi-factor analysis (treatment + time + interaction)----
dtl <- list(data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_CON_0.N,
                       X = dt1$X02w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_CON_1.N,
                       X = dt1$X02w_CON_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_UVB_0.N,
                       X = dt1$X02w_UVB_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_UVB_1.N,
                       X = dt1$X02w_UVB_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_SFN_0.N,
                       X = dt1$X02w_SFN_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_SFN_1.N,
                       X = dt1$X02w_SFN_1.X),

            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_CON_0.N,
                       X = dt1$X15w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_CON_1.N,
                       X = dt1$X15w_CON_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_UVB_0.N,
                       X = dt1$X15w_UVB_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_UVB_1.N,
                       X = dt1$X15w_UVB_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_SFN_0.N,
                       X = dt1$X15w_SFN_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_SFN_1.N,
                       X = dt1$X15w_SFN_1.X),

            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_CON_0.N,
                       X = dt1$X25w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_CON_1.N,
                       X = dt1$X25w_CON_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_UVB_0.N,
                       X = dt1$X25w_UVB_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_UVB_1.N,
                       X = dt1$X25w_UVB_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_SFN_0.N,
                       X = dt1$X25w_SFN_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_SFN_1.N,
                       X = dt1$X25w_SFN_1.X))

            # data.table(dt1[, c("chr", "pos")],
            #            N = dt1$X25t_CON_0.N,
            #            X = dt1$X25t_CON_0.X),
            # data.table(dt1[, c("chr", "pos")],
            #            N = dt1$X25t_CON_1.N,
            #            X = dt1$X25t_CON_1.X),
            # data.table(dt1[, c("chr", "pos")],
            #            N = dt1$X25t_UVB_0.N,
            #            X = dt1$X25t_UVB_0.X),
            # data.table(dt1[, c("chr", "pos")],
            #            N = dt1$X25t_UVB_1.N,
            #            X = dt1$X25t_UVB_1.X),
            # data.table(dt1[, c("chr", "pos")],
            #            N = dt1$X25t_SFN_0.N,
            #            X = dt1$X25t_SFN_0.X),
            # data.table(dt1[, c("chr", "pos")],
            #            N = dt1$X25t_SFN_1.N,
            #            X = dt1$X25t_SFN_1.X))
dtl

BSobj <- makeBSseqData(dat = dtl,
                       sampleNames = snames)
BSobj

design <- data.table(trt = rep(rep(c("Ctrl", "UVB", "SFN"),
                                   each = 2),
                               3),
                     time = rep(c(2, 15, 25),
                                each = 6))
design$trt <- factor(design$trt,
                     levels = unique(design$trt))

design$time <- factor(design$time,
                      levels = unique(design$time))
design

DMLfit <- DMLfit.multiFactor(BSobj = BSobj, 
                             design = design,
                             formula = ~ trt + time + trt*time)
summary(DMLfit)
head(DMLfit$fit$beta)

# Define contrasts----
contr <- data.table(design,
                    DMLfit$X)
contr

# Contrasts: UVB vs. Control at Week 2----
colnames(DMLfit$X)

# # Test treatment effect----
# DMLtest.UVB.Ctrl.Trt <- DMLtest.multiFactor(DMLfit, 
#                                            term = "trt")
# head(DMLtest.UVB.Ctrl.Trt)
# 
# # Test Time effect---
# DMLtest.UVB.Ctrl.Time <- DMLtest.multiFactor(DMLfit, 
#                                              term = "time")
# head(DMLtest.UVB.Ctrl.Time)
# 
# # UVB vs. Control----
# DMLtest.UVB.Ctrl <- DMLtest.multiFactor(DMLfit, 
#                                            coef = "trtUVB")
# head(DMLtest.UVB.Ctrl)

# Create a matrix of contrasts.
# Each column will correspond to a specific contrast.
# Each row will correspont to a specific column in the design matrix DMLfit$X

# We want the following comparisons (Kong email 03/19/2018):
# Well, let's go back to our discussion at the lab meeting last week. 
# So, we like to COMPARE the following, you may NOT need to include the "TIME" yet:
# We have 2 wk, 15 wk, 25 wk (adjacent normal epithelium, ANE), 
# 25 wk tumor (25T) and 25 wk tumor (whole skin - don't worry about this, may be confounding).
# Now, we want to compare
# control (no UV) vs UVB - 2 vs 15; 2 vs 25ANE, 2 vs 25T, may be also 15 vs 25ANE, 15 vs 25T
# UVB vs UVB + SFN - 2 vs 15; 2 vs 25ANE, 2 vs 25T, may be also 15 vs 25ANE, 15 vs 25T (same as above).
# Now another project down the road, we can compare the control (no UVB),
# 2 vs 25 vs 25 again effect. THIS will be later, NOT now.

# My comparisons: 
cnames <- c("Ctrl vs. UVB, Week 2",
            "Ctrl vs. UVB, Week 15",
            "Ctrl vs. UVB, Week 25",
            "UVB vs. SFN, Week 2",
            "UVB vs. SFN, Week 15",
            "UVB vs. SFN, Week 25",
            "Week 2 vs. Week 15, Control",
            "Week 2 vs. Week 25, Control",
            "Week 2 vs. Week 15, UVB",
            "Week 2 vs. Week 25, UVB",
            "Week 2 vs. Week 15, SFN",
            "Week 2 vs. Week 25, SFN",
            "Week 15 vs. Week 25, Control",
            "Week 15 vs. Week 25, UVB",
            "Week 15 vs. Week 25, SFN")

# Custom contrasts----
contr1 <- design
contr1$trttime <- paste(contr1$trt,
                        contr1$time,
                        sep = "")
xnames <- unique(contr1$trttime)
contr1[, c(xnames) := 0]
for (i in 4:12){
  contr1[, i] <- as.numeric(contr1$trttime == colnames(contr1)[i])
}
contr1
DMLfit1 <- DMLfit
DMLfit1$X <- contr1[, 4:12]
DMLfit1$X
DMLfit$X

mat1 <- matrix(0,
               nrow = length(xnames),
               ncol = length(cnames))
colnames(mat1) <- cnames
rownames(mat1) <- xnames
mat1

# 1. Ctrl vs. UVB, Week 2----
mat1[1, 1] <- -0.5
mat1[2, 1] <- 0.5
mat1
DMLtest.UVB.Ctrl.W2 <- DMLtest.multiFactor(DMLfit1,
                                           Contrast = mat1[,
                                                           1,
                                                           drop = FALSE])

DMLtest.UVB.Ctrl.W2 <- DMLtest.multiFactor(DMLfit1,
                                           coef = "UVB2")
head(DMLtest.UVB.Ctrl.W2)

DMLtest.UVB.Ctrl.W2$chr <- as.numeric(as.character(DMLtest.UVB.Ctrl.W2$chr))

DMLtest.UVB.Ctrl.W2 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                   pct[, c("X02w_CON_0",
                                                           "X02w_CON_1",
                                                           "X02w_UVB_0",
                                                           "X02w_UVB_1")]),
                                        DMLtest.UVB.Ctrl.W2,
                                        by = c("chr",
                                               "pos")))
DMLtest.UVB.Ctrl.W2$mu.ctrl2w <- (DMLtest.UVB.Ctrl.W2$X02w_CON_0 + DMLtest.UVB.Ctrl.W2$X02w_CON_0)/2
DMLtest.UVB.Ctrl.W2$mu.uvb2w <- (DMLtest.UVB.Ctrl.W2$X02w_UVB_0 + DMLtest.UVB.Ctrl.W2$X02w_UVB_1)/2

DMLtest.UVB.Ctrl.W2[5, ]

# Mean vs. difference----
DMLtest.UVB.Ctrl.W2$mu <- (DMLtest.UVB.Ctrl.W2$mu.ctrl2w + DMLtest.UVB.Ctrl.W2$mu.uvb2w)/2
DMLtest.UVB.Ctrl.W2$diff <- (DMLtest.UVB.Ctrl.W2$mu.ctrl2w - DMLtest.UVB.Ctrl.W2$mu.uvb2w)

dtp1 <- DMLtest.UVB.Ctrl.W2[!is.na(stat), ]

dtp1[dtp1$diff > 0.5]

plot(dtp1$diff ~ dtp1$mu,
     pch = ".")
points(dtp1$diff[dtp1$fdrs < 0.01] ~ dtp1$mu[dtp1$fdrs < 0.01],
       pch = "x",
       col = "red")

# 1. SFN vs. Control, Week 2----
mat1[1, 2] <- -0.5
mat1[3, 2] <- 0.5
mat1
DMLtest.SFN.Ctrl.W2 <- DMLtest.multiFactor(DMLfit1,
                                           Contrast = mat1[,
                                                           2,
                                                           drop = FALSE])
head(DMLtest.SFN.Ctrl.W2)

DMLtest.SFN.Ctrl.W2$chr <- as.numeric(as.character(DMLtest.SFN.Ctrl.W2$chr))

DMLtest.SFN.Ctrl.W2 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                   pct[, c("X02w_CON_0",
                                                           "X02w_CON_1",
                                                           "X02w_UVB_0",
                                                           "X02w_UVB_1")]),
                                        DMLtest.SFN.Ctrl.W2,
                                        by = c("chr",
                                               "pos")))
DMLtest.SFN.Ctrl.W2$mu.ctrl2w <- (DMLtest.SFN.Ctrl.W2$X02w_CON_0 + DMLtest.SFN.Ctrl.W2$X02w_CON_0)/2
DMLtest.SFN.Ctrl.W2$mu.uvb2w <- (DMLtest.SFN.Ctrl.W2$X02w_UVB_0 + DMLtest.SFN.Ctrl.W2$X02w_UVB_1)/2

DMLtest.SFN.Ctrl.W2[5, ]

# Mean vs. difference----
DMLtest.SFN.Ctrl.W2$mu <- (DMLtest.SFN.Ctrl.W2$mu.ctrl2w + DMLtest.SFN.Ctrl.W2$mu.uvb2w)/2
DMLtest.SFN.Ctrl.W2$diff <- (DMLtest.SFN.Ctrl.W2$mu.ctrl2w - DMLtest.SFN.Ctrl.W2$mu.uvb2w)

dtp1 <- DMLtest.SFN.Ctrl.W2[!is.na(stat), ]

plot(dtp1$diff ~ dtp1$mu,
     pch = ".")
points(dtp1$diff[dtp1$fdrs < 0.01] ~ dtp1$mu[dtp1$fdrs < 0.01],
       pch = "x",
       col = "red")

snames <- c("W2Ctrl1",
            "W2Ctrl2",
            "W2UVB1",
            "W2UVB2",
            "W2SFN1",
            "W2SFN2")

# Multi-factor analysis (treatment + time + interaction)----
dtl <- list(data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_CON_0.N,
                       X = dt1$X02w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_CON_1.N,
                       X = dt1$X02w_CON_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_UVB_0.N,
                       X = dt1$X02w_UVB_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_UVB_1.N,
                       X = dt1$X02w_UVB_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_SFN_0.N,
                       X = dt1$X02w_SFN_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_SFN_1.N,
                       X = dt1$X02w_SFN_1.X))
dtl

BSobj <- makeBSseqData(dat = dtl,
                       sampleNames = snames)
BSobj

design <- data.table(trt = rep(rep(c("UVB", "Ctrl", "SFN"),
                                   each = 2),
                               1))
design$trt <- factor(design$trt,
                     levels = unique(design$trt))
design

DMLfit <- DMLfit.multiFactor(BSobj = BSobj, 
                             design = design,
                             formula = ~ trt)
summary(DMLfit)
head(DMLfit$fit$beta)
colnames(DMLfit$X)

# a. Control vs. UVB at Week 2----
DMLtest.Ctrl.UVB.W2 <- DMLtest.multiFactor(DMLfit,
                                           coef = "trtCtrl")
head(DMLtest.Ctrl.UVB.W2)

DMLtest.Ctrl.UVB.W2$chr <- as.numeric(as.character(DMLtest.Ctrl.UVB.W2$chr))

DMLtest.Ctrl.UVB.W2 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                   pct[, c("X02w_CON_0",
                                                           "X02w_CON_1",
                                                           "X02w_UVB_0",
                                                           "X02w_UVB_1")]),
                                        DMLtest.Ctrl.UVB.W2,
                                        by = c("chr",
                                               "pos")))
DMLtest.Ctrl.UVB.W2$mu.ctrl2w <- (DMLtest.Ctrl.UVB.W2$X02w_CON_0 + DMLtest.Ctrl.UVB.W2$X02w_CON_1)/2
DMLtest.Ctrl.UVB.W2$mu.uvb2w <- (DMLtest.Ctrl.UVB.W2$X02w_UVB_0 + DMLtest.Ctrl.UVB.W2$X02w_UVB_1)/2

# Mean vs. difference----
DMLtest.Ctrl.UVB.W2$mu <- (DMLtest.Ctrl.UVB.W2$mu.ctrl2w + DMLtest.Ctrl.UVB.W2$mu.uvb2w)/2
DMLtest.Ctrl.UVB.W2$diff <- (DMLtest.Ctrl.UVB.W2$mu.ctrl2w - DMLtest.Ctrl.UVB.W2$mu.uvb2w)

dtp1 <- DMLtest.Ctrl.UVB.W2[!is.na(stat), ]

tiff(filename = "tmp/skin_uvb_maplot_ctrl_uvb_w2.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylation in Control vs. UVB\nat Two Weeks, FDR < 0.01")
points(dtp1$diff[dtp1$fdrs < 0.01] ~ dtp1$mu[dtp1$fdrs < 0.01],
       pch = "x",
       col = "red")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# b. SFN vs. UVB at Week 2----
DMLtest.SFN.UVB.W2 <- DMLtest.multiFactor(DMLfit,
                                          coef = "trtSFN")
head(DMLtest.SFN.UVB.W2)

DMLtest.SFN.UVB.W2$chr <- as.numeric(as.character(DMLtest.SFN.UVB.W2$chr))

DMLtest.SFN.UVB.W2 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                  pct[, c("X02w_SFN_0",
                                                          "X02w_SFN_1",
                                                          "X02w_UVB_0",
                                                          "X02w_UVB_1")]),
                                       DMLtest.SFN.UVB.W2,
                                       by = c("chr",
                                              "pos")))
DMLtest.SFN.UVB.W2$mu.sfn2w <- (DMLtest.SFN.UVB.W2$X02w_SFN_0 + DMLtest.SFN.UVB.W2$X02w_SFN_1)/2
DMLtest.SFN.UVB.W2$mu.uvb2w <- (DMLtest.SFN.UVB.W2$X02w_UVB_0 + DMLtest.SFN.UVB.W2$X02w_UVB_1)/2

# Mean vs. difference----
DMLtest.SFN.UVB.W2$mu <- (DMLtest.SFN.UVB.W2$mu.sfn2w + DMLtest.SFN.UVB.W2$mu.uvb2w)/2
DMLtest.SFN.UVB.W2$diff <- (DMLtest.SFN.UVB.W2$mu.sfn2w - DMLtest.SFN.UVB.W2$mu.uvb2w)

dtp1 <- DMLtest.SFN.UVB.W2# [!is.na(stat), ]

tiff(filename = "tmp/skin_uvb_maplot_sfn_uvb_w2.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylation in SFN vs. UVB\nat Two Weeks, FDR < 0.01")
points(dtp1$diff[dtp1$fdrs < 0.01] ~ dtp1$mu[dtp1$fdrs < 0.01],
       pch = "x",
       col = "red")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# sessionInfo()
# sink()
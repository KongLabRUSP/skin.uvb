# |--------------------------------------------------------|
# | Project: Skin UVB SKH1 moue model treated with UA/SFN. |
# | Study ID:                                              |
# | Scientist: Yuqing (Anne) Yang                          |
# | Data Analysis: Ran Yin, Renyi Wu, Davit Sargsyan       |
# | Created: 12/12/2017                                    |
# |--------------------------------------------------------|
# Header----
require(data.table)
require(knitr)
require(VennDiagram)
require(gridExtra)

# Load data----
ls.files <- dir("data")
ls.files

dt2w <- fread("data/week2 gene_exp.diff")
dt2w

dt15w <- fread("data/week15 gene_exp.diff")
dt15w

dt25w <- fread("data/week25 gene_exp.diff")
dt25w

dtf<- fread("data/WeekFinal gene_exp.diff")
dtf <- dtf[, 1:14]

# dtf[, 8:13] <- apply(dtf[, 8:13],
#                      MARGIN = 2,
#                      FUN = "as.numeric")
dtf$value_1 <- as.numeric(dtf$value_1)
dtf$value_2 <- as.numeric(dtf$value_2)
dtf$`log2(fold_change)` <- as.numeric(dtf$`log2(fold_change)`)
dtf$q_value <- as.numeric(dtf$q_value )
dtf

# Filter by log2 change and q-value----
unique(dt2w$sample_1)
unique(dt2w$sample_2)

# Week 2----
g2w.ctrl.uvb <- dt2w$gene[dt2w$sample_1 == "Control" &
                            dt2w$sample_2 == "UVB" &
                            abs(dt2w$`log2(fold_change)`) >= 2 & 
                            dt2w$q_value <= 0.1]
g2w.ctrl.uvb

g2w.uvb.ua <- dt2w$gene[dt2w$sample_1 == "UVB" &
                            dt2w$sample_2 == "UA" &
                            abs(dt2w$`log2(fold_change)`) >= 2 & 
                            dt2w$q_value <= 0.1]
g2w.uvb.ua

g2w.uvb.sfn <- dt2w$gene[dt2w$sample_1 == "UVB" &
                          dt2w$sample_2 == "SFN" &
                          abs(dt2w$`log2(fold_change)`) >= 2 & 
                          dt2w$q_value <= 0.1]
g2w.uvb.sfn

# Week 15----
g15w.ctrl.uvb <- dt15w$gene[dt15w$sample_1 == "Control" &
                              dt15w$sample_2 == "UVB" &
                              abs(dt15w$`log2(fold_change)`) >= 2 & 
                              dt15w$q_value <= 0.1]
g15w.ctrl.uvb

g15w.uvb.ua <- dt15w$gene[dt15w$sample_1 == "UVB" &
                            dt15w$sample_2 == "UA" &
                            abs(dt15w$`log2(fold_change)`) >= 2 & 
                            dt15w$q_value <= 0.1]
g15w.uvb.ua

g15w.uvb.sfn <- dt15w$gene[dt15w$sample_1 == "UVB" &
                            dt15w$sample_2 == "SFN" &
                            abs(dt15w$`log2(fold_change)`) >= 2 & 
                            dt15w$q_value <= 0.1]
g15w.uvb.sfn

# Week 25----
g25w.ctrl.uvb <- dt25w$gene[dt25w$sample_1 == "Control" &
                              dt25w$sample_2 == "UVB" &
                              abs(dt25w$`log2(fold_change)`) >= 2 & 
                              dt25w$q_value <= 0.1]
g25w.ctrl.uvb

g25w.uvb.ua <- dt25w$gene[dt25w$sample_1 == "UVB" &
                            dt25w$sample_2 == "UA" &
                            abs(dt25w$`log2(fold_change)`) >= 2 & 
                            dt25w$q_value <= 0.1]
g25w.uvb.ua

g25w.uvb.sfn <- dt25w$gene[dt25w$sample_1 == "UVB" &
                            dt25w$sample_2 == "SFN" &
                            abs(dt25w$`log2(fold_change)`) >= 2 & 
                            dt25w$q_value <= 0.1]
g25w.uvb.sfn

# Final Week----
gf.ctrl.uvb <- dtf$gene[dtf$sample_1 == "Control" &
                              dtf$sample_2 == "UVB" &
                              abs(dtf$`log2(fold_change)`) >= 2 & 
                              dtf$q_value <= 0.1]
dtf[dtf$sample_1 == "Control" &
           dtf$sample_2 == "UVB" &
           abs(dtf$`log2(fold_change)`) >= 2, ]
gf.ctrl.uvb

# Remove NAs----
sum(is.na(gf.ctrl.uvb))
gf.ctrl.uvb <- gf.ctrl.uvb[!is.na(gf.ctrl.uvb)]

gf.uvb.ua <- dtf$gene[dtf$sample_1 == "UVB" &
                            dtf$sample_2 == "UA" &
                            abs(dtf$`log2(fold_change)`) >= 2 & 
                            dtf$q_value <= 0.1]
gf.uvb.ua

gf.uvb.sfn <- dtf$gene[dtf$sample_1 == "UVB" &
                        dtf$sample_2 == "SFN" &
                        abs(dtf$`log2(fold_change)`) >= 2 & 
                        dtf$q_value <= 0.1]
gf.uvb.sfn

# g2w <- dt2w$gene[abs(dt2w$`log2(fold_change)`) >= 2 & 
#                    dt2w$q_value <= 0.05]
# g15w <- dt15w$gene[abs(dt15w$`log2(fold_change)`) >= 2 & 
#                    dt15w$q_value <= 0.05]
# g25w <- dt25w$gene[abs(dt25w$`log2(fold_change)`) >= 2 & 
#                    dt25w$q_value <= 0.05]
# gf <- dtf$gene[abs(dtf$`log2(fold_change)`) >= 2 & 
#                    dtf$q_value <= 0.05]
# 
# # Genes in common----
# gene.shared <- unique(g2w[g2w %in% g15w])
# gene.shared <- unique(gene.shared[gene.shared %in% g25w])
# gene.shared <- unique(gene.shared[gene.shared %in% gf])
# gene.shared

# Create a table of shared genes----
genes <- data.table(all.signif = unique(c(g2w.ctrl.uvb,
                                          g2w.uvb.ua,
                                          g2w.uvb.sfn,
                                          
                                          g15w.ctrl.uvb,
                                          g15w.uvb.ua,
                                          g15w.uvb.sfn,
                                          
                                          g25w.ctrl.uvb,
                                          g25w.uvb.ua,
                                          g25w.uvb.sfn,
                                          
                                          gf.ctrl.uvb,
                                          gf.uvb.ua,
                                          gf.uvb.sfn)))

genes$g2w.ctrl.uvb <- genes$all.signif %in% g2w.ctrl.uvb
genes$g2w.uvb.ua <- genes$all.signif %in% g2w.uvb.ua
genes$g2w.uvb.sfn <- genes$all.signif %in% g2w.uvb.sfn

genes$g15w.ctrl.uvb <- genes$all.signif %in% g15w.ctrl.uvb
genes$g15w.uvb.ua <- genes$all.signif %in% g15w.uvb.ua
genes$g15w.uvb.sfn <- genes$all.signif %in% g15w.uvb.sfn

genes$g25w.ctrl.uvb <- genes$all.signif %in% g25w.ctrl.uvb
genes$g25w.uvb.ua <- genes$all.signif %in% g25w.uvb.ua
genes$g25w.uvb.sfn <- genes$all.signif %in% g25w.uvb.sfn

genes$gf.ctrl.uvb <- genes$all.signif %in% gf.ctrl.uvb
genes$gf.uvb.ua <- genes$all.signif %in% gf.uvb.ua
genes$gf.uvb.sfn <- genes$all.signif %in% gf.uvb.sfn

genes

  

# Number of significant genes----
t1 <- data.table(Week = rep(c("Week2",
                              "Week15",
                              "Week25",
                              "Final"),
                            each = 3),
                 Comparison = rep(c("Control (No UVB) vs. UVB",
                                    "UVB vs. UVB+UA",
                                    "UVB vs. UVB+SFN"),
                                  4),
                 N = colSums(genes[, -1]))

t1$Week <- factor(t1$Week,
                  levels = unique(t1$Week))
t1$Comparison <- factor(t1$Comparison,
                        levels = unique(t1$Comparison))

t1 <- dcast.data.table(t1,
                       Week ~ Comparison)
t1
kable(t1)
  # |Week   | Control (No UVB) vs. UVB| UVB vs. UVB+UA| UVB vs. UVB+SFN|
  # |:------|------------------------:|--------------:|---------------:|
  # |Week2  |                      549|            102|              48|
  # |Week15 |                      379|              8|              12|
  # |Week25 |                      525|             21|              13|
  # |Final  |                      140|             13|               9|

# Venn diagrams----
# a. Common genes by comparisons
ctrl.uvb <- unique(c(g2w.ctrl.uvb,
                     g15w.ctrl.uvb,
                     g25w.ctrl.uvb,
                     gf.ctrl.uvb))
uvb.ua <- unique(c(g2w.uvb.ua,
                   g15w.uvb.ua,
                   g25w.uvb.ua,
                   gf.uvb.ua))
uvb.sfn <- unique(c(g2w.uvb.sfn,
                   g15w.uvb.sfn,
                   g25w.uvb.sfn,
                   gf.uvb.sfn))

sum(ctrl.uvb %in% uvb.ua)
# 116 genes in common

sum(ctrl.uvb %in% uvb.sfn)
# 49 genes in common

genes.compar <- data.table(all.signif = genes$all.signif)
genes.compar$ctrl.uvb <- genes$all.signif %in% ctrl.uvb
genes.compar$uvb.ua <- genes$all.signif %in% uvb.ua
genes.compar$uvb.sfn <- genes$all.signif %in% uvb.sfn
genes.compar

# Save gene list to CSV----
tmp <- genes.compar
colnames(tmp) <- c("All Significant Genes",
                   
                   "Control (No UVB) vs. UVB",
                   "UVB vs. UVB+UA",
                   "UVB vs. UVB+SFN")
tmp
write.csv(tmp,
          file = "tmp/genes_compar_log2.2_qval.0.1.csv",
          row.names = FALSE)

p1 <- venn.diagram(x = list(`Control (No UVB) vs. UVB` = ctrl.uvb,
                            `UVB vs. UVB+UA` = uvb.ua,
                            `UVB vs. UVB+SFN` = uvb.sfn),
                   col = c("red",
                           "green",
                           "blue"),
                   compression = "lzw+p",
                   cat.just = list(c(0, 8),
                                   c(1, 8),
                                   c(0.5, -7)),
                   filename = "tmp/venn1.tiff",
                   units = 'in',
                   height = 6,
                   width = 6,
                   main = "Significant genes, at least 4-fold change, q-value <= 0.05")

# b. Common genes by time
w2 <- unique(c(g2w.ctrl.uvb,
               g2w.uvb.ua,
               g2w.uvb.sfn))
w15 <- unique(c(g15w.ctrl.uvb,
                g15w.uvb.ua,
                g15w.uvb.sfn))
w25 <- unique(c(g25w.ctrl.uvb,
                g25w.uvb.ua,
                g25w.uvb.sfn))
wf <- unique(c(gf.ctrl.uvb,
               gf.uvb.ua,
               gf.uvb.sfn))

genes.time <- data.table(all.signif = genes$all.signif)
genes.time$w2 <- genes$all.signif %in% w2
genes.time$w15 <- genes$all.signif %in% w15
genes.time$w25 <- genes$all.signif %in% w25
genes.time$wf <- genes$all.signif %in% wf
genes.time

# Save gene list to CSV----
tmp <- genes.time
colnames(tmp) <- c("All Significant Genes",
                   "Week2",
                   "Week15",
                   "Week25",
                   "Final Week")
tmp
write.csv(tmp,
          file = "tmp/genes_time_log2.2_qval.0.1.csv",
          row.names = FALSE)

p2 <- venn.diagram(x = list(`Week2` = w2,
                            `Week15` = w15,
                            `Week25` = w25,
                            `Final` = wf),
                   col = c("red",
                           "green",
                           "blue",
                           "black"),
                   compression = "lzw+p",
                   filename = "tmp/venn2.tiff",
                   units = 'in',
                   height = 5,
                   width = 5,
                   main = "Significant genes, at least 4-fold change, q-value <= 0.05")

# Write raw data to CSV, significant genes ONLY
write.csv(dt2w[dt2w$gene %in% genes$all.signif, ],
          file = "tmp/dt2w_signif_genes.csv")
write.csv(dt15w[dt2w$gene %in% genes$all.signif, ],
          file = "tmp/dt15w_signif_genes.csv")
write.csv(dt2w[dt25w$gene %in% genes$all.signif, ],
          file = "tmp/dt25w_signif_genes.csv")
write.csv(dt2w[dtf$gene %in% genes$all.signif, ],
          file = "tmp/dtf_signif_genes.csv")






# c. Common genes by time within each comparison
p3 <- venn.diagram(x = list(`Control (No UVB) vs. UVB` = g2w.ctrl.uvb,
                            `UVB vs. UVB+UA` = g2w.uvb.ua,
                            `UVB vs. UVB+SFN` = g2w.uvb.sfn),
                   col = c("red",
                           "green",
                           "blue"),
                   compression = "lzw+p",
                   filename = "tmp/venn3_w2.tiff",
                   cat.just = list(c(0, 8),
                                   c(1, 8),
                                   c(0.5, -7)),
                   units = 'in',
                   height = 8,
                   width = 8,
                   main = "Significant genes, at least 4-fold change\nq-value <= 0.05, Week2 Only")

# Write raw data to CSV, significant genes ONLY
write.csv(dt2w[dt2w$gene %in% genes$all.signif, ],
          file = "tmp/dt2w_signif_genes.csv")
write.csv(dt15w[dt2w$gene %in% genes$all.signif, ],
          file = "tmp/dt15w_signif_genes.csv")
write.csv(dt2w[dt25w$gene %in% genes$all.signif, ],
          file = "tmp/dt25w_signif_genes.csv")
write.csv(dt2w[dtf$gene %in% genes$all.signif, ],
          file = "tmp/dtf_signif_genes.csv")

# sink()
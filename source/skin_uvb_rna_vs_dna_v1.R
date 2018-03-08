require(data.table)
require(ggplot2)

# Load RNA data----
dtrna <- fread("data/pc_vs_nc.csv")
dtrna

dtrna.dn <- dtrna[padj < 0.01 & log2FoldChange < -1, ]
dtrna.dn
# 274 significant genes

dtrna.up <- dtrna[padj < 0.01 & log2FoldChange > 1, ]
dtrna.up
# 188 genes

# Load DNA data----
dtdna <- fread("data/ran_dmrfinder_methylseq_data/wk2_all_results.csv_anno.csv")
dtdna


dtdna <- droplevels(dtdna[, c("gene",
                              "feature",
                              "CpG",
                              "distance",
                              "Control..UVB.diff",
                              "Control..UVB.pval")])
dtdna$padjust <- p.adjust(dtdna$Control..UVB.pval,
                       method = "fdr")

# Keep genes found in RNA data
dtdna.dn <- dtdna[gene %in% dtrna.dn$gene, ]

dtt <- merge(dtrna.dn, 
             dtdna.dn,
             by = "gene")
dtt$reg <- factor(substr(dtt$feature, 1, 3))

# dtt <- dtt[dtt$padjust < 0.05]
dtt <- dtt[dtt$Control..UVB.diff >= 0.1]

p1 <- ggplot(dtt,
             aes(x = Control..UVB.diff,
                 y = log2FoldChange,
                 colour = reg)) + 
  facet_wrap(~gene) +
  geom_point(size = 3)
  # geom_abline(slope = -5, 
  #             intercept = -3)
p1

tiff(filename = "tmp/skin_uvb_ua.vs.pc.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

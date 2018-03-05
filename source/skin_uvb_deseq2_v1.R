# Load packages----
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

require(data.table)
require(DESeq2)
require(BiocParallel)
require(glmmADMB)

set.seed(1000)

# Load data----
dt1 <- fread("data/Kong Lab-UVB skin study.csv")
dt1

dt2w <- fread("data/week2 gene_exp.diff")


### Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution (R package *DESeq2*)


# Bild a DESeq2 model with treatment*time interaction----
dds <- DESeqDataSetFromMatrix(dt1, 
                              grp,
                              ~ trt + time + trt:time)

# If all samples contain zeros, geometric means cannot be
# estimated. Change default 'type = "ratio"' to 'type = "iterate"'.
# Type '?DESeq2::estimateSizeFactors' for more details.
# dss <- estimateSizeFactors(dds, 
#                            type = "iterate")

#ALTERNATEVELY: exclude zeros from geometric mean calculation----
geoMeans <- apply(X = counts(dds), 
                  MARGIN = 1,
                  FUN = gm_mean)

dds <- estimateSizeFactors(object = dds,
                           geoMeans = geoMeans)

# Set cores for parallel processing of DESeq----
snowparam <- SnowParam(workers = snowWorkers(), 
                       type = "SOCK")
register(snowparam, 
         default = TRUE)

# Run DESeq----
dds <- DESeq(dds,
             fitType = "local",
             parallel = TRUE)
resultsNames(dds)

# Contrasts----
# a. Treatment B/Treatment A at Time 1----
resAvsBtime1 <- results(dds,
                        name = "trt_TrtB_vs_TrtA")

# b. Treatment B/Treatment A at Time 2----
resAvsBtime2 <- results(dds,
                        contrast = list(c("trt_TrtB_vs_TrtA" ,
                                          "trtTrtB.timeTime2")))

# ALTERNATIVE GROUPING----
# Create a mixed-level factor----
grp$trt_time <- factor(paste(grp$trt, 
                             grp$time,
                             sep = "_"))

dds.mix <- DESeqDataSetFromMatrix(dt1, 
                                  grp,
                                  ~ trt_time)

geoMeans.mix <- apply(X = counts(dds.mix), 
                      MARGIN = 1,
                      FUN = gm_mean)

dds.mix <- estimateSizeFactors(object = dds.mix,
                               geoMeans = geoMeans.mix)

# Run DESeq----
dds.mix <- DESeq(dds.mix,
                 fitType = "local",
                 parallel = TRUE)
resultsNames(dds.mix)

# Contrasts----
# a. Treatment B/Treatment A at Time 1----
resAvsBtime1.mix <- results(dds.mix,
                            contrast = c("trt_time",
                                         "TrtB_Time1",
                                         "TrtA_Time1"))
# Compare
data.table(do.call("cbind", 
                   resAvsBtime1@listData))
data.table(do.call("cbind", 
                   resAvsBtime1.mix@listData))

# b. Treatment B/Treatment A at Time 2----
resAvsBtime2.mix <- results(dds.mix,
                            contrast = c("trt_time",
                                         "TrtB_Time2",
                                         "TrtA_Time2"))
# Compare
data.table(do.call("cbind", 
                   resAvsBtime2@listData))
data.table(do.call("cbind", 
                   resAvsBtime2.mix@listData))
```


```{r glmmADMB, echo=TRUE, fig.height=5, fig.width=10, message=FALSE, warning=FALSE}
# Merge counts with grouping
dt.l <- data.table(grp, t(dt1))
colnames(dt.l)[5:ncol(dt.l)] <- rownames(dt1)

# Melt data table
dt.l <- melt.data.table(data = dt.l,
                        id.vars = colnames(grp),
                        measure.vars = 5:ncol(dt.l),
                        variable.name = "gene",
                        value.name = "count")

# Model
tmp <- droplevels(subset(dt.l,
                         dt.l$gene == "Gene1"))

out <- list()
for (i in 1:nlevels(dt2.KO$otu)) {
  # Subset the data to i-th OTU
  tmp.i <- droplevels(subset(dt2.KO,
                             otu == levels(otu)[i]))
  
  # Try building the model
  res <- try({
    out[[i]] <- glmmadmb(Count ~ Treatment*Week + (Treatment | MouseID) + offset(sample.gmean),
                         family = "nbinom",
                         save.dir = "tmp/glmmadmb_out",
                         # zeroInflation = TRUE,
                         data = tmp.i)
  })
  if (class(res)[1] == "try-error") {          
    out[[i]] <- NA
  }
  print(paste("Processing OTU", i, "..."))
}
names(out) <- levels(dt2.KO$otu)
out
```{r echo=FALSE}
# Subset data based on condition pair specified in RUVReport.Rmd
cond1="{{cond1}}"
cond2="{{cond2}}"
KK="{{K}}"
KK <- as.numeric(KK)
suffix <- paste0(cond2, "_vs_", cond1)
dir.create(suffix)
dir.create(paste0(suffix,"/plots"))
FigPaths <- paste0("./", suffix, "/plots/")

# sampleData clean-up ----------------------------------------------------------
### Only get the row.names (sample) of investigating condition. 
sampleData <- sampleParent[c(which(sampleParent$condition==cond1),
                             which(sampleParent$condition==cond2)), ]
### Only get the colnames (covariate) of investigating variances. 
sampleData <- sampleData[, which(covParent$take_nottake==1)]
### Relevel condition so that it's TREATMENT vs CONTROL 
sampleData$condition <- droplevels(sampleData$condition)
sampleData$condition <- factor(sampleData$condition, levels=c(cond1, cond2))

# countData clean-up -----------------------------------------------------------
countData <- countParent[, row.names(sampleData)]

# covData clean-up -------------------------------------------------------------
covData <- covParent[which(covParent$take_nottake==1),]

# Define colors to produce diagnosic plot below 
colors <- brewer.pal(8, "Set1")
colors <- colors[sampleData$condition]
```

# Differential Expression Results: `r cond2` vs `r cond1`

```{r include=FALSE}
if (all.equal(colnames(sampleData), as.character(covData$covariate)) != TRUE) {
  error1 <- "Error: Check sampleData and covData. They do not have the same covariate strings."
  print(colnames(sampleData))
  print(covData$covariate)}

if (all.equal(colnames(countData), row.names(sampleData)) != TRUE) {error2 <- "Error: Check countData and sampleData. They do not have the same column names and row names"
  } else {
    # Make DESeq object 
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=sampleData, design= ~ condition)
    dds$condition <- factor(dds$condition, levels=c(cond1, cond2))
    dds <- DESeq(dds, parallel=TRUE)
  } #-End else
```

```{r}
if(exists("error1")) { print(error1)}
if(exists("error2")) { print(error2)}
```
## DESeq2 First-Pass
DESeq2 was used to perform data normalization and differential expression analysis with an adjusted p-value threshold of 0.05 on each set of raw expression measures. This is the first pass to create a set of genes whose expressions were not differentially expressed between 2 investigated groups. These genes, refered to as empirical control genes, are then used for downstream RUVseq analysis (Removal of Unwanted Variation). 
```{r include=FALSE}
# Write results of DESeq2 first pass to a text file 
res <- results(dds)
res <- lfcShrink(dds, res=res, coef = resultsNames(dds)[2])
res <- as.data.frame(res)
res <- res[order(res$padj), ]
write.table(res, file=paste0(suffix, "/", "FirstPass_deSeq2_", cond2, "_vs_", cond1, ".txt"), quote=FALSE, sep="\t", col.names=NA)
res1 <- res
# Get negative genes for empirical control. In this case we want those genes whose padj values are more than 0.75
empirical <- row.names(res)[which(res$padj > 0.75)]
```
##### The number of empirical control genes used in this experiment is: `r length(empirical)`
```{r}
# If the number of genes are more than 100, print the last 100 genes (i.e. the highest padj values)
# if (length(empirical) > 100) {print(tail(empirical, 100))} else {print(empirical)}
filename <- paste0(date, "EmpiricalControlGenes", cond2, "_vs_", cond1, ".txt")
writeLines(empirical, sep="\n", con=filename)

# Run RUVSeq (RUVg)
set <- betweenLaneNormalization(as.matrix(countData), which="upper")
setG <- RUVg(set, empirical, k=KK)
```

## Data Visualization Before and After RUVSeq
### RLE Plot
The boxplots of relative log expression (RLE = log-ratio of read count to median read count across sample)
```{r}
# Define y-axis limit
x <- log2(as.matrix(countData + 1))
ylim <- max(apply(x, 2, mean)) 
```

```{r ,fig.path=FigPaths}
#par(mfrow=c(3,1))
plotRLE(as.matrix(countData), outline = FALSE, ylim = c(-ylim, ylim), col = colors, las = 2, cex.axis=0.8, ylab = "Relative Log Expression", main="Raw Counts")
plotRLE(set, outline = FALSE, ylim = c(-ylim, ylim), col = colors, las = 2, cex.axis=0.8, ylab = "Relative Log Expression", main="UQ Normalization")
plotRLE(as.matrix(counts(dds, normalize=TRUE)), outline = FALSE, ylim = c(-ylim, ylim), col = colors, las = 2, cex.axis=0.8, ylab = "Relative Log Expression", main="DESeq2 Normalization")
plotRLE(setG$normalizedCounts, outline = FALSE, ylim = c(-ylim, ylim), col = colors, las = 2, cex.axis=0.8, ylab = "Relative Log Expression", main=paste0("RUVg Corrected - K=", KK))
```

### Covariate to PC correlation 
The plots depict whether data are separated by biological or techinical variations. High correlation value (-log10(p-value)) indicates an association between variation in one PC and investigated covariate.
```{r include=FALSE}
COUNT <- list(as.matrix(countData), set, counts(dds, normalize=TRUE), setG$normalizedCounts)
names(COUNT) <- c("NN", "UQ", "DESeq", "RUV")
PCov <- lapply(COUNT, run_pca)
```

```{r ,fig.path=FigPaths, fig.width=5, fig.height=4}
plotPCAVar(PCov[[1]], "Raw Data")
plotPCAVar(PCov[[2]], "Upper Quartile")
plotPCAVar(PCov[[3]], "DESeq2")
plotPCAVar(PCov[[4]], "RUVg")
```

## DESeq2 Second-Pass

DESeq2 was used to perform data normalization and differential expression analysis with an adjusted p-value threshold of 0.05 on each set of raw expression measures. This is the second pass after RUVSeq. 

```{r}
sampleDataW <-cbind(sampleData, setG$W)
ddsW <- DESeqDataSetFromMatrix(countData = countData, colData = sampleDataW, design = ~ condition)
ddsW$condition <- factor(ddsW$condition, levels = c(cond1, cond2))

### compile formula for DESeq2 based on number of K:
      form_string <- ""
      for (w in 1:KK) {
        if (w == KK) {
          form_string <- paste(form_string, "W_", w, sep = "")
        } else {
          form_string <- paste(form_string, "W_", w, " + ", sep = "")
        } # end else
      } # end for
      
      design(ddsW) <- formula(paste("~", form_string, " + condition", sep = ""))

ddsW <- DESeq(ddsW, parallel = TRUE)
res <- results(ddsW)
coef <- paste0("condition_", cond2, "_vs_", cond1)
res <- lfcShrink(ddsW, res=res, coef = coef)
res.summary <- res
# Write results of DESeq2 second pass to a text file 
res <- res[order(res$padj), ]
res <- as.data.frame(res)
write.table(res, file=paste0(suffix, "/", "RUV_deSeq2_", cond2, "_vs_", cond1, ".txt"), quote=FALSE, sep="\t", col.names=NA)
res2 <- res
```

### Summary 
```{r echo=FALSE, results='asis'}
summaryString <- capture.output(summary(res.summary))

for (s in summaryString) {
  if(nchar(trimws(s)) > 0 && substr(trimws(s), 1, 1) != '['){
    cat(paste0('`', s, '` \n\n'))
  }
}
```
---
title: "Analysis for JC Virus bulk RNA-seq"
author: "Nguyen PT Huynh"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_document:
    highlight: tango
    toc: true
    toc_depth: 3
    fig_width: 5
    toc_float: true
abstract: "Script to generate figures and differential analysis results for Cui et al."
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(warning=FALSE, message=FALSE)
```

# Packages and Data Input

```{r}
# Read in libraries
library("DESeq2")
library("RColorBrewer")
library("BiocStyle")
library("knitr")
library("EDASeq")
library("VennDiagram")
library("RUVSeq")
library("ggplot2")
library("enrichR")
library("pheatmap")
library("BiocParallel")
register(MulticoreParam(8))
# color to use for heatmap and pca plots 
main_color <- scales::hue_pal()(3)
colors <- RColorBrewer::brewer.pal(8, "Set1")[3:5]
main_pal <- colorRampPalette(c("navy", "red4", "yellow"))(300)
```


```{r inputData}
# Differential Expression Test
cond <- list(test1=c("JC4", "CTR"), test2=c("JC14", "CTR"), test3=c("JC14", "JC4"))
K <- 2
# Read in data. There should be 3 txt files located under subfolder ./data
countParent <- read.table("./data/countData.txt", sep="\t", header=TRUE, row.names=1)
covParent <- read.table("./data/covData.txt", sep="\t", header=TRUE)
covParent$covariate <- as.character(covParent$covariate)
sampleParent <- read.table("./data/sampleData.txt", sep="\t", header=TRUE, row.names=1)
row.names(sampleParent) <- paste0("X", row.names(sampleParent))
# Filter data to only retain genes that have at least 2 counts in one third of the data set.
athird <- round(dim(countParent)[2]/3)
countParent <- countParent[which(apply(countParent, 1, sum) > (2*athird)), ]
survivingGenes <- dim(countParent)[1]
```

# Differential Expression Analysis with Removal of Unwated Variants

```{r results="hide"}
# Interactively run through all test conditions, computes and renders all html
src=lapply(seq_along(cond), function(i) {
  cur_cond=names(cond)[[i]]
  print(cond[[cur_cond]][2])
  print(cond[[cur_cond]][1])
  knit_expand(file="./RuvChild.Rmd", cond1=cond[[cur_cond]][2], cond2=cond[[cur_cond]][1], K=K)
})

knit_res=paste(knit_child(text=unlist(src)), collapse="\n")
```

`r paste(knitr::knit(text=knit_res))`

# Data with Batch Correction by RUV 

```{r}
# RUVg with empirical control genes from RuvChild.Rmd above
empirical <- c(readLines("EmpiricalControlGenesJC14_vs_CTR.txt"), readLines("EmpiricalControlGenesJC14_vs_JC4.txt"), readLines("EmpiricalControlGenesJC4_vs_CTR.txt"))
empirical <- unique(empirical)
countData <- read.table("./data/countData.txt")
colnames(countData) <- colnames(count)
set <- betweenLaneNormalization(as.matrix(countData), which="upper")
setG <- RUVg(set, empirical, k=2)
# Write out corrected-data for visualization
write.table(setG$normalizedCounts, "countData_RUV_K2.txt", row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
```

# Software

```{r echo=FALSE}
sessionInfo()
```

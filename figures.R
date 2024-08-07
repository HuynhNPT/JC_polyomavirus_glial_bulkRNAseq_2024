library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)
library(enrichplot)
# Encode95 table of gene name and equivalent ensemble ID for plyr::mapvalues
name2id <- read.table("hs95_name2ID.txt", stringsAsFactors = FALSE, header = TRUE)
############################# PCA ########################################
# countData_RUV_K2.txt was generated from main.md
# This is the corrected data after unwanted variation has been removed from our dataset
ruv_corrected <- read.table("RUV/data/countData_RUV_K2.txt", stringsAsFactors = FALSE)
pca <- prcomp(t(log2(ruv_corrected+1)))
pcaCoord <- pca$x

eigs <- pca$sdev^2
percentage <- c(eigs[1]/sum(eigs), eigs[2]/sum(eigs))
percentage <-  round(percentage*100, 1)
percentage <- paste0("(", percentage, "%)")
pcaCoord <- as.data.frame(pcaCoord)
pcaCoord$condition <- factor(rep(c("Day 0","Day 4","Day 14"),times = 3), levels = c("Day 0","Day 4","Day 14"))
pcaCoord$set <- rep(c("Set1", "Set2", "Set3"), each = 3)
p_pca <- ggplot(pcaCoord, aes(PC1, PC2, color=condition)) + geom_point() + theme_bw() 
p_pca <- p_pca + geom_text(aes(PC1, PC2, label=set), nudge_x = 1, nudge_y = 5)
p_pca <- p_pca + xlab(paste0("PC1 ", percentage[1])) + ylab(paste0("PC2 ", percentage[2]))
p_pca <- p_pca + guides(color = guide_legend(title = "dpi")) + scale_color_manual(values = rev(RColorBrewer::brewer.pal(3, "Set1"))) + ggtitle("PCA Plot")
############################# Venn Diagram ########################################
# Differential expression results from main.md
dres <- list(early = "RUV/JC4_vs_CTR/RUV_deSeq2_JC4_vs_CTR.txt",
             late = "RUV/JC14_vs_JC4/RUV_deSeq2_JC14_vs_JC4.txt",
             overall = "RUV/JC14_vs_CTR/RUV_deSeq2_JC14_vs_CTR.txt")
dres <- lapply(dres, function(x) read.table(x, stringsAsFactors = FALSE))
dres_filt <- lapply(dres, function(x) x[which(x$log2FoldChange != 0 & x$padj < 0.01),])
up_list <- lapply(dres_filt, function(x) row.names(x)[which(x$log2FoldChange > 0)])
dn_list <- lapply(dres_filt, function(x) row.names(x)[which(x$log2FoldChange < 0)])

vp_up <- venn.diagram(up_list, filename = NULL)
grid.newpage(); grid.draw(vp_up)

vp_dn <- venn.diagram(dn_list, filename = NULL)
grid.newpage(); grid.draw(vp_dn)

junk <- list.files(".", "VennDiagram")
file.remove(junk)
############################# Heatmaps ########################################
all_dys <- c(unlist(up_list), unlist(dn_list)) %>% unique()
filt_corrected <- log2(ruv_corrected[row.names(ruv_corrected) %in% all_dys,]+1)
annot_color <- list(condition = rev(c(RColorBrewer::brewer.pal(3, "Set1"))))
names(annot_color$condition) <- c("Day 0", "Day 4", "Day 14")
p_heat <- pheatmap(as.matrix(filt_corrected), show_rownames = FALSE, scale = 'row', treeheight_row = 0, treeheight_col = 20, color = viridis::viridis_pal(option = "B")(300), annotation_col = pcaCoord[,'condition', drop = FALSE], annotation_colors = annot_color, silent = TRUE, main = "Expression Heatmap of\nDE Genes - log2(count+1)")[[4]]

############################# Cluster Profiler ########################################
# Since late response (14 vs 4) and overall response (14 vs 0) are almost identical, we focus the gene ontology analysis on early repsonse (4 vs 0) and late (14 vs 4)
up_gid <- lapply(up_list[1:2], function(x) plyr::mapvalues(x, name2id$external_gene_name, name2id$ensembl_gene_id, warn_missing = FALSE))
dn_gid <- lapply(dn_list[1:2], function(x) plyr::mapvalues(x, name2id$external_gene_name, name2id$ensembl_gene_id, warn_missing = FALSE))
dys_gid <- c(up_gid, dn_gid)
names(dys_gid) <- c('up_early', 'up_late', 'dn_early', 'dn_late')

ck <- compareCluster(geneCluster = dys_gid,
                     fun = enrichGO,
                     OrgDb = org.Hs.eg.db,
                     keyType = 'ENSEMBL',
                     qvalueCutoff = 0.01,
                     ont='BP', readable = TRUE)

ck <- pairwise_termsim(ck)
p_cnet <- cnetplot(ck, cex_label_gene = 0.01, cex_label_category = 2) + scale_fill_manual(values = c("orange", "firebrick", "navy")) + ggtitle("Similarity in Gene Membership by GO Terms")
p_cnet2 <- cnetplot(ck, cex_label_gene = 1, cex_label_category = 1.2) + scale_fill_manual(values = c("orange", "firebrick", "navy")) + ggtitle("Similarity in Gene Membership by GO Terms with Gene Names")
p_semsim <- emapplot(ck) + scale_fill_manual(values = c("orange", "firebrick", "navy")) + ggtitle("Semantic Similarity of GO Terms")
p_egos <- dotplot(ck) + ggtitle("Enriched GO Terms\nfor Dysregulated Genes")
############################# Heatmaps of pathways ########################################
# Group terms of interest based on cnet
tid_list <- list(DNA_replication = c("DNA replication"),
                 cell_cycle = c("chromosome segregation", "mitotic nuclear division", "nuclear chromosome segregation", "nuclear division", "sister chromatid segregation"),
                 interferon_signaling = ck@compareClusterResult$Description[grep("interferon", ck@compareClusterResult$Description)],
                 sterol_biosynthetic_process = c("sterol biosynthetic process"),
                 synaptic_support = c("axon development", "regulation of nervous system development", "modulation of chemical synaptic transmission", 'regulation of trans-synaptic signaling'))
# Check for typs 
lapply(tid_list, function(x) table(x %in% ck@compareClusterResult$Description))
#
gid_mat <- list()
p_heat_top <- list()
for (i in 1:length(tid_list)) {
  tmp <- ck@compareClusterResult$geneID[ck@compareClusterResult$Description %in% tid_list[[i]]]
  tmp <- strsplit(tmp, "\\/") %>% unlist() %>% unique()
  # Top 10 genes late response since the early should be pretty similar to the late in terms of cell cycle
  tmp_dres <- dres_filt$late[row.names(dres_filt$late) %in% tmp, ] %>% top_n(-10, padj)
  gid_mat[[i]] <- log2(ruv_corrected[row.names(ruv_corrected) %in% row.names(tmp_dres),]) + 1
  # Wrap data for pheatmap 
  sets <- rep(c("Set1", "Set2", "Set3"), times = 3)
  days <- rep(c("D00", "D04", "D14"), each = 3)
  gid_mat[[i]] <- gid_mat[[i]][, paste0(sets, "_", days)]
  
  p_heat_top[[i]] <- pheatmap(gid_mat[[i]], cluster_cols = FALSE, show_colnames = FALSE, treeheight_row = 0, scale = "row", annotation_col = pcaCoord[, "condition", drop = FALSE], annotation_colors = annot_color, color = viridis::viridis_pal(option = "B")(300), silent = TRUE, annotation_legend = FALSE, main = names(tid_list)[i])[[4]]
}
# 

P_heat_terms <- plot_grid(plotlist = p_heat_top, ncol = 5)
############################# ALL TOGETHER NOW ########################################
label_size = 20
P1 <- plot_grid(p_pca, NULL, p_heat, ncol = 3, rel_widths = c(1.6, 1, 1), labels = c("A", "B", "C"), label_size = label_size)
P2 <- plot_grid(p_egos, p_cnet, ncol = 2, rel_widths = c(1.7, 2), labels = c("D", "E"), label_size = label_size)

plot_grid(P1, P2,NULL, P_heat_terms, ncol = 1, rel_heights = c(0.75, 1.5,0.1, 0.7), labels = c("", "", "F"), label_size = label_size)

plot_grid(p_cnet2, ncol = 1, labels = c("A"))
#------------------------------------------------------------------------------
# 11 - ClassifyTCGAERTumors
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(dplyr)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(ggplot2)
TCGA_BRCA_ATAC_UQ_se <- readRDS("rds/TCGA_BRCA_ATAC_UQ_se.rds")
rownames(TCGA_BRCA_ATAC_UQ_se) <- mcols(TCGA_BRCA_ATAC_UQ_se)$name

#----- deconvolute peaks -----#
#cell type specific regions
cellTypeDAR <- import.bed("ref/cellTypeDAR_hg38_scATAC.bed")
n <- names(table(cellTypeDAR$name))[-3]
cellTypeDAR <- lapply(n, function(x){
  res <- cellTypeDAR[cellTypeDAR$name == x]
  res <- subsetByOverlaps(res, cellTypeDAR[cellTypeDAR$name == "Epithelial"], invert = T)
  return(res)
})
names(cellTypeDAR) <- n
cellTypeDAR_total <- GRangesList(cellTypeDAR) %>% unlist %>% unique

#annotate peaks
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- ChIPseeker::annotatePeak(rowRanges(TCGA_BRCA_ATAC_UQ_se), tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
mcols(TCGA_BRCA_ATAC_UQ_se) <- mcols(as.GRanges(peakAnno))

#cancer distal peaks
TCGA_BRCA_ATAC_UQ_se_dst <- TCGA_BRCA_ATAC_UQ_se[grep("1kb", mcols(TCGA_BRCA_ATAC_UQ_se)$annotation, invert = T), ]
TCGA_BRCA_ATAC_UQ_se_cdp <- TCGA_BRCA_ATAC_UQ_se_dst[subsetByOverlaps(rowRanges(TCGA_BRCA_ATAC_UQ_se_dst), cellTypeDAR_total, invert = T) %>% names,]
##peaknumber :: 215920 > 160420 > 150039

#----- deconvolute peaks -----#
#ER+ tumors
idx <- which(TCGA_BRCA_ATAC_UQ_se_cdp$BRCA_scmod2 == "ER+/HER2- High Prolif" | TCGA_BRCA_ATAC_UQ_se_cdp$BRCA_scmod2 == "ER+/HER2- Low Prolif")
TCGA_BRCA_ATAC_ER_se_cdp <- TCGA_BRCA_ATAC_UQ_se_cdp[,idx]

VariablePeaks <- data.frame(peaks = mcols(TCGA_BRCA_ATAC_ER_se_cdp)$name, var = rowVars(assays(TCGA_BRCA_ATAC_ER_se_cdp)$normcounts)) %>% 
                 .[order(.$var, decreasing = T), ] %>% .[c(1:50000),]
hc1 <- hclust(dist(t(assays(TCGA_BRCA_ATAC_ER_se_cdp)$normcounts[VariablePeaks$peaks,])), method = "ward.D2")
pdf("output/Plots/S6_TCGA_hclust_ERcdp.pdf", width = 5, height = 10)
plot(hc1)
dev.off()

hc_class <- cutree(hc1, k = 3)
hc_class[hc_class == 1] <- "A"
hc_class[hc_class == 2] <- "ILC"
hc_class[hc_class == 3] <- "B"
ha1 <- HeatmapAnnotation(Cluster = hc_class, 
                         PAM50 = TCGA_BRCA_ATAC_ER_se_cdp$BRCA_pam50, 
                         Pathol = TCGA_BRCA_ATAC_ER_se_cdp$histological_type,
                         annotation_legend_param = c(Cluster = list(title = "Cluster"),
                                                     PAM50 = list(title = "PAM50", at = levels(TCGA_BRCA_ATAC_ER_se_cdp$BRCA_pam50)),
                                                     Pathol = list(title = "Pathol", at = levels(TCGA_BRCA_ATAC_ER_se_cdp$histological_type))),
                         col = list(Cluster = c("A" = "#1B9E77", "B" = "#D95F02", "ILC" = "#E7298A"),
                                    PAM50 =  c("Basal" = "red", "LumA" = "blue", "LumB" = "skyblue", "Her2" = "orange", "Normal" = "darkgray"),
                                    Pathol = c("Infiltrating Ductal Carcinoma" = "#FFFFB3", "Infiltrating Lobular Carcinoma" =  "#BEBADA", "Mixed Histology (please specify)" = "#FDB462", "Other, specify" = "black"))
                                    )
col_fun2 = colorRamp2(c(0,1,2,3,4,5), viridis(6, option = "D"))
fh = function(x) hclust(dist(x), method="ward.D2")

mtx <- assays(TCGA_BRCA_ATAC_ER_se_cdp)$normcounts[VariablePeaks$peaks,]
peak_km <- kmeans(mtx, centers = 5, iter.max = 100)
ht1 <- Heatmap(mtx, 
               name = "log2(Normalized counts)", cluster_columns = fh, cluster_rows = F, 
               show_row_names = F, show_column_names = F,
               row_split = peak_km$cluster, column_split = hc_class,
               top_annotation = ha1, col = col_fun2)
p1 <- draw(ht1)
pdf("output/Plots/S6_TCGA_HeatmapMostVariable50k_PtClass.pdf", width = 8, height = 6)
p1
dev.off()

lapply(c(1:5), function(i){
  g <- rowRanges(TCGA_BRCA_ATAC_ER_se_cdp)[names(peak_km$cluster)[peak_km$cluster == i]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/TCGAkmeansPeaks/kmeansPeaks_", i, "_output.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#PCA
pca1 <- prcomp(t(assays(TCGA_BRCA_ATAC_UQ_se_cdp)$normcounts))
df <- pca1$x[, c(1:2)] %>% data.frame
df$cluster <- as.character(TCGA_BRCA_ATAC_UQ_se_cdp$BRCA_scmod2)
df$cluster[rownames(df) %in% names(hc_class)] <- hc_class
df$pam50 <- TCGA_BRCA_ATAC_UQ_se_cdp$BRCA_pam50
df$cluster[is.na(df$cluster)] <- "NoData"
df$pam50[is.na(df$pam50)] <- "NoData"
p2 <- ggplot(df, aes(x = PC1, y = PC2, color = cluster)) + geom_point() + ArchR::theme_ArchR() +
      scale_color_manual(values = c("A" = "#1B9E77", "B" = "#D95F02", "ILC" = "#E7298A", "ER-/HER2-" = "red2", "HER2+" = "orange", "NoData" = "black"))
p3 <- ggplot(df, aes(x = PC1, y = PC2, color = pam50)) + geom_point() + ArchR::theme_ArchR() +
      scale_color_manual(values = c("Basal" = "red", "LumA" = "blue", "LumB" = "skyblue", "Her2" = "orange", "Normal" = "darkgray", "NoData" = "black"))
pdf("output/Plots/S6_TCGA_PCA_classify.pdf", width = 4, height = 4.5)
p2
p3
dev.off()

TCGA_BRCA_ATAC_ER_se_cdp$CA <- hc_class
saveRDS(TCGA_BRCA_ATAC_ER_se_cdp, "rds/TCGA_BRCA_ATAC_ER_se_cdp.rds")

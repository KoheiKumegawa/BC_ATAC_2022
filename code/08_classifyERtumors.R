#------------------------------------------------------------------------------
# 08 - classify ER+ tumors
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
se <- readRDS("rds/se2.rds")
se_cdp <- readRDS("rds/se_cdp.rds")

#----- Commonly accessible peaks -----#
peakStats <- data.frame(row.names = mcols(se_cdp)$name, 
                        median = rowMedians(assays(se_cdp)$normcounts), 
                        var =    rowVars(assays(se_cdp)$normcounts))
p1 <- ggplot(peakStats, aes(x = var, y = median)) + geom_hex(bins = 100) + ArchR::theme_ArchR() + 
  scale_fill_viridis_c(trans = "log") + labs(x = "Variance", y = "Median accessibility") + xlim(0, 2.5) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.5)
commonPeaks <- rownames(peakStats[which(peakStats$median > 3 & peakStats$var < 0.5), ])
g <- rowRanges(se_cdp)[commonPeaks]
d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
write.table(d, "output/output_bed/commonPeaks/commonPeaks_output.bed", row.names = F, col.names = F, quote = F, sep = "\t")

pdf("output/Plots/S4_commonPeaksScatter.pdf", width = 4, height = 5)
p1
dev.off()

GREAT_Com <- data.table::fread("output/great_go/commonPeaks/commonPeaks.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_Com$mlog10FDR <- -log10(GREAT_Com$FDR)
p2 <- ggplot(GREAT_Com, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 
pdf("output/Plots/S4_GREAT_commonPeaks.pdf", width = 10, height = 3)
p2
dev.off()

motif_DF <- data.table::fread("output/homer_motif/commonPeaks/commonPeaks/knownResults.txt", header = F, skip = 1)[c(1:10), c(1,3)] %>% data.frame %>% `colnames<-`(., c("Motif", "P"))
motif_DF$mlog10P <- as.numeric(gsub("1e-", "", motif_DF$P))
p3 <- ggplot(motif_DF, aes(x = mlog10P, y = reorder(Motif, mlog10P))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
      labs(x = "Motif", y = "-log10(P-value)")
pdf("output/Plots/S4_Motif_commonPeaks.pdf", width = 10, height = 3)
p3
dev.off()

#----- classify ER+ tumors -----#
se_er <- se[, which(se$Subtype == "L")]
se_cdp_er <- se_cdp[, which(se_cdp$Subtype == "L")]

VariablePeaks <- data.frame(peaks = mcols(se_cdp_er)$name, var = rowVars(assays(se_cdp_er)$normcounts)) %>% 
  .[order(.$var, decreasing = T), ] %>% .[c(1:50000),]
hc1 <- hclust(dist(t(assays(se_er)$normcounts)), method = "ward.D2")
hc2 <- hclust(dist(t(assays(se_cdp_er)$normcounts)), method = "ward.D2")
VariablePeaks <- data.frame(peaks = mcols(se_cdp_er)$name, var = rowVars(assays(se_cdp_er)$normcounts)) %>% 
                 .[order(.$var, decreasing = T), ] %>% 
                 .[c(1:50000),]
hc3 <- hclust(dist(t(assays(se_cdp_er)$normcounts[VariablePeaks$peaks,])), method = "ward.D2")
pdf("output/Plots/S4_hclust.pdf", width = 7, height = 4)
plot(hc1)
plot(hc2)
plot(hc3)
dev.off()

hc_class <- cutree(hc3, k = 3)
hc_class[hc_class == 1] <- "A"
hc_class[hc_class == 2] <- "B"
hc_class[hc_class == 3] <- "C"

col_fun1 <- colorRamp2(c(0,100), c("white", "brown"))
ha1 <- HeatmapAnnotation(Cluster = hc_class, ER = colData(se_cdp_er)$ER, PGR = colData(se_cdp_er)$PgR, HER2 = colData(se_cdp_er)$HER2, Ki67 = colData(se_cdp_er)$Ki67, Pathol = colData(se_cdp_er)$HistologicalType, Pr_Re = colData(se_cdp_er)$Primary_Relapse,
                         annotation_legend_param = c(Cluster = list(title = "Cluster"),
                                                     ER = list(title = "ER", at = levels(colData(se_cdp_er)$ER)),
                                                     PGR = list(title = "PGR", at = levels(colData(se_cdp_er)$PgR)), 
                                                     HER2 = list(title = "PGR", at = levels(colData(se_cdp_er)$HER2)),
                                                     Ki67 = list(title = "PAM50", at = c(0,25,50,75,100), col_fun = col_fun1),
                                                     Pathol = list(title = "Pathol", at = levels(colData(se_cdp_er)$HistologicalType)),
                                                     Pr_Re = list(title = "Primary_Recurrent", at = levels(colData(se_cdp_er)$Primary_Relapse))),
                         col = list(Cluster = brewer.pal(3, "Dark2") %>% `names<-`(., c("A", "B", "C")),
                                    ER = brewer.pal(length(levels(colData(se_cdp_er)$ER)), "Blues") %>% `names<-`(., levels(colData(se_cdp_er)$ER)),
                                    PGR = brewer.pal(length(levels(colData(se_cdp_er)$PgR)), "Greens") %>% `names<-`(., levels(colData(se_cdp_er)$PgR)),
                                    HER2 = brewer.pal(length(levels(colData(se_cdp_er)$HER2)), "Oranges") %>% `names<-`(., levels(colData(se_cdp_er)$HER2)),
                                    Ki67 = col_fun1,
                                    Pathol = brewer.pal(length(levels(colData(se_cdp_er)$HistologicalType)), "Set3") %>% `names<-`(., levels(colData(se_cdp_er)$HistologicalType)),
                                    Pr_Re = c("gray80", "black") %>% `names<-`(., c("P", "R")))
)
col_fun2 = colorRamp2(c(0,1,2,3,4,5,6), viridis(7, option = "D"))
fh = function(x) hclust(dist(x), method="ward.D2")

mtx <- assays(se_cdp_er)$normcounts[VariablePeaks$peaks,]
peak_km <- kmeans(mtx, centers = 5, iter.max = 100)
ht1 <- Heatmap(mtx, 
               name = "log2(Normalized counts)", cluster_columns = fh, cluster_rows = F, show_row_names = F,
               row_split = peak_km$cluster, column_split = hc_class,
               top_annotation = ha1, col = col_fun2)
p4 <- draw(ht1)
pdf("output/Plots/S4_HeatmapMostVariable50k_PtClass.pdf", width = 8, height = 8)
p4
dev.off()

lapply(c(1:5), function(i){
  g <- rowRanges(se_cdp_er)[names(peak_km$cluster)[peak_km$cluster == i]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/kmeansPeaks/kmeansPeaks_", i, "_output.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#PCA
pca1 <- prcomp(t(assays(se_cdp)$normcounts))
df <- pca1$x[, c(1:2)] %>% data.frame
df$cluster <- as.character(se_cdp$Subtype)
df$cluster[rownames(df) %in% names(hc_class)] <- hc_class
p5 <- ggplot(df, aes(x = PC1, y = PC2, color = cluster)) + geom_point() + ArchR::theme_ArchR() +
      scale_color_manual(values = c("A" = "#1B9E77", "B" = "#D95F02", "C" = "#7570B3", "TN" = "red2"))
pdf("output/Plots/S4_PCA_classify.pdf", width = 4, height = 4.5)
p5
dev.off()

se_cdp_er$CA <- hc_class
saveRDS(se_cdp_er, "rds/se_cdp_er.rds")

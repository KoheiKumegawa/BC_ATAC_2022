#------------------------------------------------------------------------------
# 06 - DeconvolutionPeaks
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridisLite)
"%ni%" <- Negate("%in%")
se <- readRDS("rds/se2.rds")
rowRanges(se) <- dropSeqlevels(rowRanges(se), "chrY")

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

#overlaps
subsetByOverlaps(rowRanges(se), cellTypeDAR_total) #19125 overlaps
overlapsCRE <- lapply(cellTypeDAR, function(x) subsetByOverlaps(rowRanges(se), x))
system("mkdir output/Tables/TMEoverlaps output/great_go")
lapply(names(overlapsCRE), function(x){
  gr <- overlapsCRE[[x]]
  d <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr))
  write.table(d, paste0("output/Tables/TMEoverlaps/", x, "_overlaps.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
})

#heatmap
cellTypeDAR <- cellTypeDAR[c(2,3,6,1,5,4)]
idx <- lapply(cellTypeDAR, function(x) findOverlaps(rowRanges(se), x) %>% queryHits %>% unique)
mtx <- assays(se)$normcounts[unlist(idx),]

col_fun1 <- colorRamp2(c(0,100), c("white", "brown"))
ha1 <- HeatmapAnnotation(Pattern = se$pattern, Subtype = se@colData$Subtype, ER = se@colData$ER, PGR = se@colData$PgR, HER2 = se@colData$HER2, Ki67 = se@colData$Ki67, Pathol = se@colData$HistologicalType, Pr_Re = se@colData$Primary_Relapse,
                         annotation_legend_param = c(Pattern = list(title = "Pattern", at = levels(se$pattern)),
                                                     Subtype = list(title = "Subtype", at = levels(se@colData$Subtype)),
                                                     ER = list(title = "ER", at = levels(se@colData$ER)),
                                                     PGR = list(title = "PGR", at = levels(se@colData$PgR)), 
                                                     HER2 = list(title = "PGR", at = levels(se@colData$HER2)),
                                                     Ki67 = list(title = "PAM50", at = c(0,25,50,75,100), col_fun = col_fun1),
                                                     Pathol = list(title = "Pathol", at = levels(se@colData$HistologicalType)),
                                                     Pr_Re = list(title = "Primary_Recurrent", at = levels(se@colData$Primary_Relapse))),
                         col = list(Pattern = brewer.pal(3, "Dark2") %>% `names<-`(., c("A", "B", "C")),
                                    Subtype = c("L" = "darkblue", "TN" = "darkred"),
                                    ER = brewer.pal(length(levels(se@colData$ER)), "Blues") %>% `names<-`(., levels(se@colData$ER)),
                                    PGR = brewer.pal(length(levels(se@colData$PgR)), "Greens") %>% `names<-`(., levels(se@colData$PgR)),
                                    HER2 = brewer.pal(length(levels(se@colData$HER2)), "Oranges") %>% `names<-`(., levels(se@colData$HER2)),
                                    Ki67 = col_fun1,
                                    Pathol = brewer.pal(length(levels(se@colData$HistologicalType)), "Set3") %>% `names<-`(., levels(se@colData$HistologicalType)),
                                    Pr_Re = c("gray80", "black") %>% `names<-`(., c("P", "R")))
)
col_fun2 = colorRamp2(c(0,1,2,3,4,5,6), viridis(7, option = "D"))
fh = function(x) hclust(dist(x), method="ward.D2")
ht1 <- Heatmap(mtx, name = "log2(Normalized counts)", cluster_columns = fh, cluster_rows = F, show_row_names = F, 
               row_split = lapply(idx, length) %>% unlist %>% `names<-`(., names(idx)) %>% rep(names(.), .),
               top_annotation = ha1, col = col_fun2)
p1 <- draw(ht1)
pdf("output/Plots/S2_HeatmapTMEpeaks.pdf", width = 8, height = 8)
p1
dev.off()

#mean signal
df <- lapply(cellTypeDAR, function(x){
  overlaps <- findOverlaps(rowRanges(se), x)
  out <- colMeans(assays(se)$normcounts[queryHits(overlaps),])
  return(out)
})
df <- as.data.frame(do.call(cbind, df))
df <- cbind(df, Subtype = as.character(se@colData$Subtype))
Subtype_col <- c("L" = "darkblue", "TN" = "darkred")
p2 <- lapply(colnames(df)[-7], function(x){
  mtx <- df[,c(x, "Subtype")] %>% `colnames<-`(., c("Mean", "Subtype"))
  p <- ggplot(mtx, aes(x = Subtype, y = Mean, fill = Subtype)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + theme_classic() +
    labs(x = "", y = "Mean ATAC signal", fill = "Subtype", title = x) + 
    scale_fill_manual(values = Subtype_col) +
    theme(axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) +
    ggsignif::geom_signif(comparisons = list(c("L", "TN")), test = "wilcox.test")
  return(p)
})
pdf("output/Plots/S2_BoxplotTMEsignalSubtype.pdf", height = 4, width = 2.5)
p2
dev.off()

#---- Cancer Distal Peaks ----#
se_dst <- se[grep("1kb", mcols(se)$annotation, invert = T), ]
se_cdp <- se_dst[subsetByOverlaps(rowRanges(se_dst), cellTypeDAR_total, invert = T) %>% names,]

saveRDS(se_cdp, "rds/se_cdp.rds")

#---- nucfree peaks ----#
nucfree_se <- readRDS("rds/nucfree_se.rds")

#peak annotation
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peakAnno <- ChIPseeker::annotatePeak(rowRanges(nucfree_se), tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
mcols(nucfree_se) <- mcols(as.GRanges(peakAnno))
nucfree_se_dst <- nucfree_se[grep("1kb", mcols(nucfree_se)$annotation, invert = T), ]
nucfree_se_cdp <- nucfree_se_dst[subsetByOverlaps(rowRanges(nucfree_se_dst), cellTypeDAR_total, invert = T) %>% names,]

saveRDS(nucfree_se_cdp, "rds/nucfree_se_cdp.rds")

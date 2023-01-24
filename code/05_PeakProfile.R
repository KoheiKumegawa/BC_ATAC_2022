#------------------------------------------------------------------------------
# 05 - PeakProfile.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(scales)

se <- readRDS("rds/se.rds")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#peak annotation
peakAnno <- annotatePeak(rowRanges(se), tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
mcols(se) <- mcols(as.GRanges(peakAnno))
pdf("output/Plots/S1_peakAnno.pdf", width = 7, height = 4)
plotAnnoPie(peakAnno) 
dev.off()

#sample correlation in promoter and distal peaks
cor_promoter <- cor(assays(se)$normcounts[grep("1kb", mcols(se)$annotation), ], method = "pearson")
cor_distal   <- cor(assays(se)$normcounts[grep("1kb", mcols(se)$annotation, invert = T), ], method = "pearson")
col_fun1 <- colorRamp2(c(0,100), c("white", "brown"))
ha1 <- HeatmapAnnotation(Subtype = colData(se)$Subtype, ER = colData(se)$ER, PGR = colData(se)$PgR, HER2 = colData(se)$HER2, Ki67 = colData(se)$Ki67, Pathol = colData(se)$HistologicalType, Pr_Re = colData(se)$Primary_Relapse,
                         annotation_legend_param = c(Subtype = list(title = "Subtype", at = levels(colData(se)$Subtype)),
                                                     ER = list(title = "ER", at = levels(colData(se)$ER)),
                                                     PGR = list(title = "PGR", at = levels(colData(se)$PgR)), 
                                                     HER2 = list(title = "PGR", at = levels(colData(se)$HER2)),
                                                     Ki67 = list(title = "PAM50", at = c(0,25,50,75,100), col_fun = col_fun1),
                                                     Pathol = list(title = "Pathol", at = levels(colData(se)$HistologicalType)),
                                                     Pr_Re = list(title = "Primary_Recurrent", at = levels(colData(se)$Primary_Relapse))),
                         col = list(Subtype = c("L" = "darkblue", "LH" = "orangered", "TN" = "darkred"),
                                    ER = brewer.pal(length(levels(colData(se)$ER)), "Blues") %>% `names<-`(., levels(colData(se)$ER)),
                                    PGR = brewer.pal(length(levels(colData(se)$PgR)), "Greens") %>% `names<-`(., levels(colData(se)$PgR)),
                                    HER2 = brewer.pal(length(levels(colData(se)$HER2)), "Oranges") %>% `names<-`(., levels(colData(se)$HER2)),
                                    Ki67 = col_fun1,
                                    Pathol = brewer.pal(length(levels(colData(se)$HistologicalType)), "Set3") %>% `names<-`(., levels(colData(se)$HistologicalType)),
                                    Pr_Re = c("gray80", "black") %>% `names<-`(., c("P", "R")))
)
col_fun2 = colorRamp2(c(0.2, 0.6, 1.0), c("blue", "white", "red"))
fh = function(x) hclust(dist(x), method="ward.D2")
ht1 <- Heatmap(cor_promoter, name = "Pearson's correlation", cluster_columns = fh, cluster_rows = fh, top_annotation = ha1, col = col_fun2)
ht2 <- Heatmap(cor_distal, name = "Pearson's correlation", cluster_columns = fh, cluster_rows = fh, top_annotation = ha1, col = col_fun2)
p1 <- draw(ht1)
p2 <- draw(ht2)
pdf("output/Plots/S1_SampleCorrelation.pdf", height = 8, width = 9)
p1
p2
dev.off()

#overlaps with TCGA-BRCA peaks
tcgapeaks <- data.table::fread("data/TCGA/BRCA_log2norm.txt")[,c(1:4)]
tcgapeaks <- GRanges(seqnames = tcgapeaks$seqnames, IRanges(start = tcgapeaks$start, end = tcgapeaks$end), name = tcgapeaks$name)

peakAnnoTCGA <- annotatePeak(tcgapeaks, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
pdf("output/Plots/S1_peakAnnoTCGA.pdf", width = 7, height = 4)
plotAnnoPie(peakAnnoTCGA) 
dev.off()

overlaps <- findOverlaps(rowRanges(se), tcgapeaks)
n <- queryHits(overlaps) %>% unique %>% length
#subjectHits(overlaps) %>% unique %>% length
df <- data.frame(data = "JFCR-BRCA", Overlaps = n, Unique = c(nrow(se) - n)) %>% reshape2::melt(.)
p5 <- ggplot(df, aes(x = data, y = value, fill = factor(variable, levels = c("Unique","Overlaps")))) + 
  geom_bar(stat = "identity", position = "fill") + theme_classic() + scale_fill_manual(values = c("orange", "darkgray")) +
  labs(x = "", y = "Percent Overlap JFCR-BRCA and TCGA-BRCA", fill = "") + scale_y_continuous(labels = percent, breaks = seq(0, 1, 0.2), expand = c(0,0)) +
  theme(axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90))
pdf("output/Plots/S1_overlapPeaks_TCGA_JFCR.pdf", height = 5, width = 2.5)
p5
dev.off()

peakAnnoOverlaps <- annotatePeak(rowRanges(se)[unique(queryHits(overlaps))], tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
peakAnnoUnique <- annotatePeak(rowRanges(se)[-unique(queryHits(overlaps))], tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
pdf("output/Plots/S1_peakAnnoOverlaps.pdf", width = 7, height = 4)
plotAnnoPie(peakAnnoOverlaps)
plotAnnoPie(peakAnnoUnique)
dev.off()

saveRDS(se, "rds/se2.rds")

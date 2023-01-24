#------------------------------------------------------------------------------
# 09 - CompareCAClusters
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(rtracklayer)
library(dplyr)
library(ggplot2)
source("code/edgeR_PairwiseFunction.R")
se <- readRDS("rds/se2.rds")
se_cdp_er <- readRDS("rds/se_cdp_er.rds")
nucfree_se_cdp <- readRDS("rds/nucfree_se_cdp.rds")

nucfree_se_cdp_er <- nucfree_se_cdp[, colnames(se_cdp_er)]
nucfree_se_cdp_er$CA <- se_cdp_er$CA

#----- TME signal -----#
cellTypeDAR <- import.bed("ref/cellTypeDAR_hg38_scATAC.bed")
n <- names(table(cellTypeDAR$name))[-3]
cellTypeDAR <- lapply(n, function(x){
  res <- cellTypeDAR[cellTypeDAR$name == x]
  res <- subsetByOverlaps(res, cellTypeDAR[cellTypeDAR$name == "Epithelial"], invert = T)
  return(res)
})
names(cellTypeDAR) <- n

df <- lapply(cellTypeDAR, function(x){
  overlaps <- findOverlaps(rowRanges(se), x)
  out <- colMeans(assays(se)$normcounts[queryHits(overlaps),])
  return(out)
})
df <- as.data.frame(do.call(cbind, df))
df <- df[colnames(se_cdp_er), ]
df$cluster <- se_cdp_er$CA

cluster_colors <- RColorBrewer::brewer.pal(3, "Dark2") %>% `names<-`(., c("A", "B", "C"))
p1 <- lapply(colnames(df)[-7], function(x){
  mtx <- df[,c(x, "cluster")] %>% `colnames<-`(., c("Mean", "Cluster"))
  p <- ggplot(mtx, aes(x = Cluster, y = Mean, fill = Cluster)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + theme_classic() +
    labs(x = "", y = "Mean ATAC signal", fill = "CA cluster", title = x) + 
    scale_fill_manual(values = cluster_colors) +
    theme(axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) +
    ggsignif::geom_signif(comparisons = list(c("A", "B"), c("B", "C"), c("A", "C")), test = "wilcox.test")
  return(p)
})
pdf("output/Plots/S5_BoxplotTMEsignalCluster.pdf", height = 4, width = 3.5)
p1
dev.off()

#----- CA-A vs CA-B -----#
#total peaks
diffPeaksAB <- edgeR_pairwise(se_cdp_er, compareCol = "CA", topGroup = "A", bottomGroup = "B")
Apeaks_AB <- rownames(diffPeaksAB)[which(assay(diffPeaksAB)[,"log2FoldChange"] > 1 & assay(diffPeaksAB)[,"FDR"] < 0.01)]
Bpeaks_AB <- rownames(diffPeaksAB)[which(assay(diffPeaksAB)[,"log2FoldChange"] < -1 & assay(diffPeaksAB)[,"FDR"] < 0.01)]
#nucfree peaks
nucfree_diffPeaksAB <- edgeR_pairwise(nucfree_se_cdp_er, compareCol = "CA", topGroup = "A", bottomGroup = "B")
nucfree_Apeaks_AB <- rownames(nucfree_diffPeaksAB)[which(assay(nucfree_diffPeaksAB)[,"log2FoldChange"] > 1 & assay(nucfree_diffPeaksAB)[,"FDR"] < 0.01)]
nucfree_Bpeaks_AB <- rownames(nucfree_diffPeaksAB)[which(assay(nucfree_diffPeaksAB)[,"log2FoldChange"] < -1 & assay(nucfree_diffPeaksAB)[,"FDR"] < 0.01)]

#export
gr <- GRangesList(Apeaks_AB = rowRanges(se_cdp_er)[Apeaks_AB],
                  Bpeaks_AB = rowRanges(se_cdp_er)[Bpeaks_AB],
                  nucfree_Apeaks_AB = rowRanges(nucfree_se_cdp_er)[nucfree_Apeaks_AB],
                  nucfree_Bpeaks_AB = rowRanges(nucfree_se_cdp_er)[nucfree_Bpeaks_AB])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/AvsB/", x, "_output.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#motif visualization
dirlist <- list.dirs("output/homer_motif/AvsB/", recursive = F, full.names = F)
motif_DF <- lapply(dirlist, function(x){
  out <- data.table::fread(paste0("output/homer_motif/AvsB/", x, "/knownResults.txt"), header = F, skip = 1)[c(1:10), c(1,3)] %>% data.frame %>% `colnames<-`(., c("Motif", "P"))
  out$mlog10P <- as.numeric(gsub("1e-", "", out$P))
  return(out)
}) %>% `names<-`(., dirlist)
p2 <- lapply(names(motif_DF), function(i){
  df <- motif_DF[[i]]
  p <- ggplot(df, aes(x = mlog10P, y = reorder(Motif, mlog10P))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
    labs(x = "Motif", y = "-log10(P-value)") + ggtitle(i)
  return(p)
})
pdf("output/Plots/S5_HomerMotif_AvsB.pdf", height = 3, width = 10)
p2
dev.off()

#GREAT visualization
GREAT_A_AB <- data.table::fread("output/great_go/AvsB/Apeaks_AB.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_B_AB <- data.table::fread("output/great_go/AvsB/Bpeaks_AB.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_A_AB$mlog10FDR <- -log10(GREAT_A_AB$FDR)
GREAT_B_AB$mlog10FDR <- -log10(GREAT_B_AB$FDR)
p3 <- ggplot(GREAT_A_AB, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 
p4 <- ggplot(GREAT_B_AB, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 
pdf("output/Plots/S5_GREAT_AvsB.pdf", height = 3, width = 10)
p3
p4
dev.off()

#volcano plot
cluster_colors2 <- c(RColorBrewer::brewer.pal(3, "Dark2"), "darkgray") %>% `names<-`(., c("A", "B", "C", "N"))

df <- assay(diffPeaksAB)[, c("log2FoldChange", "FDR")] %>% data.frame
df$mlog10FDR <- -log10(df$FDR)
df$signif <- "N"
df$signif[which(df$log2FoldChange > 1 & df$FDR < 0.01)] <- "A"
df$signif[which(df$log2FoldChange < -1 & df$FDR < 0.01)] <- "B"
df$signif <- factor(df$signif, levels = c("N", "A", "B"))
df <- df[order(df$signif),]
p5 <- ggplot(df, aes(x = log2FoldChange, y = mlog10FDR, color = signif)) + geom_point(size = 0.8) + ArchR::theme_ArchR() +
  geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + scale_color_manual(values = cluster_colors2) + 
  ggtitle("Total peaks") + labs(x = "log2FC", y = "-log10(FDR)")

df <- assay(nucfree_diffPeaksAB)[, c("log2FoldChange", "FDR")] %>% data.frame
df$mlog10FDR <- -log10(df$FDR)
df$signif <- "N"
df$signif[which(df$log2FoldChange > 1 & df$FDR < 0.01)] <- "A"
df$signif[which(df$log2FoldChange < -1 & df$FDR < 0.01)] <- "B"
df$signif <- factor(df$signif, levels = c("N", "A", "B"))
df <- df[order(df$signif),]
p6 <- ggplot(df, aes(x = log2FoldChange, y = mlog10FDR, color = signif)) + geom_point(size = 0.8) + ArchR::theme_ArchR() +
  geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + scale_color_manual(values = cluster_colors2) + 
  ggtitle("Nucleosome-free peaks") + labs(x = "log2FC", y = "-log10(FDR)")
pdf("output/Plots/S5_Volcano_AvsB.pdf", width = 5, height = 5.5)
p5
p6
dev.off()

nearby_Apeaks_AB <- rowRanges(se_cdp_er)[Apeaks_AB]$SYMBOL %>% table %>% sort(., decreasing = T)
nearby_Bpeaks_AB <- rowRanges(se_cdp_er)[Bpeaks_AB]$SYMBOL %>% table %>% sort(., decreasing = T)
write.csv(data.frame(gene = names(nearby_Apeaks_AB), number = as.numeric(nearby_Apeaks_AB)), "output/nearbygenes/AvsB/nearby_Apeaks.csv", row.names = F, quote = F)
write.csv(data.frame(gene = names(nearby_Bpeaks_AB), number = as.numeric(nearby_Bpeaks_AB)), "output/nearbygenes/AvsB/nearby_Bpeaks.csv", row.names = F, quote = F)

df <- data.frame(gene = names(nearby_Apeaks_AB), number = as.numeric(nearby_Apeaks_AB), rank = 1:length(nearby_Apeaks_AB))
p6.2 <- ggplot(df, aes(x = rank, y = number)) + geom_point(size = 0.5) + ArchR::theme_ArchR() + ggtitle("CA-A specific") +
        scale_y_continuous(breaks = seq(0,10,1))
df <- data.frame(gene = names(nearby_Bpeaks_AB), number = as.numeric(nearby_Bpeaks_AB), rank = 1:length(nearby_Bpeaks_AB))
p6.3 <- ggplot(df, aes(x = rank, y = number)) + geom_point(size = 0.5) + ArchR::theme_ArchR() + ggtitle("CA-B specific") +
        scale_y_continuous(breaks = seq(0,9,1))
pdf("output/Plots/S5_RankPlot_nearbygenes_AvsB.pdf", width = 3, height = 3)
p6.2
p6.3
dev.off()

#----- CA-B vs CA-C -----#
#total peaks
diffPeaksBC <- edgeR_pairwise(se_cdp_er, compareCol = "CA", topGroup = "B", bottomGroup = "C")
Bpeaks_BC <- rownames(diffPeaksBC)[which(assay(diffPeaksBC)[,"log2FoldChange"] > 1 & assay(diffPeaksBC)[,"FDR"] < 0.01)]
Cpeaks_BC <- rownames(diffPeaksBC)[which(assay(diffPeaksBC)[,"log2FoldChange"] < -1 & assay(diffPeaksBC)[,"FDR"] < 0.01)]
#nucfree peaks
nucfree_diffPeaksBC <- edgeR_pairwise(nucfree_se_cdp_er, compareCol = "CA", topGroup = "B", bottomGroup = "C")
nucfree_Bpeaks_BC <- rownames(nucfree_diffPeaksBC)[which(assay(nucfree_diffPeaksBC)[,"log2FoldChange"] > 1 & assay(nucfree_diffPeaksBC)[,"FDR"] < 0.01)]
nucfree_Cpeaks_BC <- rownames(nucfree_diffPeaksBC)[which(assay(nucfree_diffPeaksBC)[,"log2FoldChange"] < -1 & assay(nucfree_diffPeaksBC)[,"FDR"] < 0.01)]

#export
gr <- GRangesList(Bpeaks_BC = rowRanges(se_cdp_er)[Bpeaks_BC],
                  Cpeaks_BC = rowRanges(se_cdp_er)[Cpeaks_BC],
                  nucfree_Bpeaks_BC = rowRanges(nucfree_se_cdp_er)[nucfree_Bpeaks_BC],
                  nucfree_Cpeaks_BC = rowRanges(nucfree_se_cdp_er)[nucfree_Cpeaks_BC])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/BvsC/", x, "_output.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#motif visualization
dirlist <- list.dirs("output/homer_motif/BvsC/", recursive = F, full.names = F)
motif_DF <- lapply(dirlist, function(x){
  out <- data.table::fread(paste0("output/homer_motif/BvsC/", x, "/knownResults.txt"), header = F, skip = 1)[c(1:10), c(1,3)] %>% data.frame %>% `colnames<-`(., c("Motif", "P"))
  out$mlog10P <- as.numeric(gsub("1e-", "", out$P))
  return(out)
}) %>% `names<-`(., dirlist)
p7 <- lapply(names(motif_DF), function(i){
  df <- motif_DF[[i]]
  p <- ggplot(df, aes(x = mlog10P, y = reorder(Motif, mlog10P))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
    labs(x = "Motif", y = "-log10(P-value)") + ggtitle(i)
  return(p)
})
pdf("output/Plots/S5_HomerMotif_BvsC.pdf", height = 3, width = 10)
p7
dev.off()

#GREAT visualization
GREAT_C_BC <- data.table::fread("output/great_go/BvsC/Cpeaks_BC.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_C_BC$mlog10FDR <- -log10(GREAT_C_BC$FDR)
p8 <- ggplot(GREAT_C_BC, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 
pdf("output/Plots/S5_GREAT_BvsC.pdf", height = 3, width = 10)
p8
dev.off()

#volcano plot
df <- assay(diffPeaksBC)[, c("log2FoldChange", "FDR")] %>% data.frame
df$mlog10FDR <- -log10(df$FDR)
df$signif <- "N"
df$signif[which(df$log2FoldChange > 1 & df$FDR < 0.01)] <- "B"
df$signif[which(df$log2FoldChange < -1 & df$FDR < 0.01)] <- "C"
df$signif <- factor(df$signif, levels = c("N", "B", "C"))
df <- df[order(df$signif),]
p9 <- ggplot(df, aes(x = log2FoldChange, y = mlog10FDR, color = signif)) + geom_point(size = 0.8) + ArchR::theme_ArchR() +
  geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + scale_color_manual(values = cluster_colors2) + 
  ggtitle("Total peaks") + labs(x = "log2FC", y = "-log10(FDR)")

df <- assay(nucfree_diffPeaksBC)[, c("log2FoldChange", "FDR")] %>% data.frame
df$mlog10FDR <- -log10(df$FDR)
df$signif <- "N"
df$signif[which(df$log2FoldChange > 1 & df$FDR < 0.01)] <- "B"
df$signif[which(df$log2FoldChange < -1 & df$FDR < 0.01)] <- "C"
df$signif <- factor(df$signif, levels = c("N", "B", "C"))
df <- df[order(df$signif),]
p10 <- ggplot(df, aes(x = log2FoldChange, y = mlog10FDR, color = signif)) + geom_point(size = 0.8) + ArchR::theme_ArchR() +
  geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + scale_color_manual(values = cluster_colors2) + 
  ggtitle("Nucleosome-free peaks") + labs(x = "log2FC", y = "-log10(FDR)")
pdf("output/Plots/S5_Volcano_BvsC.pdf", width = 5, height = 5.5)
p9
p10
dev.off()

#----- CA-A vs CA-C -----#
#total peaks
diffPeaksAC <- edgeR_pairwise(se_cdp_er, compareCol = "CA", topGroup = "A", bottomGroup = "C")
Apeaks_AC <- rownames(diffPeaksAC)[which(assay(diffPeaksAC)[,"log2FoldChange"] > 1 & assay(diffPeaksAC)[,"FDR"] < 0.01)]
Cpeaks_AC <- rownames(diffPeaksAC)[which(assay(diffPeaksAC)[,"log2FoldChange"] < -1 & assay(diffPeaksAC)[,"FDR"] < 0.01)]
#nucfree peaks
nucfree_diffPeaksAC <- edgeR_pairwise(nucfree_se_cdp_er, compareCol = "CA", topGroup = "A", bottomGroup = "C")
nucfree_Apeaks_AC <- rownames(nucfree_diffPeaksAC)[which(assay(nucfree_diffPeaksAC)[,"log2FoldChange"] > 1 & assay(nucfree_diffPeaksAC)[,"FDR"] < 0.01)]
nucfree_Cpeaks_AC <- rownames(nucfree_diffPeaksAC)[which(assay(nucfree_diffPeaksAC)[,"log2FoldChange"] < -1 & assay(nucfree_diffPeaksAC)[,"FDR"] < 0.01)]

#export
gr <- GRangesList(Apeaks_AC = rowRanges(se_cdp_er)[Apeaks_AC],
                  Cpeaks_AC = rowRanges(se_cdp_er)[Cpeaks_AC],
                  nucfree_Apeaks_AC = rowRanges(nucfree_se_cdp_er)[nucfree_Apeaks_AC],
                  nucfree_Cpeaks_AC = rowRanges(nucfree_se_cdp_er)[nucfree_Cpeaks_AC])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/AvsC//", x, "_output.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#motif visualization
dirlist <- list.dirs("output/homer_motif/AvsC/", recursive = F, full.names = F)
motif_DF <- lapply(dirlist, function(x){
  out <- data.table::fread(paste0("output/homer_motif/AvsC/", x, "/knownResults.txt"), header = F, skip = 1)[c(1:10), c(1,3)] %>% data.frame %>% `colnames<-`(., c("Motif", "P"))
  out$mlog10P <- as.numeric(gsub("1e-", "", out$P))
  return(out)
}) %>% `names<-`(., dirlist)
p11 <- lapply(names(motif_DF), function(i){
  df <- motif_DF[[i]]
  p <- ggplot(df, aes(x = mlog10P, y = reorder(Motif, mlog10P))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
    labs(x = "Motif", y = "-log10(P-value)") + ggtitle(i)
  return(p)
})
pdf("output/Plots/S5_HomerMotif_AvsC.pdf", height = 3, width = 10)
p11
dev.off()

#GREAT visualization
GREAT_A_AC <- data.table::fread("output/great_go/AvsC/Apeaks_AC.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_C_AC <- data.table::fread("output/great_go/AvsC/Cpeaks_AC.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_A_AC$mlog10FDR <- -log10(GREAT_A_AC$FDR)
GREAT_C_AC$mlog10FDR <- -log10(GREAT_C_AC$FDR)

p12 <- ggplot(GREAT_A_AC, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 
p13 <- ggplot(GREAT_C_AC, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 
pdf("output/Plots/S5_GREAT_AvsC.pdf", height = 3, width = 10)
p12
p13
dev.off()

#volcano plot
df <- assay(diffPeaksAC)[, c("log2FoldChange", "FDR")] %>% data.frame
df$mlog10FDR <- -log10(df$FDR)
df$signif <- "N"
df$signif[which(df$log2FoldChange > 1 & df$FDR < 0.01)] <- "A"
df$signif[which(df$log2FoldChange < -1 & df$FDR < 0.01)] <- "C"
df$signif <- factor(df$signif, levels = c("N", "A", "C"))
df <- df[order(df$signif),]
p14 <- ggplot(df, aes(x = log2FoldChange, y = mlog10FDR, color = signif)) + geom_point(size = 0.8) + ArchR::theme_ArchR() +
  geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + scale_color_manual(values = cluster_colors2) + 
  ggtitle("Total peaks") + labs(x = "log2FC", y = "-log10(FDR)")

df <- assay(nucfree_diffPeaksAC)[, c("log2FoldChange", "FDR")] %>% data.frame
df$mlog10FDR <- -log10(df$FDR)
df$signif <- "N"
df$signif[which(df$log2FoldChange > 1 & df$FDR < 0.01)] <- "A"
df$signif[which(df$log2FoldChange < -1 & df$FDR < 0.01)] <- "C"
df$signif <- factor(df$signif, levels = c("N", "A", "C"))
df <- df[order(df$signif),]
p15 <- ggplot(df, aes(x = log2FoldChange, y = mlog10FDR, color = signif)) + geom_point(size = 0.8) + ArchR::theme_ArchR() +
  geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + scale_color_manual(values = cluster_colors2) + 
  ggtitle("Nucleosome-free peaks") + labs(x = "log2FC", y = "-log10(FDR)")
pdf("output/Plots/S5_Volcano_AvsC.pdf", width = 5, height = 5.5)
p14
p15
dev.off()

#----- chromVAR motif score -----#
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
#add GC bias
se_cdp_er <- addGCBias(se_cdp_er, genome = BSgenome.Hsapiens.UCSC.hg38)
#homer motif
motifs <- c(homer_pwms)
motif_ix <- matchMotifs(motifs, se_cdp_er, genome = BSgenome.Hsapiens.UCSC.hg38)
motif_ps <- matchMotifs(motifs, se_cdp_er, genome = BSgenome.Hsapiens.UCSC.hg38, out = "positions")
#calculate deviations
bg <- getBackgroundPeaks(object = se_cdp_er)
dev <- computeDeviations(object = se_cdp_er, annotations = motif_ix, background_peaks = bg)

#nucleosome-free peaks
nucfree_se_cdp_er2 <- nucfree_se_cdp_er[which(rowSums(assays(nucfree_se_cdp_er)$counts) != 0), ]
nucfree_se_cdp_er2 <- addGCBias(nucfree_se_cdp_er2, genome = BSgenome.Hsapiens.UCSC.hg38)
nucfree_motif_ix <- matchMotifs(motifs, nucfree_se_cdp_er2, genome = BSgenome.Hsapiens.UCSC.hg38)
nucfree_motif_ps <- matchMotifs(motifs, nucfree_se_cdp_er2, genome = BSgenome.Hsapiens.UCSC.hg38, out = "positions")
nucfree_bg <- getBackgroundPeaks(object = nucfree_se_cdp_er2)
nucfree_dev <- computeDeviations(object = nucfree_se_cdp_er2, annotations = nucfree_motif_ix, background_peaks = nucfree_bg)

interest_motifs <- c("AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
                     "CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer",
                     "ERE(NR),IR3/MCF7-ERa-ChIP-Seq(Unpublished)/Homer",
                     "FOXA1(Forkhead)/MCF7-FOXA1-ChIP-Seq(GSE26831)/Homer",
                     "Sox3(HMG)/NPC-Sox3-ChIP-Seq(GSE33059)/Homer",
                     "PGR(NR)/EndoStromal-PGR-ChIP-Seq(GSE69539)/Homer")
p16 <- lapply(interest_motifs, function(x){
  d <- dev[x,]
  df <- data.frame(t(d@assays@data$z)) %>% `colnames<-`(., "MotifScore")
  df$CA <- d@colData$CA
  p <- ggplot(df, aes(x = CA, y = MotifScore, fill = CA)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    theme_classic() + theme(legend.position = "none") + 
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2") %>% `names<-`(., c("A", "B", "C"))) + 
    ggsignif::geom_signif(comparisons = list(c("A", "B"), c("A", "C"), c("B", "C")), test = "t.test") + 
    ggtitle(as.character(d@elementMetadata$name))
  return(p)
})
p17 <- lapply(interest_motifs, function(x){
  d <- nucfree_dev[x,]
  df <- data.frame(t(d@assays@data$z)) %>% `colnames<-`(., "MotifScore")
  df$CA <- d@colData$CA
  p <- ggplot(df, aes(x = CA, y = MotifScore, fill = CA)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    theme_classic() + theme(legend.position = "none") + 
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2") %>% `names<-`(., c("A", "B", "C"))) + 
    ggsignif::geom_signif(comparisons = list(c("A", "B"), c("A", "C"), c("B", "C")), test = "t.test") + 
    ggtitle(paste0(as.character(d@elementMetadata$name), ": Nucleosome-free peaks"))
  return(p)
})
pdf("output/Plots//S5_ChromVARscore_keyTF.pdf", width = 3, height = 6)
p16
p17
dev.off()

saveRDS(dev, "rds/dev.rds")
saveRDS(nucfree_dev, "rds/nucfree_dev.rds")
saveRDS(motif_ps, "rds/motif_ps.rds")
saveRDS(nucfree_motif_ps, "rds/nucfree_motif_ps.rds")

#----- Compare CA-A vs CA-B, only IDC -----#
se_cdp_er_IDC <- se_cdp_er[, se_cdp_er$HistologicalType == "IDC"]
nucfree_se_cdp_er_IDC <- nucfree_se_cdp_er[, nucfree_se_cdp_er$HistologicalType == "IDC"]

table(se_cdp_er_IDC$CA, se_cdp_er_IDC$HistologicalType)

#total peaks
IDC_diffPeaksAB <- edgeR_pairwise(se_cdp_er_IDC, compareCol = "CA", topGroup = "A", bottomGroup = "B")
IDC_Apeaks_AB <- rownames(IDC_diffPeaksAB)[which(assay(IDC_diffPeaksAB)[,"log2FoldChange"] > 1 & assay(IDC_diffPeaksAB)[,"FDR"] < 0.05)]
IDC_Bpeaks_AB <- rownames(IDC_diffPeaksAB)[which(assay(IDC_diffPeaksAB)[,"log2FoldChange"] < -1 & assay(IDC_diffPeaksAB)[,"FDR"] < 0.05)]
#nucfree peaks
IDC_nucfree_diffPeaksAB <- edgeR_pairwise(nucfree_se_cdp_er_IDC, compareCol = "CA", topGroup = "A", bottomGroup = "B")
IDC_nucfree_Apeaks_AB <- rownames(IDC_nucfree_diffPeaksAB)[which(assay(IDC_nucfree_diffPeaksAB)[,"log2FoldChange"] > 1 & assay(IDC_nucfree_diffPeaksAB)[,"FDR"] < 0.05)]
IDC_nucfree_Bpeaks_AB <- rownames(IDC_nucfree_diffPeaksAB)[which(assay(IDC_nucfree_diffPeaksAB)[,"log2FoldChange"] < -1 & assay(IDC_nucfree_diffPeaksAB)[,"FDR"] < 0.05)]

#export
gr <- GRangesList(IDC_Apeaks_AB = rowRanges(se_cdp_er_IDC)[IDC_Apeaks_AB],
                  IDC_Bpeaks_AB = rowRanges(se_cdp_er_IDC)[IDC_Bpeaks_AB],
                  IDC_nucfree_Apeaks_AB = rowRanges(nucfree_se_cdp_er_IDC)[IDC_nucfree_Apeaks_AB],
                  IDC_nucfree_Bpeaks_AB = rowRanges(nucfree_se_cdp_er_IDC)[IDC_nucfree_Bpeaks_AB])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/AvsB_IDC/", x, "_output.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#motif visualization
dirlist <- list.dirs("output/homer_motif/AvsB_IDC/", recursive = F, full.names = F)
motif_DF <- lapply(dirlist, function(x){
  out <- data.table::fread(paste0("output/homer_motif/AvsB_IDC/", x, "/knownResults.txt"), header = F, skip = 1)[c(1:10), c(1,3)] %>% data.frame %>% `colnames<-`(., c("Motif", "P"))
  out$mlog10P <- as.numeric(gsub("1e-", "", out$P))
  return(out)
}) %>% `names<-`(., dirlist)
p18 <- lapply(names(motif_DF), function(i){
  df <- motif_DF[[i]]
  p <- ggplot(df, aes(x = mlog10P, y = reorder(Motif, mlog10P))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
    labs(x = "Motif", y = "-log10(P-value)") + ggtitle(i)
  return(p)
})
pdf("output/Plots/S5_HomerMotif_AvsB_IDC.pdf", height = 3, width = 10)
p18
dev.off()


#----- Compare CA by clinical features -----#
#hscore
hscore <- read.csv("ref/20220606_Hscore.csv") %>% `colnames<-`(., c("Patient", "Cluster", "ER", "PGR"))
df <- as.data.frame(hscore)
rownames(df) <- paste0("P", df$Patient)
df <- df[colnames(se_cdp_er),]
df$Cluster <- se_cdp_er$CA

p19 <- ggplot(df, aes(x = Cluster, y = ER, fill = Cluster)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  theme_classic() + theme(legend.position = "none") + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2") %>% `names<-`(., c("A", "B", "C"))) + 
  ggsignif::geom_signif(comparisons = list(c("A", "B"), c("A", "C"), c("B", "C")), test = "t.test") + 
  ggtitle("ER H-score")
p20 <- ggplot(df, aes(x = Cluster, y = PGR, fill = Cluster)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  theme_classic() + theme(legend.position = "none") + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2") %>% `names<-`(., c("A", "B", "C"))) + 
  ggsignif::geom_signif(comparisons = list(c("A", "B"), c("A", "C"), c("B", "C")), test = "t.test") + 
  ggtitle("PGR H-score")
pdf("output/Plots//S5_Hscore_ER_PGR.pdf", width = 3, height = 6)
p19
p20
dev.off()

#other features
df <- colData(se_cdp_er)[, c("Age", "ER", "PgR", "HER2", "Ki67")] %>% data.frame
df <- as.data.frame(apply(apply(df, 2, as.character), 2, as.numeric))
df$Cluster <- se_cdp_er$CA
p21 <- lapply(colnames(df), function(x){
  df2 <- df[,c(x, "Cluster")]
  colnames(df2) <- c("value", "Cluster")
  p <- ggplot(df2, aes(x = Cluster, y = value, fill = Cluster)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    theme_classic() + theme(legend.position = "none") + 
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2") %>% `names<-`(., c("A", "B", "C"))) + 
    ggsignif::geom_signif(comparisons = list(c("A", "B"), c("A", "C"), c("B", "C")), test = "t.test") + 
    ggtitle(x)
  return(p)
})
pdf("output/Plots//S5_clinicalInfo.pdf", width = 3, height = 6)
p21
dev.off()

saveRDS(nucfree_se_cdp_er, "rds/nucfree_se_cdp_er.rds")

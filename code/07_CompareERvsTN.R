#------------------------------------------------------------------------------
# 07 - CompareERvsTN
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(rtracklayer)
library(dplyr)
library(ggplot2)
source("code/edgeR_PairwiseFunction.R")
se_cdp <- readRDS("rds/se_cdp.rds")
nucfree_se_cdp <- readRDS("rds/nucfree_se_cdp.rds")

#total peaks
diffPeaks <- edgeR_pairwise(se_cdp, compareCol = "Subtype", topGroup = "L", bottomGroup = "TN")
ERpeaks <- rownames(diffPeaks)[which(assay(diffPeaks)[,"log2FoldChange"] > 1 & assay(diffPeaks)[,"FDR"] < 0.01)]
TNpeaks <- rownames(diffPeaks)[which(assay(diffPeaks)[,"log2FoldChange"] < -1 & assay(diffPeaks)[,"FDR"] < 0.01)]

#nucfree peaks
nucfree_diffPeaks <- edgeR_pairwise(nucfree_se_cdp, compareCol = "Subtype", topGroup = "L", bottomGroup = "TN")
nucfree_ERpeaks <- rownames(nucfree_diffPeaks)[which(assay(nucfree_diffPeaks)[,"log2FoldChange"] > 1 & assay(nucfree_diffPeaks)[,"FDR"] < 0.01)]
nucfree_TNpeaks <- rownames(nucfree_diffPeaks)[which(assay(nucfree_diffPeaks)[,"log2FoldChange"] < -1 & assay(nucfree_diffPeaks)[,"FDR"] < 0.01)]

#export
gr <- GRangesList(ERpeaks = rowRanges(se_cdp)[ERpeaks],
                  TNpeaks = rowRanges(se_cdp)[TNpeaks],
                  nucfree_ERpeaks = rowRanges(nucfree_se_cdp)[nucfree_ERpeaks],
                  nucfree_TNpeaks = rowRanges(nucfree_se_cdp)[nucfree_TNpeaks])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/ERvsTN/", x, "_output.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

## GREAT analysis & HOMER motif analysis, shell script
#rowRanges(se_cdp)[ERpeaks]$SYMBOL %>% table %>% sort(., decreasing = T)
#rowRanges(se_cdp)[TNpeaks]$SYMBOL %>% table %>% sort(., decreasing = T)

#GREAT visualization
GREAT_ER <- data.table::fread("output/great_go/ERvsTN/ERpeaks.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_TN <- data.table::fread("output/great_go/ERvsTN/TMpeaks.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_ER$mlog10FDR <- -log10(GREAT_ER$FDR)
GREAT_TN$mlog10FDR <- -log10(GREAT_TN$FDR)

p1 <- ggplot(GREAT_ER, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 
p2 <- ggplot(GREAT_TN, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 

#motif visualization
dirlist <- list.dirs("output/homer_motif/ERvsTN/", recursive = F, full.names = F)
motif_DF <- lapply(dirlist, function(x){
  out <- data.table::fread(paste0("output/homer_motif/ERvsTN/", x, "/knownResults.txt"), header = F, skip = 1)[c(1:10), c(1,3)] %>% data.frame %>% `colnames<-`(., c("Motif", "P"))
  out$mlog10P <- as.numeric(gsub("1e-", "", out$P))
  return(out)
}) %>% `names<-`(., dirlist)

p3 <- lapply(names(motif_DF), function(i){
  df <- motif_DF[[i]]
  p <- ggplot(df, aes(x = mlog10P, y = reorder(Motif, mlog10P))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
    labs(x = "Motif", y = "-log10(P-value)") + ggtitle(i)
  return(p)
})

#volcano
Subtype_col <- c("ER" = "darkblue", "TN" = "darkred", "N" = "darkgray")

df <- assay(diffPeaks)[, c("log2FoldChange", "FDR")] %>% data.frame
df$mlog10FDR <- -log10(df$FDR)
df$signif <- "N"
df$signif[which(df$log2FoldChange > 1 & df$FDR < 0.01)] <- "ER"
df$signif[which(df$log2FoldChange < -1 & df$FDR < 0.01)] <- "TN"
df$signif <- factor(df$signif, levels = c("N", "ER", "TN"))
df <- df[order(df$signif),]
p4 <- ggplot(df, aes(x = log2FoldChange, y = mlog10FDR, color = signif)) + geom_point() + ArchR::theme_ArchR() +
      geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + scale_color_manual(values = Subtype_col) + 
      ggtitle("Total peaks") + labs(x = "log2FC", "y = -log10(FDR)")

df <- assay(nucfree_diffPeaks)[, c("log2FoldChange", "FDR")] %>% data.frame
df$mlog10FDR <- -log10(df$FDR)
df$signif <- "N"
df$signif[which(df$log2FoldChange > 1 & df$FDR < 0.01)] <- "ER"
df$signif[which(df$log2FoldChange < -1 & df$FDR < 0.01)] <- "TN"
df$signif <- factor(df$signif, levels = c("N", "ER", "TN"))
df <- df[order(df$signif),]
p5 <- ggplot(df, aes(x = log2FoldChange, y = mlog10FDR, color = signif)) + geom_point() + ArchR::theme_ArchR() +
      geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + scale_color_manual(values = Subtype_col) + 
      ggtitle("Nucleosome-free peaks") + labs(x = "log2FC", y = "-log10(FDR)")

#output
pdf("output/Plots/S3_GREAT_ER.pdf", width = 10, height = 3)
p1
dev.off()
pdf("output/Plots/S3_GREAT_TN.pdf", width = 10, height = 3)
p2
dev.off()
pdf("output/Plots/S3_HomerMotif_ERvsTN.pdf", width = 10, height = 3)
p3
dev.off()
pdf("output/Plots/S3_Volcano_ERvsTN.pdf", width = 5, height = 5.5)
p4
dev.off()
pdf("output/Plots/S3_Volcano_ERvsTN_nucfree.pdf", width = 5, height = 5.5)
p5
dev.off()

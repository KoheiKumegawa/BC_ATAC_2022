#------------------------------------------------------------------------------
# 12 - CompareTCGACAClusters
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
source("code/edgeR_PairwiseFunction.R")
TCGA_BRCA_ATAC_ER_se_cdp <- readRDS("rds/TCGA_BRCA_ATAC_ER_se_cdp.rds")

#----- chromVAR motif score -----#
#add GC bias
TCGA_BRCA_ATAC_ER_se_cdp <- addGCBias(TCGA_BRCA_ATAC_ER_se_cdp, genome = BSgenome.Hsapiens.UCSC.hg38)
#homer motif
motifs <- c(homer_pwms)
motif_ix <- matchMotifs(motifs, TCGA_BRCA_ATAC_ER_se_cdp, genome = BSgenome.Hsapiens.UCSC.hg38)
motif_ps <- matchMotifs(motifs, TCGA_BRCA_ATAC_ER_se_cdp, genome = BSgenome.Hsapiens.UCSC.hg38, out = "positions")
#calculate deviations
bg <- getBackgroundPeaks(object = TCGA_BRCA_ATAC_ER_se_cdp)
dev <- computeDeviations(object = TCGA_BRCA_ATAC_ER_se_cdp, annotations = motif_ix, background_peaks = bg)

interest_motifs <- c("AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
                     "CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer",
                     "ERE(NR),IR3/MCF7-ERa-ChIP-Seq(Unpublished)/Homer",
                     "FOXA1(Forkhead)/MCF7-FOXA1-ChIP-Seq(GSE26831)/Homer",
                     "Sox3(HMG)/NPC-Sox3-ChIP-Seq(GSE33059)/Homer",
                     "PGR(NR)/EndoStromal-PGR-ChIP-Seq(GSE69539)/Homer")
p1 <- lapply(interest_motifs, function(x){
  d <- dev[x,]
  df <- data.frame(t(d@assays@data$z)) %>% `colnames<-`(., "MotifScore")
  df$CA <- d@colData$CA
  p <- ggplot(df, aes(x = CA, y = MotifScore, fill = CA)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    theme_classic() + theme(legend.position = "none") + 
    scale_fill_manual(values = c("A" = "#1B9E77", "B" = "#D95F02", "ILC" = "#E7298A")) + 
    ggsignif::geom_signif(comparisons = list(c("A", "B"), c("A", "ILC"), c("B", "ILC")), test = "t.test") + 
    ggtitle(as.character(d@elementMetadata$name))
  return(p)
})
pdf("output/Plots/S7_TCGA_ChromVARscore_keyTF.pdf", width = 3, height = 6)
p1
dev.off()

saveRDS(dev, "rds/TCGA_dev.rds")
saveRDS(motif_ps, "rds/TCGA_motif_ps.rds")

#----- CA-A vs CA-B -----#
#total peaks
diffPeaksAB <- edgeR_pairwise(TCGA_BRCA_ATAC_ER_se_cdp, compareCol = "CA", topGroup = "A", bottomGroup = "B")
Apeaks_AB <- rownames(diffPeaksAB)[which(assay(diffPeaksAB)[,"log2FoldChange"] > 1 & assay(diffPeaksAB)[,"FDR"] < 0.01)]
Bpeaks_AB <- rownames(diffPeaksAB)[which(assay(diffPeaksAB)[,"log2FoldChange"] < -1 & assay(diffPeaksAB)[,"FDR"] < 0.01)]

#export
gr <- GRangesList(Apeaks_AB = rowRanges(TCGA_BRCA_ATAC_ER_se_cdp)[Apeaks_AB], Bpeaks_AB = rowRanges(TCGA_BRCA_ATAC_ER_se_cdp)[Bpeaks_AB])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/TCGA_AvsB/", x, "_output.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#motif visualization
dirlist <- list.dirs("output/homer_motif/TCGA_AvsB/", recursive = F, full.names = F)
motif_DF <- lapply(dirlist, function(x){
  out <- data.table::fread(paste0("output/homer_motif/TCGA_AvsB/", x, "/knownResults.txt"), header = F, skip = 1)[c(1:10), c(1,3)] %>% data.frame %>% `colnames<-`(., c("Motif", "P"))
  out$mlog10P <- as.numeric(gsub("1e-", "", out$P))
  return(out)
}) %>% `names<-`(., dirlist)
p2 <- lapply(names(motif_DF), function(i){
  df <- motif_DF[[i]]
  p <- ggplot(df, aes(x = mlog10P, y = reorder(Motif, mlog10P))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
    labs(x = "Motif", y = "-log10(P-value)") + ggtitle(i)
  return(p)
})
pdf("output/Plots/S7_TCGA_HomerMotif_AvsB.pdf", height = 3, width = 10)
p2
dev.off()

#GREAT visualization
GREAT_A_AB <- data.table::fread("output/great_go/TCGA_AvsB/Apeaks_AB.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_B_AB <- data.table::fread("output/great_go/TCGA_AvsB/Bpeaks_AB.tsv", header = F, skip = 1)[c(1:10), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
GREAT_A_AB$mlog10FDR <- -log10(GREAT_A_AB$FDR)
GREAT_B_AB$mlog10FDR <- -log10(GREAT_B_AB$FDR)
p3 <- ggplot(GREAT_A_AB, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 
p4 <- ggplot(GREAT_B_AB, aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() 
pdf("output/Plots/S7_TCGA_GREAT_AvsB.pdf", height = 3, width = 10)
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
pdf("output/Plots/S7_Volcano_AvsB.pdf", width = 5, height = 5.5)
p5
dev.off()

nearby_Apeaks_AB <- rowRanges(TCGA_BRCA_ATAC_ER_se_cdp)[Apeaks_AB]$SYMBOL %>% table %>% sort(., decreasing = T)
nearby_Bpeaks_AB <- rowRanges(TCGA_BRCA_ATAC_ER_se_cdp)[Bpeaks_AB]$SYMBOL %>% table %>% sort(., decreasing = T)
write.csv(data.frame(gene = names(nearby_Apeaks_AB), number = as.numeric(nearby_Apeaks_AB)), "output/nearbygenes/TCGA_AvsB/nearby_Apeaks.csv", row.names = F, quote = F)
write.csv(data.frame(gene = names(nearby_Bpeaks_AB), number = as.numeric(nearby_Bpeaks_AB)), "output/nearbygenes/TCGA_AvsB/nearby_Bpeaks.csv", row.names = F, quote = F)

df <- data.frame(gene = names(nearby_Apeaks_AB), number = as.numeric(nearby_Apeaks_AB), rank = 1:length(nearby_Apeaks_AB))
p6 <- ggplot(df, aes(x = rank, y = number)) + geom_point(size = 0.5) + ArchR::theme_ArchR() + ggtitle("CA-A specific") +
  scale_y_continuous(breaks = seq(0,20,5))
df <- data.frame(gene = names(nearby_Bpeaks_AB), number = as.numeric(nearby_Bpeaks_AB), rank = 1:length(nearby_Bpeaks_AB))
p7 <- ggplot(df, aes(x = rank, y = number)) + geom_point(size = 0.5) + ArchR::theme_ArchR() + ggtitle("CA-B specific") +
  scale_y_continuous(breaks = seq(0,60,10))
pdf("output/Plots/S7_RankPlot_nearbygenes_AvsB.pdf", width = 3, height = 3)
p6
p7
dev.off()

#----- RNA-seq -----#
TCGA_BRCA_RNA_PT_se <- readRDS("rds/TCGA_BRCA_RNA_PT_se.rds")
TCGA_BRCA_RNA_ER_se <- TCGA_BRCA_RNA_PT_se[, TCGA_BRCA_ATAC_ER_se_cdp$submitter_id]
TCGA_BRCA_RNA_ER_se$CA <- TCGA_BRCA_ATAC_ER_se_cdp$CA

df <- assays(TCGA_BRCA_RNA_ER_se)$tpm_unstrand[which(mcols(TCGA_BRCA_RNA_ER_se)$gene_name %in% c("ESR1", "FOXA1")),]
df <- log2(df+1) %>% t %>% data.frame
colnames(df) <- c("ESR1", "FOXA1")
df$CA <- TCGA_BRCA_RNA_ER_se$CA

p8 <- lapply(c("ESR1", "FOXA1"), function(x){
  mtx <- df[,c(x, "CA")]
  colnames(mtx) <- c("value", "CA")
  p <- ggplot(mtx, aes(x = CA, y = value, fill = CA)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    theme_classic() + theme(legend.position = "none") + 
    scale_fill_manual(values = c("A" = "#1B9E77", "B" = "#D95F02", "ILC" = "#E7298A")) + 
    ggsignif::geom_signif(comparisons = list(c("A", "B"), c("A", "ILC"), c("B", "ILC")), test = "wilcox.test") + 
    ggtitle(x)
  return(p)
})
pdf("output/Plots/S7_TCGA_GEx_keyTF.pdf", width = 3, height = 6)
p8
dev.off()

#differential
assays(TCGA_BRCA_RNA_ER_se)$counts <- assays(TCGA_BRCA_RNA_ER_se)$unstranded
GEx_DEG_AvsB <- edgeR_pairwise(TCGA_BRCA_RNA_ER_se, compareCol = "CA", topGroup = "A", bottomGroup = "B")
A_DEG <- rownames(GEx_DEG_AvsB)[which(assay(GEx_DEG_AvsB)[,"log2FoldChange"] > 1 & assay(GEx_DEG_AvsB)[,"FDR"] < 0.01)]
B_DEG <- rownames(GEx_DEG_AvsB)[which(assay(GEx_DEG_AvsB)[,"log2FoldChange"] < -1 & assay(GEx_DEG_AvsB)[,"FDR"] < 0.01)]

DEG_mtx <- assays(GEx_DEG_AvsB)$differential %>% data.frame
DEG_mtx$gene <- mcols(GEx_DEG_AvsB)$gene_name
write.csv(DEG_mtx, "output/Tables/TCGA_GEx_DEG_AvsB.csv", quote = F)
#mcols(TCGA_BRCA_RNA_ER_se)[B_DEG, "gene_name"] %>% sort
df <- assay(GEx_DEG_AvsB)[, c("log2FoldChange", "FDR")] %>% data.frame
df$mlog10FDR <- -log10(df$FDR)
df$signif <- "N"
df$signif[which(df$log2FoldChange > 1 & df$FDR < 0.01)] <- "A"
df$signif[which(df$log2FoldChange < -1 & df$FDR < 0.01)] <- "B"
df$signif <- factor(df$signif, levels = c("N", "A", "B"))
df <- df[order(df$signif),]
p9 <- ggplot(df, aes(x = log2FoldChange, y = mlog10FDR, color = signif)) + geom_point(size = 0.8) + ArchR::theme_ArchR() +
  geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + scale_color_manual(values = cluster_colors2) + 
  ggtitle("Total peaks") + labs(x = "log2FC", y = "-log10(FDR)") + scale_y_continuous(breaks = seq(0,10,2))
pdf("output/Plots/S7_GEx_Volcano_AvsB.pdf", width = 5, height = 5.5)
p9
dev.off()

#----- prognosis -----#
library(survminer)
library(survival)
surv <- data.frame(sample = colnames(TCGA_BRCA_RNA_ER_se),
                   daysFO = TCGA_BRCA_RNA_ER_se$days_to_last_follow_up,
                   daysDE = TCGA_BRCA_RNA_ER_se$days_to_death,
                   status = TCGA_BRCA_RNA_ER_se$vital_status,
                   group  = TCGA_BRCA_RNA_ER_se$CA)
surv$days <- surv$daysFO
surv$days[is.na(surv$days)] <- surv$daysDE[is.na(surv$days)]
surv$status2 <- 0
surv$status2[surv$status == "Dead"] <- 1 

sf <- survfit(Surv(surv$days, surv$status2)~surv$group)
p10 <- ggsurvplot(fit = sf, data = surv,
                palette = c("#1B9E77", "#D95F02", "#E7298A"),
                pval = TRUE, pval.method = TRUE,
                risk.table = TRUE, conf.int = FALSE,
                ncensor.plot = FALSE, size = 1.5, #linetype = c(1, 3),
                title = "CA",
                legend.title = "",
                log.rank.weights = "1",
                risk.table.title = "",
                risk.table.y.text.col = TRUE,
                risk.table.y.text = FALSE)
pdf("output/Plots/S7_KMplot_CA.pdf", width = 6, height = 6)
p10
dev.off()

#----- sarrogate marker -----#
#mtx <- assay(GEx_DEG_AvsB)[which(assay(GEx_DEG_AvsB)[,"log2FoldChange"] < -1 & assay(GEx_DEG_AvsB)[,"FDR"] < 0.01),]
#mtx <- mtx[order(mtx[, "FDR"]),] %>% data.frame
#mtx$gene <- mcols(GEx_DEG_AvsB)[rownames(mtx), "gene_name"]
TCGA_BRCA_RNA_ERALL_se <- TCGA_BRCA_RNA_PT_se[, which(TCGA_BRCA_RNA_PT_se$ER == "Positive" & TCGA_BRCA_RNA_PT_se$HER2 == "Negative")]

markerValue <- assays(TCGA_BRCA_RNA_ERALL_se)$tpm_unstrand[B_DEG,]
markerValue <- log2(markerValue+1) %>% colMeans
plot(sort(markerValue))
summary(markerValue)

df <- data.frame(number = sort(markerValue), rank = c(1:length(markerValue)))
p11 <- ggplot(df, aes(x = rank, y = number)) + geom_point(size = 0.8) + ArchR::theme_ArchR() + ggtitle("Surrogate maker exp") +
       labs(x = "Rank", y = "Average Expression [log2(TPM)+1]") + scale_y_continuous(breaks = seq(0, 1, 0.2)) + geom_hline(yintercept = 0.4, col = "red", lty = "dashed")
pdf("output/Plots/S7_GEx_Surrogate_AverageExp.pdf", width = 4, height = 4)
p11
dev.off()

#cutoff 0.4
BOTTOM <- names(markerValue)[markerValue <= 0.4]
TOP <- names(markerValue)[markerValue > 0.4]

surv <- data.frame(sample = colnames(TCGA_BRCA_RNA_ERALL_se),
                   daysFO = TCGA_BRCA_RNA_ERALL_se$days_to_last_follow_up,
                   daysDE = TCGA_BRCA_RNA_ERALL_se$days_to_death,
                   status = TCGA_BRCA_RNA_ERALL_se$vital_status)
surv <- surv[surv$sample %in% c(BOTTOM, TOP), ]
surv$group <- "N"
surv$group[surv$sample %in% BOTTOM] <- "BOTTOM"
surv$group[surv$sample %in% TOP] <- "TOP"

surv <- surv[!duplicated(surv$sample),]

surv$days <- surv$daysFO
surv$days[is.na(surv$days)] <- surv$daysDE[is.na(surv$days)]
surv$status2 <- 0
surv$status2[surv$status == "Dead"] <- 1 
sf <- survfit(Surv(surv$days, surv$status2)~surv$group)
p12 <- ggsurvplot(fit = sf, data = surv,
                  palette = c("blue", "red"),
                  pval = TRUE, pval.method = TRUE,
                  risk.table = TRUE, conf.int = FALSE,
                  ncensor.plot = FALSE, size = 1.5, #linetype = c(1, 3),
                  title = "GEx Surrogate",
                  legend.title = "",
                  log.rank.weights = "1",
                  risk.table.title = "",
                  risk.table.y.text.col = TRUE,
                  risk.table.y.text = FALSE)
pdf("output/Plots/S7_KMplot_GExSurrogate.pdf", width = 6, height = 6)
p12
dev.off()

#----- METABRIC data -----#
# BiocManager::install("cBioPortalData")
# library(cBioPortalData)
# METABRIC <- cBioDataPack("brca_metabric",use_cache = TRUE,
#                          names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene"), 
#                          ask = TRUE)
# METABRIC.assay <- assays(METABRIC)
# METAClin <- colData(METABRIC)
# 
# #ER+/HER2-
# METAClin_ER <- METAClin[which(METAClin$ER_STATUS == "Positive" & METAClin$HER2_STATUS == "Negative"),]
# METAOS_HT <- as.data.frame(METAClin_ER[which(METAClin_ER$HORMONE_THERAPY == "YES"), c("OS_STATUS","OS_MONTHS")])
# METAOS_NHT <- as.data.frame(METAClin_ER[which(METAClin_ER$HORMONE_THERAPY == "NO"), c("OS_STATUS","OS_MONTHS")])
# 
# #Zscore matrix
# mtx <- METABRIC.assay$mrna_agilent_microarray
# ht_mtx <- mtx[, intersect(rownames(METAOS_HT), colnames(mtx))]
# nht_mtx <- mtx[, intersect(rownames(METAOS_NHT), colnames(mtx))]
# 
# #markergene
# markers <- read.table("output/Tables/TCGA_GEx_DEG_AvsB_CA-B.txt")[,1] #114 genes
# markers <- intersect(rownames(mtx), markers) #33 genes
# write.table(markers, "output/Tables/TCGA_GEx_DEG_AvsB_CA-B_overlapMETABRIC.txt", quote = F, row.names = F, col.names = F)
# 
# sort(colMeans(ht_mtx[markers, ])) %>% plot
# 
# #survival, hormone treatment
# markerValue <- colMeans(ht_mtx[markers, ]) #968 patients
# 
# #BOTTOM <- names(sort(markerValue))[c(1:319)]
# #TOP    <- names(sort(markerValue, decreasing = T))[c(1:319)]
# BOTTOM <- names(markerValue)[markerValue < 5.8]
# TOP    <- names(markerValue)[markerValue >= 5.8]
# 
# surv <- data.frame(sample = rownames(METAOS_HT),
#                    OS_MONTHS = METAOS_HT$OS_MONTHS,
#                    OS_STATUS = stringr::str_split(METAOS_HT$OS_STATUS, ":", simplify = T)[,1] %>% as.numeric)
# surv <- surv[surv$sample %in% c(BOTTOM, TOP), ]
# surv$group <- "N"
# surv$group[surv$sample %in% BOTTOM] <- "BOTTOM"
# surv$group[surv$sample %in% TOP] <- "TOP"
# 
# sf <- survfit(Surv(surv$OS_MONTHS, surv$OS_STATUS)~surv$group)
# p13 <- ggsurvplot(fit = sf, data = surv,
#                   palette = c("blue", "red"),
#                   pval = TRUE, pval.method = TRUE,
#                   risk.table = TRUE, conf.int = FALSE,
#                   ncensor.plot = FALSE, size = 1.5, #linetype = c(1, 3),
#                   title = "GEx Surrogate",
#                   legend.title = "",
#                   log.rank.weights = "1",
#                   risk.table.title = "",
#                   risk.table.y.text.col = TRUE,
#                   risk.table.y.text = FALSE)
# pdf("output/Plots/S7_KMplot_GExSurrogate.pdf", width = 6, height = 6)
# p12
# dev.off()

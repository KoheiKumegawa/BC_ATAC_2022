#------------------------------------------------------------------------------
# 13 - CompareCREs_FOXA1
#------------------------------------------------------------------------------
library(data.table)
library(SummarizedExperiment)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(ggplot2)

fr1  <- fread("output/output_bed/AvsB/Apeaks_AB_output.bed") %>% data.frame(.)
fr2  <- fread("output/output_bed/TCGA_AvsB/Apeaks_AB_output.bed") %>% data.frame(.)
fr3  <- fread("output/output_bed/AvsB/Bpeaks_AB_output.bed") %>% data.frame(.)
fr4  <- fread("output/output_bed/TCGA_AvsB/Bpeaks_AB_output.bed") %>% data.frame(.)
fr5  <- fread("ref/ERX008600.05.FOXA1.MCF7.bed") %>% data.frame(.)
fr6  <- fread("ref/ERX008605.05.FOXA1.T47D.bed") %>% data.frame(.)

JFCR_Apeaks <- GRanges(seqnames = fr1[,1], IRanges(start = fr1[,2]+1, end = fr1[,3]))
TCGA_Apeaks <- GRanges(seqnames = fr2[,1], IRanges(start = fr2[,2]+1, end = fr2[,3]))
JFCR_Bpeaks <- GRanges(seqnames = fr3[,1], IRanges(start = fr3[,2]+1, end = fr3[,3]))
TCGA_Bpeaks <- GRanges(seqnames = fr4[,1], IRanges(start = fr4[,2]+1, end = fr4[,3]))
MCF7_FOXA1 <- GRanges(seqnames = fr5[,1], IRanges(start = fr5[,2]+1, end = fr5[,3]))
T47D_FOXA1 <- GRanges(seqnames = fr6[,1], IRanges(start = fr6[,2]+1, end = fr6[,3]))

pdf("output/Plots/S7_overlaps_CA-A_peaks.pdf", width = 5, height = 5)
vennplot(list(JFCR = JFCR_Apeaks, TCGA = TCGA_Apeaks))
vennplot(list(JFCR = JFCR_Apeaks, MCF7_FOXA1 = MCF7_FOXA1, T47D_FOXA1 = T47D_FOXA1))
vennplot(list(TCGA = TCGA_Apeaks, MCF7_FOXA1 = MCF7_FOXA1, T47D_FOXA1 = T47D_FOXA1))
dev.off()

pdf("output/Plots/S7_overlaps_CA-B_peaks.pdf", width = 5, height = 5)
vennplot(list(JFCR = JFCR_Bpeaks, TCGA = TCGA_Bpeaks))
vennplot(list(JFCR = JFCR_Bpeaks, MCF7_FOXA1 = MCF7_FOXA1, T47D_FOXA1 = T47D_FOXA1))
vennplot(list(TCGA = TCGA_Bpeaks, MCF7_FOXA1 = MCF7_FOXA1, T47D_FOXA1 = T47D_FOXA1))
dev.off()

### bedtools fisher
pval_JFCR_T47D <- data.frame(CA = c("CA-A", "CA-B"), 
                             Signif = c(phyper(1119 - 1, 2226, 4426665 - 2226, 39642, lower.tail=F, log.p = T) / 2.303,
                                        phyper(175 - 1, 4293, 4426876 - 4293, 39757, lower.tail=F, log.p = T)/ 2.303))
pval_JFCR_MCF7 <- data.frame(CA = c("CA-A", "CA-B"), 
                             Signif = c(phyper(1250 - 1, 2226, 4583003 - 2226, 61862, lower.tail=F, log.p = T)/ 2.303,
                                        phyper(172 - 1, 4293, 4583234 - 4293, 62110, lower.tail=F, log.p = T)/ 2.303))
pval_TCGA_T47D <- data.frame(CA = c("CA-A", "CA-B"), 
                             Signif = c(phyper(3222 - 1, 5269, 4420707 - 5269, 39749, lower.tail=F, log.p = T)/ 2.303,
                                        phyper(148 - 1, 9830, 4420770 - 9830, 39763, lower.tail=F, log.p = T)/ 2.303))
pval_TCGA_MCF7 <- data.frame(CA = c("CA-A", "CA-B"), 
                             Signif = c(phyper(3096 - 1, 5269, 4576616 - 5269, 62069, lower.tail=F, log.p = T)/ 2.303,
                                        phyper(289 - 1, 9830, 4576739 - 9830, 62121, lower.tail=F, log.p = T)/ 2.303))

p1 <- ggplot(pval_JFCR_T47D, aes(x = CA, y = -Signif, fill = CA)) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) + ggtitle("JFCR peaks vs T47D FOXA1") + labs(x = "CA", y = "-log10(P-value)")
p2 <- ggplot(pval_JFCR_MCF7, aes(x = CA, y = -Signif, fill = CA)) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) + ggtitle("JFCR peaks vs MCF7 FOXA1") + labs(x = "CA", y = "-log10(P-value)")
p3 <- ggplot(pval_TCGA_T47D, aes(x = CA, y = -Signif, fill = CA)) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) + ggtitle("TCGA peaks vs T47D FOXA1") + labs(x = "CA", y = "-log10(P-value)")
p4 <- ggplot(pval_TCGA_MCF7, aes(x = CA, y = -Signif, fill = CA)) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) + ggtitle("TCGA peaks vs MCF7 FOXA1") + labs(x = "CA", y = "-log10(P-value)")

pdf("output/Plots/S7_overlaps_CA_significant_barplot.pdf", width = 2, height = 4)
p1
p2
p3
p4
dev.off()

# JFCR and TCGA CA-A
phyper(405 - 1, 5269, 3193319 - 5269, 2226, lower.tail=F, log.p = T) / 2.303
phyper(360 - 1, 9830, 3193319 - 9830, 4293, lower.tail=F, log.p = T) / 2.303


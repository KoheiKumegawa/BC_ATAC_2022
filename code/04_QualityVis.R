#------------------------------------------------------------------------------
# 04 - QualityVis.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)

QualityMetricsATAC <- readRDS("rds/QualityMetricsATAC.rds")

#tss score
tssscore <- lapply(QualityMetricsATAC, function(i) i$enrichScore) %>% unlist
tssscore <- data.frame(score = tssscore, tumor = factor(names(tssscore), levels = names(tssscore)))
p1 <- ggplot(tssscore, aes(x = tumor, y = score)) + geom_bar(stat = "identity") + ArchR::theme_ArchR() +
      geom_hline(yintercept = 5, color = "red", lty = "dashed") + labs(x = "Tumor", y = "TSS enrichment score") +
      theme(axis.text.x = element_text(angle = 90))

#normalized insertions
df <- QualityMetricsATAC$P1$plotDF
p2 <- ggplot(df, aes(x = distance, y = normInsert, color = sample)) + geom_line(size = 1) +
  theme_classic() + labs(x = "Distance From Center (bp)", y = "Normalized Insertion Profile") +
  scale_color_manual(values = "black") + ArchR::theme_ArchR() +
  scale_y_continuous(limits = c(0, max(df$normInsert)*1.05), expand = c(0,0)) +
  scale_x_continuous(limits = c(min(df$distance), max(df$distance)), expand = c(0,0)) + theme(legend.position = "none")

#fragment width
df <- data.frame(l = QualityMetricsATAC$P1$fragmentWidth, pt = "P1")
p3 <- ggplot(df, aes(x = l, color = pt)) + geom_line(stat = "density", size = 1) + ArchR::theme_ArchR() +
  scale_color_manual(values = "black") + xlim(0, 600) + ylab("Density") + xlab("Size of fragments (bp)") + theme(legend.position = "none")

pdf("output/Plots/QC_TSSenrichment.pdf", width = 6, height = 3)
p1
dev.off()
pdf("output/Plots/QC_NormalizedInsertions.pdf", width = 3, height = 3)
p2
dev.off()
pdf("output/Plots/QC_FragmentWidth.pdf", width = 3, height = 3)
p3
dev.off()
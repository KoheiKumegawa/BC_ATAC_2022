#------------------------------------------------------------------------------
# 10 - PrepareTCGADataSet
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(data.table)
library(rtracklayer)

#----- ATAC-seq -----#
counts   <- fread("data/TCGA/BRCA_raw_counts.txt")
log2norm <- fread("data/TCGA/BRCA_log2norm.txt")
peaks <- GRanges(seqnames = counts$seqnames, IRanges(start = counts$start, end = counts$end), name = counts$name)
counts <- as.matrix(counts[,-c(1:5)])
log2norm <- as.matrix(log2norm[,-c(1:5)])
TCGA_BRCA_ATAC_PT_se <- SummarizedExperiment(assays = list(counts = counts, normcounts = log2norm), rowRanges = peaks)
  
#patient info from Science paper
tcgaPatient <- read.csv("ref/TCGA_patient.csv", row.names = 1)
colnames(TCGA_BRCA_ATAC_PT_se) <- stringr::str_split(string = colnames(TCGA_BRCA_ATAC_PT_se), pattern = "_P", simplify = T)[,1]
colData(TCGA_BRCA_ATAC_PT_se) <- DataFrame(tcgaPatient[colnames(TCGA_BRCA_ATAC_PT_se),])
TCGA_BRCA_ATAC_UQ_se <- TCGA_BRCA_ATAC_PT_se[, TCGA_BRCA_ATAC_PT_se$Technical_Replicate_Number == "T1"]

#patient info from Xena
xena <- fread("ref/TCGA.BRCA.sampleMap-BRCA_clinicalMatrix")
xenaDF <- data.frame(row.names = xena$sampleID, histological_type = xena$histological_type)
TCGA_BRCA_ATAC_UQ_se$histological_type <- xenaDF[paste0(stringr::str_split(string = TCGA_BRCA_ATAC_UQ_se$Tissue_Barcode, pattern = "-01", simplify = T)[,1], "-01"),]

#----- Download RNA-seq data -----#
setwd("data/")
query <- GDCquery(project="TCGA-BRCA",
                  data.category="Transcriptome Profiling",
                  data.type="Gene Expression Quantification",
                  workflow.type="STAR - Counts")
GDCdownload(query)
TCGA_BRCA_RNA_se <- GDCprepare(query)
TCGA_BRCA_RNA_PT_se <- TCGA_BRCA_RNA_se[, which(TCGA_BRCA_RNA_se$sample_type == "Primary Tumor")]

#----- Download clincal info -----#
query2 <-  GDCquery(project = "TCGA-BRCA", 
                    data.category = "Clinical",
                    data.type = "Clinical Supplement", 
                    data.format = "BCR Biotab")
GDCdownload(query2)
clinical_BCRtab <- GDCprepare(query2)

IHC_status <- data.frame(row.names = clinical_BCRtab$clinical_patient_brca$bcr_patient_barcode,
                         ER = clinical_BCRtab$clinical_patient_brca$er_status_by_ihc,
                         PR = clinical_BCRtab$clinical_patient_brca$pr_status_by_ihc,
                         HER2 = clinical_BCRtab$clinical_patient_brca$her2_status_by_ihc)
IHC_status <- IHC_status[-c(1,2),]
colData(TCGA_BRCA_RNA_PT_se) <- cbind(colData(TCGA_BRCA_RNA_PT_se), IHC_status[colData(TCGA_BRCA_RNA_PT_se)$patient,])
setwd("../")

#----- saveRDS -----#
saveRDS(TCGA_BRCA_ATAC_UQ_se, "rds/TCGA_BRCA_ATAC_UQ_se.rds")
saveRDS(TCGA_BRCA_RNA_PT_se, "rds/TCGA_BRCA_RNA_PT_se.rds")
saveRDS(clinical_BCRtab, "rds/clinical_BCRtab.rds")





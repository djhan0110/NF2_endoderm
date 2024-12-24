# Set working environment -------------------------------------------------
rm(list = ls())
cat("\014")
setwd("C:/Workspace/R/NF2_Endo/RNA-Seq")

# Load packages -----------------------------------------------------------
library(readr)
library(tidyverse)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggforce)

# Load count matrix -------------------------------------------------------
NF2_Endo_RCs <- read.csv("00_Raw_data/NF2_Endo_RCs.csv")
clean_df <- janitor::clean_names(NF2_Endo_RCs)

clean_df <- clean_df %>%
  column_to_rownames("id") %>%
  rename(ctrl_n1 = "WT-1", ctrl_n2 = "WT-2", ctrl_n3 = "WT-3",
         ko_n1 = "KO-1", ko_n2 = "KO-2", ko_n3 = "KO-3",
         kod_n1 = "KOD-1", kod_n2 = "KOD-2", kod_n3 = "KOD-3")

head(clean_df)

genotype <- factor(rep(c("WT", "KO", "KOD"), each = 3)) %>% relevel(ref = "WT")

info_df <- data.frame(genotype)
rownames(info_df) <- colnames(clean_df)
info_df

clean_df[rownames(clean_df) %in% "ENSG00000186575", ]

# DESeq2 -------------------------------------------------------------
ddsm <- DESeqDataSetFromMatrix(countData = clean_df,
                               colData = info_df,
                               design = ~ genotype)

keep_rows <- rowSums(counts(ddsm)) >= ncol(ddsm)
keep_ddsm <- ddsm[keep_rows, ]
dds1 <- DESeq(keep_ddsm)
vst_dds1 <- vst(dds1, blind = F)
annotation_vst1 <- as.data.frame(assay(vst_dds1)) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(vst_dds1), keytype = "ENSEMBL", column = "SYMBOL")) 
annotation_vst1[annotation_vst1$hgnc_symbol %in% "NF2", ]
write.csv(annotation_vst1, "01_Export/NF2_Endo_VST.csv")

# Ref = WT ----------------------------------------------------------------
resultsNames(dds1)
resultsNames(dds1)[2]

LFC1 <- lfcShrink(dds1, coef = resultsNames(dds1)[2], type = "apeglm")
annotation_LFC1 <- as.data.frame(LFC1) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(LFC1), keytype = "ENSEMBL", column = "SYMBOL"))
annotation_LFC1[annotation_LFC1$hgnc_symbol %in% "NF2", ]
write.csv(annotation_LFC1, "01_Export/NF2_Endo_LFC_KO_vs_WT.csv")

resultsNames(dds1)[3]
LFC2 <- lfcShrink(dds1, coef = resultsNames(dds1)[3], type = "apeglm")
annotation_LFC2 <- as.data.frame(LFC2) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(LFC2), keytype = "ENSEMBL", column = "SYMBOL"))
annotation_LFC2[annotation_LFC2$hgnc_symbol %in% "NF2", ]
write.csv(annotation_LFC2, "01_Export/NF2_Endo_LFC_KOD_vs_WT.csv")

# Ref = KO ----------------------------------------------------------------
dds1$genotype <- relevel(genotype, "KO")
dds1$genotype
dds2 <- DESeq(dds1)
resultsNames(dds2)
resultsNames(dds2)[3]

LFC3 <- lfcShrink(dds2, coef = resultsNames(dds2)[3], type = "apeglm")
annotation_LFC3 <- as.data.frame(LFC3) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(LFC3), keytype = "ENSEMBL", column = "SYMBOL"))
annotation_LFC3[annotation_LFC3$hgnc_symbol %in% "NF2", ]
write.csv(annotation_LFC3, "01_Export/NF2_Endo_LFC_KOD_vs_KO.csv")

# LRT ---------------------------------------------------------------------
dds_LRT <- DESeq(keep_ddsm, test="LRT", reduced = ~ 1)
vst_LRT <- vst(dds_LRT, blind = FALSE)
annotated_vst_LRT <- as.data.frame(assay(vst_LRT))
annotated_vst_LRT <- annotated_vst_LRT %>%
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(annotated_vst_LRT), keytype = "ENSEMBL", column = "SYMBOL"))
write.csv(annotated_vst_LRT, "01_Export/NF2_Endo_VST_LRT.csv")
# saveRDS(deseq, file = "01_Export/NF2_Endo_LRT.rds")
# deseq <- readRDS("01_Export/NF2_Endo_LRT.rds")

resultsNames(dds_LRT)
resultsNames(dds_LRT)[2]
LFC_LRT2 <- lfcShrink(dds_LRT, coef = resultsNames(dds_LRT)[2], type = "apeglm")
annotated_LFC_LRT2 <- as.data.frame(LFC_LRT2)
annotated_LFC_LRT2 <- annotated_LFC_LRT2 %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(annotated_LFC_LRT2), keytype = "ENSEMBL", column = "SYMBOL"))
write.csv(annotated_LFC_LRT2, "01_Export/NF2_Endoderm_LFC_LRT_KO_vs_WT.csv")

resultsNames(dds_LRT)[3]
LFC_LRT3 <- lfcShrink(dds_LRT, coef = resultsNames(dds_LRT)[3], type = "apeglm")
annotated_LFC_LRT3 <- as.data.frame(LFC_LRT3)
annotated_LFC_LRT3 <- annotated_LFC_LRT3 %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(annotated_LFC_LRT3), keytype = "ENSEMBL", column = "SYMBOL"))
write.csv(annotated_LFC_LRT3, "01_Export/NF2_Endoderm_LFC_LRT_KOD_vs_WT.csv")

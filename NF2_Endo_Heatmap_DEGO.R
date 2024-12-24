# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/NF2_Endo/RNA-Seq")

# Load packages -----------------------------------------------------------
library(readr)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)

# Define the functions ----------------------------------------------------
clean_data <- function(data) {
  cleaned_data <- data %>%
    janitor::clean_names() %>%
    na.omit() %>%
    mutate(pattern = case_when(
      padj < 0.05 & log2fold_change > 0 ~ "Up",
      padj < 0.05 & log2fold_change < 0 ~ "Down",
      TRUE ~ "ns"),
      nlogpadj = -log10(padj)) %>%
    arrange(desc(log2fold_change))
  return(cleaned_data)
}

# Load count matrix -------------------------------------------------------
NF2_Endo_LFC <- read_csv("01_Export/NF2_Endo_LFC.csv")
NF2_Endo_VST <- read_csv("01_Export/NF2_Endo_VST.csv")
LFC_clean <- clean_data(NF2_Endo_LFC)
VST_signif <- LFC_clean %>%
  left_join(NF2_Endo_VST, by = c("x1" = "...1")) %>% 
  filter(pattern == "Up")

NF2_Endo_BiNGO <- read_excel("01_Export/NF2_Endo_BiNGO.xlsx")
NF2_Endo_BiNGO$GO_id <- paste0("GO:", NF2_Endo_BiNGO$GO_id)

# Retrieve GO terms and gene information ----------------------------------
# GO:0007399 Nervous system development
# GO:0048738 Cardiac muscle tissue development
# GO:0007507 Heart development
# GO:0060537 Muscle tissue development
 
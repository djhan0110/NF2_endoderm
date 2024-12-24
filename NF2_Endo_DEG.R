# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
project_name <- "NF2_Endo"
setwd("C:/Workspace/R/NF2_Endo/RNA-Seq")

# Load packages -----------------------------------------------------------
library(readr)
library(tidyverse)
library(scales)
library(ggforce) 
library(ggh4x)

# Load count matrix -------------------------------------------------------
NF2_Endo_VST <- read_csv("01_Export/NF2_Endo_VST.csv")
NF2_Endo_LFC_KOD_vs_WT <- read_csv("01_Export/NF2_Endo_LFC_KOD_vs_WT.csv")
NF2_Endo_LFC_KOD_vs_KO <- read_csv("01_Export/NF2_Endo_LFC_KOD_vs_KO.csv")
NF2_Endo_LFC_KO_vs_WT <- read_csv("01_Export/NF2_Endo_LFC_KO_vs_WT.csv")

# Define Functions --------------------------------------------------------
padj_cutoff = 0.05
lfc_cutoff = 1

clean_data <- function(LFC, compare_set) {
  LFC %>% janitor::clean_names() %>%  na.omit() %>%
    mutate(pattern = case_when(padj < padj_cutoff & log2fold_change > lfc_cutoff ~ "Up",
                            padj < padj_cutoff & log2fold_change < lfc_cutoff ~ "Down",
                            TRUE ~ "ns"),
           nlogpadj = ifelse(padj != 0, -log10(padj), -log10(.Machine$double.xmin)),
           compare_set = compare_set) %>% 
    arrange(., desc(log2fold_change))
}

# Apply the functions to the data -----------------------------------------
KOD_WT <- clean_data(NF2_Endo_LFC_KOD_vs_WT, "KOD_WT")
KOD_KO <- clean_data(NF2_Endo_LFC_KOD_vs_KO, "KOD_KO")
KO_WT <- clean_data(NF2_Endo_LFC_KO_vs_WT, "KO_WT")

table(KOD_WT$pattern)
table(KOD_KO$pattern)
table(KO_WT$pattern)

split_KOD_WT <- split(KOD_WT, KOD_WT$pattern)
split_KOD_KO <- split(KOD_KO, KOD_KO$pattern) 
split_KO_WT <- split(KO_WT, KO_WT$pattern) 

topN_KOD_WT <- rbind(head(split_KOD_WT$Up, 10), tail(split_KOD_WT$Down, 10))
topN_KOD_KO <- rbind(head(split_KOD_KO$Up, 10), tail(split_KOD_KO$Down, 10))
topN_KOD_KO <- rbind(head(split_KOD_KO$Up, 10), tail(split_KOD_KO$Down, 10))
#Top_genes$pattern <- factor(Top_genes$pattern, levels = c("Up", "Down"))

# Create the plot ---------------------------------------------------------
DEG_plot <- ggplot(topN_KOD_WT, aes(x = log2fold_change, y = reorder(hgnc_symbol, abs(log2fold_change)))) +
  geom_col(aes(fill = pattern), position = position_dodge(0.5), width = 0.1, linewidth = 0.2, color = "black", alpha = 1) +
  geom_point(aes(size = nlogpadj, fill = pattern), position = position_dodge(0.5), shape = 21, color = "black", alpha = 1) +
  geom_vline(xintercept = 0) +
  scale_size_continuous(breaks = pretty_breaks(n = 5), range = c(0.5, 4)) +
  scale_fill_manual(values = c("steelblue", "firebrick")) +
  facet_wrap(~ pattern, scales = "free", ncol = 2, nrow = 1) +
  facetted_pos_scales(x = list(pattern == "Up" ~ scale_x_continuous(expand = expansion(mult = c(0, 0.2)),
                                                                    breaks = pretty_breaks(n = 3),
                                                                    labels = label_number(accuracy = 1)),
                               pattern == "Down" ~ scale_x_continuous(expand = expansion(mult = c(0.2, 0)),
                                                                      breaks = pretty_breaks(n = 3),
                                                                      labels = label_number(accuracy = 1))),
                      y = list(pattern == "Up" ~ scale_y_discrete(position = "right"), pattern == "Down" ~ scale_y_discrete(position = "left"))) +
  labs(x = expression(bold("log"["2"]*"FC")), y = element_blank(), 
       fill = "", size = expression(bold("-Log(Adj."*bolditalic(p)*")"))) +
  guides(fill = "none") +
  theme_classic() +
  theme(panel.spacing = unit(0, "lines"),
        text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(face = "bold", color = "black"),
        axis.text.y = element_text(face = "bold.italic", color = "black"),
        legend.position = "top",
        legend.text = element_text(color = "black", face = "bold"),
        legend.title = element_text(color = "black",face = "bold"),
        legend.margin = margin(0, 20, 0, 0),
        legend.key.spacing = unit(1, "mm"),
        legend.box.spacing = unit(1, "mm"),
        strip.background = element_blank(),
        strip.text = element_blank()) 

DEG_plot
plot_ratio = 1
ggsave("02_Figure/NF2_Endo_DEG.tiff", DEG_plot, units = "in", device = "tiff", dpi = 300,
       width = 3*plot_ratio, height = 3*plot_ratio)  

rescue_join <- KOD_KO %>% 
  left_join(KO_WT)
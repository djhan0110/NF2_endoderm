# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/NF2_Endo/RNA-Seq")

# Load required packages --------------------------------------------------
library(readr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(ggforce) 
library(ggh4x)
library(AnnotationDbi)
library(ComplexHeatmap)
rstatix 
# Define the functions ----------------------------------------------------
kegg_enrichment <- function(data) {
  enrichment <- enrichKEGG(gene = data$entrez_ids,
                           organism = "hsa",
                           keyType = "kegg",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
  result <- data.frame(enrichment@result) %>% 
    mutate(GR_converted = sapply(GeneRatio, function(x) eval(parse(text = x))),
           BG_converted = sapply(BgRatio, function(x) eval(parse(text = x))),
           nlogq = -log10(qvalue)) %>%
    arrange(desc(GR_converted))
  return(result)
}

# Load data and processing ------------------------------------------------
NF2_Endo_LFC_KO_vs_WT <- read_csv("01_Export/NF2_Endo_LFC_KO_vs_WT.csv")
LFC_clean <- NF2_Endo_LFC_KO_vs_WT %>% 
  janitor::clean_names() %>%
  na.omit() %>%
  mutate(pattern = case_when(log2fold_change > 1 & padj <= 0.05 ~ "Up",
                             log2fold_change < -1 & padj <= 0.05 ~ "Down",
                             TRUE ~ "ns"),
         nlogpadj = -log10(padj),
         entrez_ids = mapIds(org.Hs.eg.db, keys = x1, keytype = "ENSEMBL", column = "ENTREZID")) %>%
  arrange(desc(log2fold_change))

KEGG_up <- LFC_clean %>% 
  filter(pattern == "Up")
KEGG_up <- kegg_enrichment(KEGG_up) %>% 
  mutate(pattern = "Up") 

KEGG_down <- LFC_clean %>% 
  filter(pattern == "Down")
KEGG_down <- kegg_enrichment(KEGG_down) %>% 
  mutate(pattern = "Down")

KEGG_bind <- rbind(head(KEGG_up, 6), head(KEGG_down, 6))
KEGG_bind$pattern <- factor(KEGG_bind$pattern, levels = c("Up", "Down"))
KEGG_bind$labels <- str_wrap(KEGG_bind$Description, width = 50)
KEGG_bind <- KEGG_bind %>% 
  mutate(group = tidytext::reorder_within(labels, GR_converted, within = pattern))


# Plot --------------------------------------------------------------------
KEGG_plot <- ggplot(KEGG_bind, aes(x = GR_converted, y = group, fill = nlogq)) +
  facet_grid(pattern ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_col() +
  geom_text(aes(x = 0, y = group, label = str_wrap(Description, width = 45)),
            hjust = 0, nudge_x = 0.001, lineheight = 0.8) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_distiller(palette = "YlOrBr") +
  labs(title = expression(bold("KEGG pathways")),
    x = "Gene Ratio", y = "",
    fill = expression("-Log(" * italic("q") * "-value)")) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.text = element_text(color = "black"),
    strip.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.key.height = unit(0.25, "cm"),
    legend.position = "top"
  )
  
KEGG_plot

plot_ratio = 2
ggsave("02_Figure/NF2_Endo_KEGG.tiff", KEGG_plot, units = "in", device = "tiff", dpi = 300,
       width = 1.8*plot_ratio, height = 2.2*plot_ratio)


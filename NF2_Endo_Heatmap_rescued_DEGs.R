# Set working environment -------------------------------------------------
rm(list = ls())
setwd("C:/Workspace/R/NF2_Endo/RNA-seq")

# Load required packages --------------------------------------------------
library(readr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(seriation)
library(janitor)
library(reshape2)

# Filtering induced, suppressed gene list ---------------------------------
DEG_Endo_rescued <- read_csv("01_Export/NF2_Endo_rescued_DEGs.csv")
DEG_1 <- as.data.frame(clean_names(DEG_Endo_rescued))
head(DEG_1)
ind_1 <- DEG_1 %>% 
  filter(attr == 'ind') %>% 
  mutate(dist = abs(log2fold_change_x - log2fold_change_y))
ind_1 <- ind_1[order(ind_1$dist, decreasing = T), ]
top_ind <- head(ind_1, 200)

sup_1 <- DEG_1 %>% 
  filter(attr == 'sup') %>% 
  mutate(dist = abs(log2fold_change_x - log2fold_change_y))
sup_1 <- sup_1[order(sup_1$dist, decreasing = T), ]
top_sup <- head(sup_1, 200)

glist_1 <- data.frame(top_ind$ensg, "ind")
glist_2 <- data.frame(top_sup$ensg, "sup")

# Trimming transformed count number ---------------------------------------     
TFs_Endo_all <- read_csv("01_Export/NF2_Endo_TFs_all.csv")
head(TFs_Endo_all)
TFs_1 <- as.data.frame(clean_names(TFs_Endo_all)) %>% 
  na.omit(TFs_1)
head(TFs_1)
sel_1 <- TFs_1 %>% 
  filter(TFs_1$x1 %in% glist_1$top_ind.ensg)
df_1 <- sel_1[, c(2:10)]
row.names(df_1) <- sel_1$hgnc_symbol
head(df_1)
sel_2 <- TFs_1 %>% 
  filter(TFs_1$x1 %in% glist_2$top_sup.ensg)
df_2 <- sel_2[, c(2:10)]
row.names(df_2) <- sel_2$hgnc_symbol
head(df_2)
col_names <- c("Wi-1", "Wi-2", "Wi-3", 
               "-Dox-1", "-Dox-2", "-Dox-3", 
               "+Dox-1", "+Dox-2", "+Dox-3")
col_names <- factor(c("Wi-1", "Wi-2", "Wi-3", 
                      "-Dox-1", "-Dox-2", "-Dox-3", 
                      "+Dox-1", "+Dox-2", "+Dox-3"),
                    levels = c("-Dox-1", "-Dox-2", "-Dox-3", 
                    "+Dox-1", "+Dox-2", "+Dox-3",
                    "Wi-1", "Wi-2", "Wi-3"))
colnames(df_1) <- col_names
colnames(df_2) <- col_names

head(df_1)
head(df_2)

mat_1 <- as.matrix(df_1)
mat_2 <- as.matrix(df_2)

order_ind <- seriate(dist(mat_1), method = "TSP")
order_sup <- seriate(dist(mat_2), method = "TSP")

ind_cols <- colorRamp2(breaks = c(1, 5, 10), colors = c("gray80", "white", "red"))
sup_cols <- colorRamp2(breaks = c(1, 5, 10), colors = c("blue", "white", "gray80"))

fa_1 = factor(rep(letters[1:3], each = 3), levels = c("b", "c", "a"))
fa_2 = factor(rep(letters[1:3], each = 3), levels = c("b", "a", "c"))

# Plotting ----------------------------------------------------------------
hm_ind <- Heatmap(mat_1, col = ind_cols,
                 
                  show_row_dend = F,
                  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                  row_names_side = "left",
                  row_dend_reorder = TRUE,
                  
                  show_column_dend = F,
                  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
                  column_split = fa_1,
                  
                  heatmap_legend_param = list(title = expression(bold("rlog")),
                                              direction = "vertical",
                                              color_bar = "continuous",
                                              legend_height = unit(3, "cm"),
                                              title_gp = gpar(fontsize = 12, fontface = "bold"),
                                              labels_gp = gpar(fontsize = 12, fontface = "bold"))
                  )
hm_ind
tiff("hm_ind.tiff", width = 6, height = 6, units = "in", res = 300)
draw(hm_ind)
dev.off()

hm_sup <- Heatmap(mat_2, col = sup_cols,
                  
                  show_row_dend = F,
                  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                  row_names_side = "left",
                  row_dend_reorder = TRUE,
                  
                  show_column_dend = F,
                  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
                  column_split = fa_2,
                  
                  heatmap_legend_param = list(title = expression(bold("rlog")),
                                              direction = "vertical",
                                              color_bar = "continuous",
                                              legend_height = unit(3, "cm"),
                                              title_gp = gpar(fontsize = 12, fontface = "bold"),
                                              labels_gp = gpar(fontsize = 12, fontface = "bold"))
)

hm_sup
tiff("hm_sup.tiff", width = 6, height = 6, units = "in", res = 300)
draw(hm_sup)
dev.off()

# Mast cell analysis - CellChat analysis
# This script performs cell-cell communication analysis using CellChat
# Author: Based on analysis from README.md  
# Date: 2024

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(ggpubr)
library(ggimage)

# Set working directory (modify as needed)
# setwd("path/to/your/analysis")

# Load data
# sce.in <- readRDS("path/to/totalCITE-seq.rds")
# MC <- readRDS("path/to/MC_object.rds")

# Define color palette
allcolour <- c(
  "#D0AFC4", "#89558D", "#AFC2D9", "#435B95", "#79B99D",
  "#D55640", "#E69F84", "#6CB8D2", "#479D88", "#415284",
  "#C6367A", "#ECDC52", "#D1352B", "#9B5B33"
)
names(allcolour) <- c('Mcpt9 high MC', 'Mcpt9 medium MC', 'Cycling MC', 'Nr4a1 high MC', 'Lrmda+ MC', 
                      "IE-hpD07", "IE-hpD14", "LP-hpD0", "LP-hpD07", "LP-hpD14", 
                      "BMCP", "GMP", "DN(P)", "DN(Q)")

# 1. Epithelial cells annotation ----

# Subset epithelial cells
epi <- subset(sce.in, celltype == "Epithelium")

# Standard Seurat processing pipeline
epi <- NormalizeData(epi) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

# Batch correction with Harmony
epi <- RunHarmony(epi, group.by.var = "group")

# Clustering
epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:7)
epi <- FindClusters(epi, resolution = 0.1)
table(epi@meta.data$seurat_clusters)

# UMAP embedding
epi <- RunUMAP(epi, reduction = "harmony", dims = 1:7, return.model = TRUE)

# Find epithelial cell markers
Idents(epi) <- "seurat_clusters"
epi_markers <- FindAllMarkers(
  object = epi,
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25,
  slot = "data"
)

epi_markers_filtered <- epi_markers %>% 
  dplyr::select(gene, everything()) %>% 
  subset(p_val < 0.05)

top20_epi <- epi_markers_filtered %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)

# 2. Cx3cr1 expression analysis in mast cells ----

# Basic violin plot for Cx3cr1 by group
p1 <- VlnPlot(MC, pt.size = 0, group.by = "group", features = "Cx3cr1") + 
  scale_fill_manual(values = allcolour) +
  geom_boxplot(width = 0.2, color = "black", 
               outlier.size = 0, 
               position = position_dodge(width = 0.75),
               aes(fill = NULL), 
               alpha = 0.2,
               size = 1) +
  ggtitle("Cx3cr1 expression by group")

# Violin plot by celltype
Idents(MC) <- "celltype"
p2 <- VlnPlot(MC, pt.size = 0, features = "Cx3cr1") + 
  scale_fill_manual(values = allcolour) +
  ggtitle("Cx3cr1 expression by celltype")

# Violin plot by tissue
Idents(MC) <- "tissue"
p3 <- VlnPlot(MC, pt.size = 0, features = "Cx3cr1") + 
  scale_fill_manual(values = allcolour) +
  ggtitle("Cx3cr1 expression by tissue")

print(p1)
print(p2)
print(p3)

# 3. Statistical analysis with t-tests ----

# Define comparisons for statistical testing
my_comparisons1 <- list(c("IE-hpD07", "IE-hpD14"))
my_comparisons2 <- list(c("LP-hpD07", "IE-hpD07"))
my_comparisons3 <- list(c("LP-hpD07", "LP-hpD14"))
my_comparisons4 <- list(c("LP-hpD0", "LP-hpD14"))

# Set identity for statistical testing
Idents(MC) <- "group"

# Create violin plot with statistical comparisons
p_stat <- VlnPlot(MC, features = "Cx3cr1", pt.size = 0) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = allcolour) +
  ggtitle("Cx3cr1 expression with statistical tests")

# Add statistical comparisons
p_stat_final <- p_stat + 
  stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  stat_compare_means(comparisons = my_comparisons2, method = "t.test") +
  stat_compare_means(comparisons = my_comparisons3, method = "t.test") +
  stat_compare_means(comparisons = my_comparisons4, method = "t.test")

print(p_stat_final)

# 4. Additional analysis functions ----

# Function to prepare data for CellChat analysis
prepare_cellchat_data <- function(seurat_obj, group_by = "celltype") {
  # Extract expression data
  data_input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  
  # Extract metadata
  labels <- seurat_obj@meta.data[[group_by]]
  names(labels) <- colnames(seurat_obj)
  
  return(list(data = data_input, labels = labels))
}

# Function to calculate Cx3cr1 expression statistics
analyze_cx3cr1_expression <- function(seurat_obj) {
  # Extract Cx3cr1 expression
  cx3cr1_expr <- GetAssayData(seurat_obj, features = "Cx3cr1", slot = "data")
  
  # Add to metadata
  seurat_obj$Cx3cr1_expression <- as.numeric(cx3cr1_expr)
  
  # Calculate statistics by group
  stats_by_group <- seurat_obj@meta.data %>%
    group_by(group) %>%
    summarise(
      mean_expr = mean(Cx3cr1_expression),
      median_expr = median(Cx3cr1_expression),
      sd_expr = sd(Cx3cr1_expression),
      n_cells = n(),
      .groups = 'drop'
    )
  
  # Calculate statistics by celltype  
  stats_by_celltype <- seurat_obj@meta.data %>%
    group_by(celltype) %>%
    summarise(
      mean_expr = mean(Cx3cr1_expression),
      median_expr = median(Cx3cr1_expression),
      sd_expr = sd(Cx3cr1_expression),
      n_cells = n(),
      .groups = 'drop'
    )
  
  return(list(
    by_group = stats_by_group,
    by_celltype = stats_by_celltype
  ))
}

# Run Cx3cr1 analysis
cx3cr1_stats <- analyze_cx3cr1_expression(MC)
print("Cx3cr1 expression statistics by group:")
print(cx3cr1_stats$by_group)
print("Cx3cr1 expression statistics by celltype:")
print(cx3cr1_stats$by_celltype)

# 5. Save results ----

# Save processed objects and results
save(list = c("epi", "epi_markers_filtered", "top20_epi", "cx3cr1_stats"),
     file = "cellchat_analysis_results.RData")

# Save plots
plots_list <- list(
  cx3cr1_by_group = p1,
  cx3cr1_by_celltype = p2,
  cx3cr1_by_tissue = p3,
  cx3cr1_with_stats = p_stat_final
)

save(plots_list, file = "cellchat_plots.RData")

print("CellChat analysis completed successfully!")
print(paste("Number of epithelial cells:", ncol(epi)))
print(paste("Number of epithelial clusters:", length(unique(epi$seurat_clusters))))
print(paste("Cx3cr1 expressing cells:", sum(MC$Cx3cr1_expression > 0)))
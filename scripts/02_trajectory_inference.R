# Mast cell analysis - Trajectory Inference
# This script contains various trajectory inference methods
# Author: Based on analysis from README.md
# Date: 2024

# Load required libraries
library(Seurat)
library(SingleCellExperiment)
library(destiny)
library(ggplot2)
library(plotly)
library(SeuratDisk)
library(dyno)
library(tidyverse)
library(Matrix)
library(slingshot)
library(monocle)
library(magrittr)
library(RColorBrewer)
library(reshape2)
library(Biobase)
library(ggsci)
library(ggpubr)
library(data.table)

# Set working directory (modify as needed)
# setwd("C:/Users/Taotao/OneDrive/桌面/data/2024.8.14 harmony analysis")

# Load the processed Mast cell object
# MC <- readRDS(file = "./MC_by_trace.rds")

# 2.1.1a) Diffusion Map by R ----

# Subset to relevant samples
MC_subset <- subset(MC, group %in% c("LP-hpD07", "LP-hpD14", "IE-hpD07", "IE-hpD14"))
Idents(MC_subset) <- "celltype"
sce <- as.SingleCellExperiment(MC_subset)

# Run DiffusionMap (may take a long time)
dm <- DiffusionMap(sce)

# Prepare data for visualization
cellLabels <- sce$ident
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  DC4 = eigenvectors(dm)[, 4],
                  Clusters = cellLabels)

# 2D diffusion map plot
ggplot(tmp, aes(x = DC1, y = DC2, colour = Clusters)) +
  geom_point()  + 
  xlab("DC1") + 
  ylab("DC2") +
  theme_classic()+scale_color_manual(values = allcolour)+theme_dr()+ 
          guides(color = guide_legend(override.aes = list(size=5)))+ theme(panel.grid=element_blank(),
          legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.key.size=unit(1,'cm'))+ggtitle("Diffusion Map for Mast Cell")+theme(plot.title = element_text(size = 20, hjust = 0.5))

# Interactive 3D diffusion map
p = plot_ly(x=tmp$DC1, y=tmp$DC2, z=tmp$DC3, type="scatter3d", mode="markers", color=tmp$Clusters, marker = list(size = 3 ), colors = allcolour) %>%
  layout(scene = list(
    aspectmode = 'cube',
    xaxis = list(
      title = list(text = 'DC1', font = list(size = 12, color = 'black')),
      tickvals = NULL,
      showgrid = TRUE,
      gridcolor = 'lightgrey',
      zeroline = FALSE,
      showline = TRUE,
      linecolor = 'black',
      linewidth = 4
    ),
    yaxis = list(
      title = list(text = 'DC2', font = list(size = 12, color = 'black')),
      tickvals = NULL,
      showgrid = TRUE,
      gridcolor = 'lightgrey',
      zeroline = FALSE,
      showline = TRUE,
      linecolor = 'black',
      linewidth = 4
    ),
    zaxis = list(
      title = list(text = 'DC3', font = list(size = 12, color = 'black')),
      tickvals = NULL,
      showgrid = TRUE,
      gridcolor = 'lightgrey',
      zeroline = FALSE,
      showline = TRUE,
      linecolor = 'black',
      linewidth = 4
    ),
    legend = list(
      font = list(size = 16)
  )))

htmlwidgets::saveWidget(as_widget(p), "Interactive3D.html", title = "Diffusion map")

# 2.1.1b) Diffusion Map by Python ----

# Convert Seurat object to h5ad format for Python
SaveH5Seurat(MC, filename="test.h5seurat", overwrite = TRUE)
Convert("test.h5seurat", dest = "h5ad", overwrite = TRUE)

# 2.1.2) dyno Analysis ----

# Prepare data for dyno
dataset <- wrap_expression(
  counts = t(MC@assays$RNA@counts),
  expression = t(MC@assays$RNA@data)
)

# Add prior information
dataset <- add_prior_information(
  dataset,
  start_id = "LP-hpD07_cite_GTGCAAGGT_AAGGTGGTA_AACAGGAAC"
)

# Add clustering information
dataset <- add_grouping(
   dataset,
   MC$celltype
)

# Choose the best method for dataset
guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected

# Run different trajectory inference methods
model_slingshot <- infer_trajectory(dataset, methods_selected[1])
model_paga_tree <- infer_trajectory(dataset, methods_selected[2])
model_scorpius <- infer_trajectory(dataset, methods_selected[3])
model_angle <- infer_trajectory(dataset, methods_selected[4])

# Visualize trajectory
model <- model_scorpius
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)

# 2.1.3) Slingshot Analysis ----

sce_slingshot <- as.SingleCellExperiment(MC, assay = "RNA")

sce_slingshot1 <- slingshot(sce_slingshot,
                     reducedDim = 'UMAP',
                     clusterLabels = sce_slingshot$celltype,
                     start.clus = 'Cycling MC',
                     approx_points = 150)

SlingshotDataSet(sce_slingshot1)

# 2.1.4) Monocle2 Analysis ----

set.seed(12345)

# Prepare data for Monocle2
expr_matrix <- MC@assays$RNA@counts
sample_sheet <- MC@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(MC))
rownames(gene_annotation) <- rownames(MC)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

# Normalization
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Filter low expression genes
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

# Find trajectory genes
diff_celltype <- differentialGeneTest(cds[expressed_genes,], 
                                     fullModelFormulaStr = "~celltype", 
                                     reducedModelFormulaStr = "~group", 
                                     relative_expr = TRUE,
                                     cores = 4)

# Select ordering genes
ordering_genes <- row.names(subset(diff_celltype, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)

# Reduce dimension and order cells
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

# Plot trajectory
plot_cell_trajectory(cds, color_by = "celltype")
plot_cell_trajectory(cds, color_by = "Pseudotime")

# Save results
# saveRDS(list(
#   diffusion_map = dm,
#   dyno_models = list(slingshot = model_slingshot, scorpius = model_scorpius),
#   slingshot = sce_slingshot1,
#   monocle2 = cds
# ), "trajectory_results.rds")
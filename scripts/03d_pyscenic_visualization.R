# Mast cell analysis - pySCENIC visualization
# This script visualizes the results from pySCENIC analysis
# Author: Based on analysis from README.md
# Date: 2024

# Load required libraries
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
library(ggrepel)

# Set working directory (modify as needed)
# setwd("path/to/your/analysis")

# Load Seurat object
# sce <- readRDS("./MC_by_trace.rds")

# Define color palette
allcolour <- c(
  "#D0AFC4", "#89558D", "#AFC2D9", "#435B95", "#79B99D",
  "#D55640", "#E69F84", "#6CB8D2", "#479D88", "#415284",
  "#C6367A", "#ECDC52", "#D1352B", "#9B5B33"
)
names(allcolour) <- c('Mcpt9 high MC', 'Mcpt9 medium MC', 'Cycling MC', 'Nr4a1 high MC', 'Lrmda+ MC', 
                      "IE-hpD07", "IE-hpD14", "LP-hpD0", "LP-hpD07", "LP-hpD14", 
                      "BMCP", "GMP", "DN(P)", "DN(Q)")

# 1. Extract pySCENIC loom information ----

loom <- open_loom("./sample_SCENIC.loom") 

# Extract regulon information
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)  

close_loom(loom)

# Match cells between Seurat and SCENIC
sub_regulonAUC <- regulonAUC[, match(colnames(sce), colnames(regulonAUC))]
dim(sub_regulonAUC)

# Verify cell order consistency
if (!identical(colnames(sub_regulonAUC), colnames(sce))) {
  stop("Cell order mismatch between SCENIC and Seurat objects")
}

# Prepare cell metadata
cellClusters <- data.frame(row.names = colnames(sce), 
                          seurat_clusters = as.character(sce$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(sce), 
                       celltype = sce$celltype)

# Save data for further analysis
save(sub_regulonAUC, cellTypes, cellClusters, sce,
     file = 'for_rss_and_visual.Rdata')

# 2. Visualize key transcription factors ----

# Focus on key Mast cell transcription factors: GATA2(+) and SOX4(+)
regulonsToPlot <- c('Gata2(+)', 'Sox4(+)')

# Check if regulons exist
if (!all(regulonsToPlot %in% rownames(sub_regulonAUC))) {
  missing_regulons <- regulonsToPlot[!regulonsToPlot %in% rownames(sub_regulonAUC)]
  warning("Missing regulons: ", paste(missing_regulons, collapse = ", "))
  regulonsToPlot <- regulonsToPlot[regulonsToPlot %in% rownames(sub_regulonAUC)]
}

# Add regulon activities to Seurat metadata
sce@meta.data <- cbind(sce@meta.data, t(assay(sub_regulonAUC[regulonsToPlot, ])))

# Create visualizations
p1 <- DotPlot(sce, features = unique(regulonsToPlot)) + RotatedAxis()
p2 <- RidgePlot(sce, features = regulonsToPlot, ncol = 2) 
p3 <- VlnPlot(sce, features = regulonsToPlot, pt.size = 0)
p4 <- FeaturePlot(sce, features = regulonsToPlot)

# Combine plots
combined_plot <- wrap_plots(p1, p2, p3, p4)
print(combined_plot)

# 3. TF activity analysis by cell groups ----

# Calculate average TF activities by cell type
selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes), cellTypes[, selectedResolution])

# Remove extended regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), ] 

# Calculate average expression per group
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                 function(cells) 
                                   rowMeans(getAUC(sub_regulonAUC)[, cells]))

# Scale expression (z-score normalization)
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                         center = TRUE, scale = TRUE))
regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

# Create heatmap
heatmap_colors <- colorRampPalette(c("#3d54a1", "#8b83bc", "#e4d9eb", 
                                    "#f8ddd3", "#e58065", "#df553f"))(2000)

tf_heatmap <- Heatmap(
  regulonActivity_byGroup_Scaled,
  name = "z-score",
  col = heatmap_colors,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE
)

draw(tf_heatmap)

# 4. RSS (Regulon Specificity Score) analysis ----

# Calculate RSS
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), 
               cellAnnotation = cellTypes[colnames(sub_regulonAUC), selectedResolution])
rss <- na.omit(rss)

# Interactive RSS plot
rssPlot <- plotRSS(rss)
print(plotly::ggplotly(rssPlot$plot))

# Enhanced RSS plot
rssPlot_enhanced <- plotRSS(
  rss,
  thr = 0.01,
  zThreshold = 1.2,
  cluster_columns = FALSE,
  order_rows = TRUE,
  varName = "cellType",
  col.low = '#330066',
  col.mid = '#66CC66',
  col.high = '#FFCC33'
)
print(rssPlot_enhanced)

# 5. Waterfall plots for top regulons ----

topN <- 5
plot_list <- list()

for (i in colnames(rss)) {
  df.i <- data.frame(SpecificityScore = rss[, i], labels = rownames(rss))
  df.i <- arrange(df.i, desc(SpecificityScore))
  df.i$Regulons <- 1:nrow(df.i)
  df.i$color <- ifelse(df.i$Regulons <= topN, "#E9E55A", "grey")
  df.i$labels[df.i$Regulons > topN] <- NA
  
  p <- ggplot(df.i, aes(Regulons, SpecificityScore)) +
    geom_point(size = 3, color = df.i$color) +
    geom_label_repel(aes(label = labels), size = 3.5) +
    scale_x_continuous(limits = c(0, 125)) +
    ggtitle(i) + 
    xlab("Regulons Rank") + 
    ylab("Specificity Score") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  plot_list[[i]] <- p
}

# Combine waterfall plots
waterfall_combined <- wrap_plots(plot_list, ncol = 3)
print(waterfall_combined)

# 6. Binary regulon analysis ----

# Calculate activity thresholds
regulonAUC_matrix <- assay(sub_regulonAUC)
bin.T <- AUCell_exploreThresholds(regulonAUC_matrix,
                                 smallestPopPercent = 0.25,
                                 assignCells = TRUE,
                                 plotHist = FALSE,
                                 verbose = FALSE)

# Binary conversion
regulonBin <- lapply(rownames(regulonAUC_matrix), function(reg) {
  as.numeric(colnames(regulonAUC_matrix)) %in% bin.T[[reg]][["assignment"]]
})

regulonBin <- do.call("rbind", regulonBin)
dimnames(regulonBin) <- list(rownames(regulonAUC_matrix), colnames(regulonAUC_matrix))

# 7. Differential analysis with Seurat ----

# Add regulon data as new assays
sce[["Regulon"]] <- CreateAssayObject(counts = regulonAUC_matrix)
sce[["binRegulon"]] <- CreateAssayObject(counts = regulonBin)
DefaultAssay(sce) <- "Regulon"

# Scale data and find markers
sce <- ScaleData(sce, features = rownames(sce))
Idents(sce) <- "celltype"
deg <- FindAllMarkers(sce, only.pos = TRUE, logfc.threshold = 0)
top_regulons <- group_by(deg, cluster) %>% 
  top_n(10, avg_log2FC) %>% 
  pull(gene) %>% 
  unique()

# Create heatmap of top regulons
regulon_heatmap <- DoHeatmap(sce, features = top_regulons, label = FALSE,
                            group.colors = allcolour, group.by = "celltype", angle = 45) +
  scale_fill_gradientn(colors = c("#88558D", "black", "#E9E55A")) + 
  theme(text = element_text(size = 14, family = "arial")) +
  guides(fill = NULL)

print(regulon_heatmap)

# Save key results
save(list = c("sce", "sub_regulonAUC", "rss", "regulonActivity_byGroup_Scaled", 
              "deg", "top_regulons"),
     file = "pyscenic_results.RData")

print("pySCENIC visualization completed successfully!")
print(paste("Number of regulons analyzed:", nrow(sub_regulonAUC)))
print(paste("Number of cells:", ncol(sub_regulonAUC)))
print(paste("Top regulons by specificity:", head(rownames(rss)[order(rowSums(rss), decreasing = TRUE)], 10)))
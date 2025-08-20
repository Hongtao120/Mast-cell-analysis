# Mast cell analysis - Extract cell information for RNA velocity
# This script extracts cell information from Seurat object for RNA velocity analysis
# Author: Based on analysis from README.md
# Date: 2024

# Load required libraries
library(Seurat)

# Set working directory (modify as needed)
# setwd("C:/Users/Taotao/OneDrive/桌面/data/2024.8.14 harmony analysis/RNA velocity")

# Load the processed Mast cell object
# MC <- readRDS("../MC_by_trace.rds")

# Extract diffusion map coordinates (from previous analysis)
# diffusionmap <- tmp[,1:2] # See Diffusion map section
# write.csv(diffusionmap, file = "./diffusionmap.csv")

# Extract cell IDs
df <- data.frame(Cells=Cells(MC))
write.csv(df$Cells, file = "./cellID_obs.csv", row.names = FALSE)

# Extract UMAP coordinates
cell_embeddings <- Embeddings(MC, reduction = "umap")
rownames(cell_embeddings) <- df$Cells
write.csv(cell_embeddings, file = "./cell_embeddings.csv") 

# Extract cell type information
clusters_obs <- MC$celltype
names(clusters_obs) <- df$Cells
write.csv(clusters_obs, file = "./clusters_obs.csv")

print("Cell information extracted successfully!")
print(paste("Number of cells:", nrow(df)))
print(paste("Cell types:", length(unique(clusters_obs))))
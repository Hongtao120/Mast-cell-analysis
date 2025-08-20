# Mast cell analysis - Data preprocessing for pySCENIC
# This script prepares Seurat data for pySCENIC analysis
# Author: Based on analysis from README.md
# Date: 2024

# Load required libraries
library(Seurat)

# Set working directory (modify as needed)
# setwd("path/to/your/analysis")

# Load the Seurat object
# seurat_obj <- readRDS("path/to/your/seurat_object.rds")

# Function to prepare Seurat object for pySCENIC
prepare_seurat_for_pyscenic <- function(seurat_obj) {
  
  # Extract raw count matrix
  expr_matrix <- GetAssayData(seurat_obj, slot = "counts")
  
  # Filter out genes with zero expression across all cells
  genes_to_keep <- rowSums(expr_matrix) > 0
  expr_matrix_filtered <- expr_matrix[genes_to_keep, ]
  
  # Create a new Seurat object with filtered data
  seurat_obj_filtered <- CreateSeuratObject(
    counts = expr_matrix_filtered,
    meta.data = seurat_obj@meta.data
  )
  
  # Copy normalized data if available
  if ("data" %in% slotNames(seurat_obj@assays$RNA)) {
    norm_data_filtered <- GetAssayData(seurat_obj, slot = "data")[genes_to_keep, ]
    seurat_obj_filtered <- SetAssayData(seurat_obj_filtered, slot = "data", new.data = norm_data_filtered)
  }
  
  # Copy scaled data if available
  if (length(seurat_obj@assays$RNA@scale.data) > 0) {
    scaled_data_filtered <- GetAssayData(seurat_obj, slot = "scale.data")[genes_to_keep, ]
    seurat_obj_filtered <- SetAssayData(seurat_obj_filtered, slot = "scale.data", new.data = scaled_data_filtered)
  }
  
  # Copy reduction information
  for (reduction in names(seurat_obj@reductions)) {
    seurat_obj_filtered@reductions[[reduction]] <- seurat_obj@reductions[[reduction]]
  }
  
  return(seurat_obj_filtered)
}

# Prepare the data
seurat_obj_filtered <- prepare_seurat_for_pyscenic(seurat_obj)

# Check the filtered object
print(seurat_obj_filtered)

# Export expression matrix for pySCENIC (cells as rows, genes as columns)
write.csv(t(as.matrix(seurat_obj_filtered@assays$RNA@counts)), file = "sce_exp.csv")

# Export cell information
cellInfo <- seurat_obj_filtered@meta.data[, c("celltype", "nCount_RNA", "nFeature_RNA")]
colnames(cellInfo) <- c('celltype', 'nGene', 'nUMI')
head(cellInfo)
write.csv(cellInfo, file = "cellInfo.csv")

print("Data exported successfully!")
print(paste("Expression matrix dimensions:", nrow(seurat_obj_filtered), "cells x", ncol(seurat_obj_filtered), "genes"))
print(paste("Cell info saved with", nrow(cellInfo), "cells"))
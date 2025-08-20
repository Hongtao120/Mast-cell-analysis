# Mast Cell Analysis - Environment Setup
# This script sets up the required environment and libraries
# Author: Based on analysis from README.md
# Date: 2024

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      print(paste("Installing", pkg))
      if (pkg %in% c("Seurat", "dplyr", "ggplot2", "stringr")) {
        install.packages(pkg)
      } else if (pkg %in% c("SingleCellExperiment", "destiny", "org.Mm.eg.db", "AUCell", "SCENIC")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else if (pkg == "dyno") {
        devtools::install_github("dynverse/dyno")
      } else if (pkg == "irGSEA") {
        devtools::install_github("chuiqin/irGSEA")
      } else {
        install.packages(pkg)
      }
    }
  }
}

# Required R packages ----
r_packages <- c(
  # Core single-cell packages
  "Seurat", "SingleCellExperiment",
  
  # Data manipulation
  "dplyr", "stringr", "magrittr", "data.table",
  
  # Visualization
  "ggplot2", "ggpubr", "ggimage", "patchwork", "plotly",
  "RColorBrewer", "ComplexHeatmap", "circlize", "ggrepel",
  
  # Trajectory inference
  "destiny", "dyno", "slingshot", "monocle", 
  
  # Batch correction
  "harmony",
  
  # TF analysis
  "SCopeLoomR", "AUCell", "SCENIC",
  
  # RNA velocity
  "SeuratDisk",
  
  # Gene set analysis
  "irGSEA", "UCell", "Nebulosa",
  
  # Annotation
  "org.Mm.eg.db",
  
  # Other utilities
  "tidydr", "clustree", "KernSmooth", "BiocParallel", 
  "grid", "scRNAseq", "reshape2", "Biobase", "ggsci"
)

print("Installing/checking R packages...")
install_if_missing(r_packages)

# Load core packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Custom theme function (used in visualizations)
theme_dr <- function() {
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.text = element_text(size = 12)
  )
}

# Python packages required (install separately) ----
python_requirements <- c(
  "scanpy>=1.8.0",
  "scvelo>=0.2.4", 
  "loompy>=3.0.0",
  "pandas>=1.3.0",
  "numpy>=1.21.0",
  "matplotlib>=3.4.0",
  "anndata>=0.8.0"
)

cat("Required Python packages:\n")
cat(paste(python_requirements, collapse = "\n"))
cat("\n\nInstall with: pip install", paste(python_requirements, collapse = " "))

# pySCENIC reference files ----
cat("\n\nRequired pySCENIC reference files:\n")
cat("1. allTFs_mm.txt - Mouse transcription factors list\n")
cat("2. mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather - Motif rankings\n") 
cat("3. motifs-v9-nr.mgi-m0.001-o0.0.tbl - Motif annotations\n")
cat("\nDownload from: https://resources.aertslab.org/cistarget/\n")

# Data requirements ----
cat("\nRequired input data:\n")
cat("1. Seurat object with mast cells (e.g., 3cell.rds)\n")
cat("2. Velocyto loom files for RNA velocity:\n")
cat("   - IE-hpD07_cite.loom\n")
cat("   - IE-hpD14_cite.loom\n") 
cat("   - LP-hpD0_cite.loom\n")
cat("   - LP-hpD07_cite.loom\n")
cat("   - LP-hpD14_cite.loom\n")

print("Environment setup completed!")
print(paste("R version:", R.version.string))
print(paste("Seurat version:", packageVersion("Seurat")))

# Check for common issues
if (!requireNamespace("Seurat", quietly = TRUE)) {
  warning("Seurat not found. Please install with: install.packages('Seurat')")
}

if (!requireNamespace("harmony", quietly = TRUE)) {
  warning("harmony not found. Please install with: install.packages('harmony')")
}
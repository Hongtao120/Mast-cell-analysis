# Basic Analysis

This directory contains the basic clustering and annotation analysis of mast cells.

## Contents

- **src/01_basic_analysis.R** - Main analysis script
- **Output files** - Generated plots and visualizations

## Script Overview

The `01_basic_analysis.R` script performs:

1. **Mast cell extraction** from the full dataset
2. **Cell cycle scoring** using Seurat's built-in markers
3. **Clustering** with Harmony batch correction
4. **Cell type annotation** based on marker genes:
   - Cycling MC
   - Mcpt9 medium MC  
   - Mcpt9 high MC
   - Nr4a1 high MC
   - Lrmda+ MC
5. **Visualization** including UMAP plots by celltype, group, and tissue
6. **Cell proportion analysis** across samples and celltypes

## Usage

```r
# Load your Seurat object first
MC <- readRDS("path/to/your/3cell.rds")

# Run the analysis
source("src/01_basic_analysis.R")
```

## Outputs

Generated visualizations include:
- Dimplot for mast cell subclusters.svg
- Dimplot for Mast cell origins.svg  
- Dimplot for mast cell by group.svg
- FeaturePlot.svg
- Dotplot.svg
- Dotplot by celltype.svg
- Cell proportion plots
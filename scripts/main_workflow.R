# Mast Cell Analysis - Main Workflow
# This script provides an overview of the complete analysis pipeline
# Author: Based on analysis from README.md
# Date: 2024

# This analysis pipeline performs comprehensive analysis of mast cells from cite-seq data
# The pipeline includes clustering, trajectory inference, transcription factor analysis,
# cell-cell communication analysis, and gene set enrichment analysis.

# Prerequisites ----
# 1. Seurat object with mast cells: "./3cell.rds" or similar
# 2. pySCENIC reference files (TF list, motif rankings, annotations)
# 3. Velocyto loom files for RNA velocity analysis

# Analysis Pipeline Overview ----

# Step 1: Basic Analysis
# Script: 01_basic_analysis.R
# Purpose: Clustering, annotation, and basic visualization of mast cells
# Input: Raw Seurat object
# Output: Processed mast cell object with subtypes

# Step 2: Trajectory Inference  
# Scripts: 02_trajectory_inference.R, 02c_diffusion_map_python.py
# Purpose: Infer developmental trajectories using multiple methods
# Methods: Diffusion maps, dyno, Slingshot, Monocle2
# Input: Processed mast cell object
# Output: Trajectory models and visualizations

# Step 3: RNA Velocity Analysis
# Scripts: 02a_extract_cell_info_for_rna_velocity.R, 02b_rna_velocity_analysis.py
# Purpose: Analyze RNA velocity to understand cell dynamics
# Input: Seurat object + Velocyto loom files
# Output: Velocity vectors and visualizations

# Step 4: Transcription Factor Analysis
# Scripts: 03a_prepare_data_for_pyscenic.R, 03b_create_loom_file.py, 
#          03c_run_pyscenic.sh, 03d_pyscenic_visualization.R
# Purpose: Identify active transcription factors and regulatory networks
# Input: Processed mast cell object
# Output: Regulon activities and TF networks

# Step 5: Cell-Cell Communication
# Script: 04_cellchat_analysis.R
# Purpose: Analyze cell-cell communication patterns
# Input: Processed objects (mast cells + epithelial cells)
# Output: Communication networks and Cx3cr1 analysis

# Step 6: Gene Set Enrichment
# Script: 05_gene_set_scoring.R
# Purpose: Functional annotation using multiple gene set databases
# Input: Processed mast cell object
# Output: Pathway enrichment scores and visualizations

# Running the Analysis ----

# 1. Set up environment
print("Setting up analysis environment...")
source("scripts/00_setup_environment.R")

# 2. Run basic analysis
print("Step 1: Running basic analysis...")
source("scripts/01_basic_analysis.R")

# 3. Run trajectory inference
print("Step 2: Running trajectory inference...")
source("scripts/02_trajectory_inference.R")

# Note: Python scripts should be run separately:
# python scripts/02c_diffusion_map_python.py
# python scripts/02b_rna_velocity_analysis.py

# 4. Run TF analysis (requires cluster computing for pySCENIC)
print("Step 4: Preparing for TF analysis...")
source("scripts/03a_prepare_data_for_pyscenic.R")
# Note: Run 03b_create_loom_file.py and 03c_run_pyscenic.sh before visualization
# source("scripts/03d_pyscenic_visualization.R")

# 5. Run CellChat analysis
print("Step 5: Running CellChat analysis...")
source("scripts/04_cellchat_analysis.R")

# 6. Run gene set scoring
print("Step 6: Running gene set scoring...")
source("scripts/05_gene_set_scoring.R")

print("Analysis pipeline completed!")

# Output Summary ----
# The analysis generates the following outputs:
# - Basic analysis: Cell clusters, annotations, UMAP plots
# - Trajectory: Diffusion maps, pseudotime trajectories
# - RNA velocity: Velocity vectors, dynamical analysis
# - TF analysis: Regulon activities, transcription factor networks
# - CellChat: Cell communication patterns, Cx3cr1 expression
# - GSEA: Pathway enrichment scores, functional annotations

# All visualizations are saved in their respective analysis directories:
# - Basic analysis/
# - Diffusionmap/
# - RNA velocity/figures/
# - TF analysis/
# - Cellchat/
# - GSEA/
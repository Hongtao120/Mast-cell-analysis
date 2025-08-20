# TF Analysis

This directory contains transcription factor analysis using pySCENIC.

## Contents

- **src/** - Source code for TF analysis
  - `03a_prepare_data_for_pyscenic.R` - Data preparation
  - `03b_create_loom_file.py` - Convert to loom format
  - `03c_run_pyscenic.sh` - pySCENIC pipeline (for cluster)
  - `03d_pyscenic_visualization.R` - Results visualization
- **Output files** - Generated plots and analysis results

## Pipeline Overview

### 1. Data Preparation (`03a_prepare_data_for_pyscenic.R`)
- Filters genes with zero expression
- Exports expression matrix and cell metadata
- Prepares data in format required by pySCENIC

### 2. Loom File Creation (`03b_create_loom_file.py`)
- Converts CSV expression matrix to loom format
- Required for pySCENIC input

### 3. pySCENIC Analysis (`03c_run_pyscenic.sh`)
- **GRN inference** using GRNBoost2
- **Regulon prediction** using cisTarget
- **Activity scoring** using AUCell
- Designed for cluster computing with SLURM

### 4. Visualization (`03d_pyscenic_visualization.R`)
- Regulon activity heatmaps
- RSS (Regulon Specificity Score) analysis
- Focus on key mast cell TFs: GATA2(+), SOX4(+)
- Waterfall plots for top regulons

## Usage

```r
# 1. Prepare data
source("src/03a_prepare_data_for_pyscenic.R")

# 2. Create loom file
# python src/03b_create_loom_file.py

# 3. Run pySCENIC (submit to cluster)
# sbatch src/03c_run_pyscenic.sh

# 4. Visualize results
source("src/03d_pyscenic_visualization.R")
```

## Requirements

### pySCENIC Reference Files
Download from: https://resources.aertslab.org/cistarget/
- `allTFs_mm.txt`
- `mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `motifs-v9-nr.mgi-m0.001-o0.0.tbl`

### Computational Resources
pySCENIC analysis is computationally intensive and benefits from:
- Multiple CPU cores (10+ recommended)
- High memory (32GB+ recommended)
- Cluster computing environment

## Outputs

Key transcription factors analyzed:
- **GATA2(+)** - Key mast cell regulator
- **SOX4(+)** - Mast cell development factor

Generated visualizations:
- Heatmap for MC TF by celltype.svg
- Heatmap for MC TF by group.svg
- RSS plots showing TF specificity
- Top regulon analysis
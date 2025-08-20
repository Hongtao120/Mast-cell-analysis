# Mast Cell Analysis Pipeline

This repository contains a comprehensive analysis pipeline for mast cell analysis based on total cite-seq dataset. The analysis has been reorganized from a single large README into modular, reusable scripts.

## 📁 Repository Structure

```
Mast-cell-analysis/
├── scripts/                          # Main analysis scripts
│   ├── 00_setup_environment.R        # Environment setup and dependencies
│   ├── 01_basic_analysis.R           # Clustering, annotation, visualization
│   ├── 02_trajectory_inference.R     # Trajectory analysis (R methods)
│   ├── 02a_extract_cell_info_for_rna_velocity.R  # RNA velocity preparation
│   ├── 02b_rna_velocity_analysis.py  # RNA velocity analysis (Python)
│   ├── 02c_diffusion_map_python.py  # Diffusion map (Python)
│   ├── 03a_prepare_data_for_pyscenic.R  # pySCENIC data preparation
│   ├── 03b_create_loom_file.py       # Convert to loom format
│   ├── 03c_run_pyscenic.sh           # pySCENIC pipeline (cluster script)
│   ├── 03d_pyscenic_visualization.R  # pySCENIC results visualization
│   ├── 04_cellchat_analysis.R        # Cell-cell communication analysis
│   ├── 05_gene_set_scoring.R         # Gene set enrichment analysis
│   └── main_workflow.R               # Main workflow overview
├── Basic analysis/                    # Basic analysis outputs
│   └── src/                          # Source code (future use)
├── Diffusionmap/                     # Diffusion map outputs
│   └── src/                          # Source code (future use)
├── TF analysis/                      # Transcription factor analysis outputs
│   └── src/                          # Source code (future use)
├── RNA velocity/                     # RNA velocity analysis outputs
│   ├── figures/                      # Velocity plots
│   └── src/                          # Source code (future use)
├── Cellchat/                         # Cell communication analysis outputs
│   └── src/                          # Source code (future use)
├── GSEA/                             # Gene set enrichment outputs
│   └── src/                          # Source code (future use)
├── Projection/                       # Projection analysis outputs
└── README.md                         # This file
```

## 🚀 Quick Start

### 1. Environment Setup

First, set up your R environment with required packages:

```r
source("scripts/00_setup_environment.R")
```

### 2. Data Requirements

Ensure you have the following input files:

**R Data:**
- Seurat object with mast cells (e.g., `3cell.rds`)
- Total cite-seq dataset (e.g., `totalCITE-seq.rds`)

**RNA Velocity (Optional):**
- Velocyto loom files for each sample:
  - `IE-hpD07_cite.loom`
  - `IE-hpD14_cite.loom`
  - `LP-hpD0_cite.loom`
  - `LP-hpD07_cite.loom`
  - `LP-hpD14_cite.loom`

**pySCENIC Reference Files (Optional):**
- `allTFs_mm.txt` - Mouse transcription factors
- `mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `motifs-v9-nr.mgi-m0.001-o0.0.tbl`

Download pySCENIC files from: https://resources.aertslab.org/cistarget/

### 3. Run Analysis

You can run individual analysis modules or follow the complete workflow:

```r
# Load your data
MC <- readRDS("path/to/your/3cell.rds")

# Run basic analysis
source("scripts/01_basic_analysis.R")

# Run trajectory inference
source("scripts/02_trajectory_inference.R")

# For Python scripts, run separately:
# python scripts/02b_rna_velocity_analysis.py
# python scripts/02c_diffusion_map_python.py

# Run gene set scoring
source("scripts/05_gene_set_scoring.R")
```

## 📊 Analysis Modules

### 1. Basic Analysis (`01_basic_analysis.R`)
- **Purpose:** Clustering, annotation, and visualization of mast cells
- **Methods:** Seurat workflow with Harmony batch correction
- **Outputs:** 
  - Cell clusters and annotations
  - UMAP plots by celltype, group, tissue
  - Cell proportion analysis
  - Marker gene identification

### 2. Trajectory Inference (`02_trajectory_inference.R`)
- **Purpose:** Infer developmental trajectories
- **Methods:** 
  - Diffusion maps (R & Python)
  - dyno (multiple TI methods)
  - Slingshot
  - Monocle2
- **Outputs:** Trajectory plots and pseudotime analysis

### 3. RNA Velocity Analysis (`02a_*.R`, `02b_*.py`)
- **Purpose:** Analyze RNA velocity dynamics
- **Methods:** scVelo with dynamical modeling
- **Outputs:** 
  - Velocity vector fields
  - Cell cycle analysis
  - Gene expression dynamics

### 4. Transcription Factor Analysis (`03a-d_pyscenic_*`)
- **Purpose:** Identify active regulatory networks
- **Methods:** pySCENIC workflow
- **Outputs:**
  - Regulon activities
  - TF specificity scores (RSS)
  - Regulatory network visualization

### 5. Cell Communication (`04_cellchat_analysis.R`)
- **Purpose:** Analyze cell-cell communication
- **Methods:** Cx3cr1 expression analysis
- **Outputs:** Communication network plots

### 6. Gene Set Scoring (`05_gene_set_scoring.R`)
- **Purpose:** Functional pathway analysis
- **Methods:** irGSEA with multiple scoring methods
- **Outputs:**
  - Hallmark pathway scores
  - KEGG pathway enrichment
  - GO-BP functional analysis

## 🔧 Dependencies

### R Packages
Core packages automatically installed by setup script:
- **Single-cell:** Seurat, SingleCellExperiment
- **Visualization:** ggplot2, ComplexHeatmap, patchwork
- **Trajectory:** dyno, slingshot, monocle, destiny
- **Batch correction:** harmony
- **TF analysis:** SCENIC, AUCell, SCopeLoomR
- **Gene sets:** irGSEA, UCell
- **Annotation:** org.Mm.eg.db

### Python Packages
Install separately with pip:
```bash
pip install scanpy scvelo loompy pandas numpy matplotlib anndata
```

## 💡 Key Features

### ✅ Modular Design
- Each analysis type in separate, focused scripts
- Easy to run individual modules
- Clear input/output relationships

### ✅ Comprehensive Documentation
- Detailed comments in all scripts
- Function documentation
- Clear variable naming

### ✅ Reproducible Analysis
- Version control friendly
- Consistent coding style
- Saved intermediate results

### ✅ Flexible Workflow
- Run complete pipeline or individual modules
- Easy parameter modification
- Multiple visualization options

## 🎯 Output Summary

The analysis generates comprehensive outputs including:

- **Clustering:** Cell type annotations and cluster markers
- **Trajectories:** Developmental pathways and pseudotime
- **Regulation:** Active transcription factors and networks
- **Communication:** Cell-cell interaction patterns
- **Function:** Pathway enrichment and biological processes
- **Visualization:** High-quality plots and interactive graphics

## 📝 Citation

If you use this analysis pipeline, please cite the original methods:

- **Seurat:** Hao et al. (2021) Nature Biotechnology
- **Harmony:** Korsunsky et al. (2019) Nature Methods  
- **scVelo:** Bergen et al. (2020) Nature Biotechnology
- **pySCENIC:** Aibar et al. (2017) Nature Methods
- **irGSEA:** Qin et al. (2022) bioRxiv

## 🤝 Contributing

To contribute to this analysis pipeline:
1. Fork the repository
2. Create analysis modules following the existing structure
3. Add comprehensive documentation
4. Submit pull request with clear description

## 📞 Support

For questions about this analysis pipeline:
- Check script comments for detailed explanations
- Review the original methods papers
- Create GitHub issues for specific problems

---

**Note:** This organized structure replaces the previous single README.md file containing all code. The original analysis logic is preserved but now distributed across focused, reusable scripts.
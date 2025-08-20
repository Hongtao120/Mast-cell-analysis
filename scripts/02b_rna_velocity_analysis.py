# Mast cell analysis - RNA velocity analysis
# This script performs RNA velocity analysis using scVelo
# Author: Based on analysis from README.md
# Date: 2024

import scvelo as scv
import numpy as np
import scanpy as sc
import pandas as pd
import anndata
import matplotlib as plt
import re

# Set verbosity level
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')

def transform_string_IED7(s):
    """Transform cell ID for IED7 sample"""
    match = re.match(r"HP7day-IEMC:(.*)x", s)
    if match:
        extracted_part = match.group(1)
        transformed = f"IE-hpD07_cite_{extracted_part}"
        return transformed
    return s

def transform_string_IED14(s):
    """Transform cell ID for IED14 sample"""
    match = re.match(r"HP14day-IEMC:(.*)x", s)
    if match:
        extracted_part = match.group(1)
        transformed = f"IE-hpD14_cite_{extracted_part}"
        return transformed
    return s

def transform_string_LPD0(s):
    """Transform cell ID for LPD0 sample"""
    match = re.match(r"HP0day-LPMC:(.*)x", s)
    if match:
        extracted_part = match.group(1)
        transformed = f"LP-hpD0_cite_{extracted_part}"
        return transformed
    return s

def transform_string_LPD7(s):
    """Transform cell ID for LPD7 sample"""
    match = re.match(r"HP7day-LPMC:(.*)x", s)
    if match:
        extracted_part = match.group(1)
        transformed = f"LP-hpD07_cite_{extracted_part}"
        return transformed
    return s

def transform_string_LPD14(s):
    """Transform cell ID for LPD14 sample"""
    match = re.match(r"HP14day-LPMC:(.*)x", s)
    if match:
        extracted_part = match.group(1)
        transformed = f"LP-hpD14_cite_{extracted_part}"
        return transformed
    return s

# Load loom files
print("Loading loom files...")
sample_IED7 = anndata.read_loom("./IE-hpD07_cite.loom")
sample_IED14 = anndata.read_loom("./IE-hpD14_cite.loom")
sample_LPD0 = anndata.read_loom("./LP-hpD0_cite.loom")
sample_LPD07 = anndata.read_loom("./LP-hpD07_cite.loom")
sample_LPD14 = anndata.read_loom("./LP-hpD14_cite.loom")

# Transform cell IDs for each sample
print("Transforming cell IDs...")
sample_IED7.obs['obs_names'] = sample_IED7.obs['obs_names'].apply(transform_string_IED7)
sample_IED14.obs['obs_names'] = sample_IED14.obs['obs_names'].apply(transform_string_IED14)
sample_LPD0.obs['obs_names'] = sample_LPD0.obs['obs_names'].apply(transform_string_LPD0)
sample_LPD07.obs['obs_names'] = sample_LPD07.obs['obs_names'].apply(transform_string_LPD7)
sample_LPD14.obs['obs_names'] = sample_LPD14.obs['obs_names'].apply(transform_string_LPD14)

# Load cell information from R
print("Loading cell information from R...")
sample_obs = pd.read_csv("./cellID_obs.csv")
umap = pd.read_csv("./cell_embeddings.csv")
cell_clusters = pd.read_csv("./clusters_obs.csv")
dm = pd.read_csv("./diffusionmap.csv")

# Filter cells in each sample
print("Filtering cells...")
sample_IED7 = sample_IED7[np.isin(sample_IED7.obs['obs_names'], sample_obs["x"])]
sample_IED14 = sample_IED14[np.isin(sample_IED14.obs['obs_names'], sample_obs["x"])]
sample_LPD0 = sample_LPD0[np.isin(sample_LPD0.obs['obs_names'], sample_obs["x"])]
sample_LPD07 = sample_LPD07[np.isin(sample_LPD07.obs['obs_names'], sample_obs["x"])]
sample_LPD14 = sample_LPD14[np.isin(sample_LPD14.obs['obs_names'], sample_obs["x"])]

# Concatenate all samples
print("Concatenating samples...")
adata = sample_IED7.concatenate(sample_IED14, sample_LPD0, sample_LPD07, sample_LPD14)

# Add UMAP coordinates
print("Adding UMAP coordinates...")
adata_index = pd.DataFrame(adata.obs['obs_names'])
adata_index = adata_index.rename(columns={'obs_names': 'Cell ID'})

umap = umap.rename(columns={'Unnamed: 0': 'Cell ID'})
umap = umap[np.isin(umap["Cell ID"], adata_index["Cell ID"])]
umap_ordered = adata_index.merge(umap, on="Cell ID")
umap_ordered = umap_ordered.iloc[:, 1:]
adata.obsm['X_umap'] = umap_ordered.values

# Add diffusion map coordinates
print("Adding diffusion map coordinates...")
dm = dm.rename(columns={'Unnamed: 0': 'Cell ID'})
dm = dm[np.isin(dm["Cell ID"], adata_index["Cell ID"])]
dm_ordered = adata_index.merge(dm, on="Cell ID")
dm_ordered = dm_ordered.iloc[:, 1:]
adata.obsm['dm'] = dm_ordered.values

# Add cell type information
print("Adding cell type information...")
cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'Cell ID'})
cell_clusters_ordered = adata_index.merge(cell_clusters, on="Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:, 1:]
adata.obs['celltype'] = cell_clusters_ordered.values

# Preprocessing
print("Preprocessing data...")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Calculate RNA velocity
print("Calculating RNA velocity...")
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata)

# Save the processed data
print("Saving processed data...")
adata.write('data.h5ad')

# Define color palette for visualization
palette = {
    "Mcpt9 high MC": "#D0AFC4",
    "Mcpt9 medium MC": "#89558D",
    "Cycling MC": "#AFC2D9", 
    "Lrmda+ MC": "#435B95",
    "Nr4a1 high MC": "#79B99D", 
    "IE-hpD07": "#D55640",
    "IE-hpD14": "#E69F84",
    "LP-hpD0": "#6CB8D2", 
    "LP-hpD07": "#479D88",
    "LP-hpD14": "#415284"
}

# Generate visualizations
print("Generating visualizations...")

# Grid plots
scv.pl.velocity_embedding_grid(adata, basis='umap', color='celltype', 
                               save='embedding_grid.svg', title='RNA Velocity', 
                               scale=0.25, palette=palette, figsize=(4.5, 4))

scv.pl.velocity_embedding_grid(adata, basis='dm', color='celltype', 
                               save='embedding_grid_dm.svg', title='RNA Velocity', 
                               scale=0.25, palette=palette, figsize=(4.5, 4))

# Stream plots
scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype', 
                                 save='embedding_stream.svg', title='RNA Velocity', 
                                 palette=palette, figsize=(4.5, 4))

scv.pl.velocity_embedding_stream(adata, basis='dm', color='celltype', 
                                 save='embedding_stream_dm.svg', title='RNA Velocity', 
                                 palette=palette, figsize=(4.5, 4), legend_loc='right')

# Arrow plot
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120, 
                          color='celltype', save='arrow_stream.svg', 
                          title='RNA Velocity', palette=palette, figsize=(4.5, 4))

# Set gene names as index
adata.var.index = adata.var['var_names']

# Gene velocity plots
scv.pl.velocity(adata, ['Top2a', 'Lrmda', 'Nr4a1', 'Mcpt9'], ncols=2, 
                color='celltype', save='gene.svg', figsize=(4.5, 4))
scv.pl.velocity(adata, ['Top2a', 'Lrmda', 'Nr4a1', 'Mcpt9'], ncols=2, 
                color='celltype', save='gene_dm.svg', basis='dm', figsize=(4.5, 4))

# Cell cycle analysis
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, 
               perc=[5, 95], save='cycling progenitors.svg', figsize=(4.5, 4))
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, 
               perc=[5, 95], basis="dm", save='cycling progenitors_dm.svg', figsize=(4.5, 4))

# Cell cycle markers
s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
s_genes = scv.get_df(adata[:, s_genes], 'spearmans_score', sort_values=True).index
g2m_genes = scv.get_df(adata[:, g2m_genes], 'spearmans_score', sort_values=True).index

kwargs = dict(frameon=False, ylabel='cell cycle genes')
scv.pl.scatter(adata, list(s_genes[:2]) + list(g2m_genes[:3]), 
               save='cycling marker.svg', color='celltype', **kwargs)

# Specific gene plots
scv.pl.velocity(adata, ['Hells', 'Top2a'], ncols=2, add_outline=True, 
                color='celltype', save='hells and top2a.svg', figsize=(4.5, 4))
scv.pl.velocity(adata, ['Hells', 'Top2a'], ncols=2, add_outline=True, 
                color='celltype', basis='dm', save='hells and top2a_dm.svg', figsize=(4.5, 4))

print("RNA velocity analysis completed successfully!")
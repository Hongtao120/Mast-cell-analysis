# Mast cell analysis - Diffusion Map in Python  
# This script performs diffusion map analysis using scanpy
# Author: Based on analysis from README.md
# Date: 2024

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set scanpy settings
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

def run_diffusion_map_analysis(h5ad_file):
    """
    Run diffusion map analysis on h5ad file
    
    Parameters:
    h5ad_file: Path to h5ad file (converted from Seurat)
    """
    
    print(f"Loading data from {h5ad_file}...")
    adata = sc.read(h5ad_file)
    
    print(f"Data shape: {adata.shape}")
    print(f"Available obsm keys: {list(adata.obsm.keys())}")
    print(f"Available obs columns: {list(adata.obs.columns)}")
    
    # Build neighborhood graph using harmony coordinates if available
    if 'X_harmony' in adata.obsm.keys():
        print("Using Harmony coordinates for neighbor calculation...")
        sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50, use_rep='X_harmony', method='gauss')
    else:
        print("Using PCA coordinates for neighbor calculation...")
        sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)
    
    # Calculate diffusion map
    print("Calculating diffusion map...")
    sc.tl.diffmap(adata)
    
    # Store diffusion map coordinates
    adata.obsm['X_diffmap_'] = adata.obsm['X_diffmap'][:, 1:]
    
    # Generate visualizations
    print("Generating visualizations...")
    
    # Plot by celltype if available
    if 'celltype' in adata.obs.columns:
        sc.pl.embedding(adata, 'diffmap', color='celltype', save='_diffmap_celltype.svg')
        sc.pl.embedding(adata, 'diffmap', color='celltype', save='_diffmap_celltype.pdf')
    
    # Plot by group if available  
    if 'group' in adata.obs.columns:
        sc.pl.embedding(adata, 'diffmap', color='group', save='_diffmap_group.svg')
        sc.pl.embedding(adata, 'diffmap', color='group', save='_diffmap_group.pdf')
    
    # Plot by tissue if available
    if 'tissue' in adata.obs.columns:
        sc.pl.embedding(adata, 'diffmap', color='tissue', save='_diffmap_tissue.svg')
        sc.pl.embedding(adata, 'diffmap', color='tissue', save='_diffmap_tissue.pdf')
    
    # Save results
    adata.write('diffmap_results.h5ad')
    
    # Export diffusion map coordinates for R
    diffmap_coords = pd.DataFrame(adata.obsm['X_diffmap'][:, :2], 
                                 columns=['DC1', 'DC2'],
                                 index=adata.obs.index)
    diffmap_coords.to_csv('diffusion_map_coordinates.csv')
    
    print("Diffusion map analysis completed!")
    print(f"Diffusion map shape: {adata.obsm['X_diffmap'].shape}")
    
    return adata

def plot_diffusion_components(adata, n_components=4):
    """
    Plot multiple diffusion components
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    for i in range(min(n_components, 4)):
        if i < 3:  # First 3 components
            dc_x = adata.obsm['X_diffmap'][:, 0]
            dc_y = adata.obsm['X_diffmap'][:, i+1]
            axes[i].scatter(dc_x, dc_y, c=adata.obs['celltype'], s=20, alpha=0.7)
            axes[i].set_xlabel(f'DC1')
            axes[i].set_ylabel(f'DC{i+2}')
            axes[i].set_title(f'Diffusion Components: DC1 vs DC{i+2}')
        else:
            # DC2 vs DC3
            dc_x = adata.obsm['X_diffmap'][:, 1]
            dc_y = adata.obsm['X_diffmap'][:, 2]
            axes[i].scatter(dc_x, dc_y, c=adata.obs['celltype'], s=20, alpha=0.7)
            axes[i].set_xlabel('DC2')
            axes[i].set_ylabel('DC3')
            axes[i].set_title('Diffusion Components: DC2 vs DC3')
    
    plt.tight_layout()
    plt.savefig('diffusion_components_overview.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    # Main analysis
    input_file = "test.h5ad"  # Created from Seurat object
    
    try:
        # Run diffusion map analysis
        adata = run_diffusion_map_analysis(input_file)
        
        # Plot components overview
        plot_diffusion_components(adata)
        
        print("Analysis completed successfully!")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        print("Please ensure the h5ad file exists and was properly converted from Seurat.")
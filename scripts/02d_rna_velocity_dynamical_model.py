# RNA Velocity - Dynamical Model Analysis
# This script performs advanced RNA velocity analysis using the dynamical model
# Author: Based on analysis from README.md
# Date: 2024

import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set scvelo settings
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')

def run_dynamical_model(adata_file='data.h5ad'):
    """
    Run dynamical model analysis on preprocessed RNA velocity data
    
    Parameters:
    adata_file: Path to preprocessed h5ad file from basic velocity analysis
    """
    
    print(f"Loading preprocessed data from {adata_file}...")
    adata = scv.read(adata_file)
    
    print(f"Data shape: {adata.shape}")
    print("Running dynamical model analysis...")
    
    # Define color palette
    palette = {
        "Mcpt9 high MC": "#D0AFC4",
        "Mcpt9 medium MC": "#89558D",
        "Cycling MC": "#AFC2D9", 
        "Lrmda+ MC": "#435B95",
        "Nr4a1 high MC": "#79B99D"
    }
    
    # 1. Recover dynamics and run dynamical model ----
    print("Recovering dynamics...")
    scv.tl.recover_dynamics(adata)
    
    print("Computing dynamical velocity...")
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)
    
    # Set gene names as index
    adata.var.index = adata.var['var_names']
    
    # Save results (dynamical model can take a while)
    print("Saving dynamical model results...")
    adata.write('adata_dynamical.h5ad')
    
    # 2. Visualize dynamical velocity ----
    print("Generating dynamical velocity visualizations...")
    
    # Stream plots with dynamical model
    scv.pl.velocity_embedding_stream(adata, basis='umap', 
                                     save='dynamicalstream.svg', 
                                     title='RNA Velocity Dynamical Model', 
                                     color='celltype', palette=palette, 
                                     figsize=(4.5, 4))
    
    scv.pl.velocity_embedding_stream(adata, basis='dm', 
                                     save='dynamicalstream_dm.svg', 
                                     title='RNA Velocity Dynamical Model', 
                                     color='celltype', palette=palette, 
                                     figsize=(4.5, 4))
    
    # 3. Analyze kinetic rate parameters ----
    print("Analyzing kinetic parameters...")
    
    # Filter for high-confidence velocity genes
    df = adata.var
    df = df[(df['fit_likelihood'] > 0.1) & (df['velocity_genes'] == True)]
    
    # Plot kinetic parameters
    kwargs = dict(xscale='log', fontsize=16)
    with scv.GridSpec(ncols=3) as pl:
        pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
        pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', 
                xticks=[.1, .4, 1], **kwargs)
        pl.hist(df['fit_gamma'], xlabel='degradation rate', 
                xticks=[.1, .4, 1], **kwargs)
    
    # Display kinetic parameters summary
    kinetic_summary = scv.get_df(adata, 'fit*', dropna=True).head()
    print("Kinetic parameters summary:")
    print(kinetic_summary)
    
    # 4. Latent time analysis ----
    print("Computing latent time...")
    
    scv.tl.latent_time(adata)
    scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', 
                   size=80, save='Latent time.svg', figsize=(4.5, 4))
    
    # Heatmap of top genes by latent time
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', 
                   col_color='celltype', n_convolve=100, save='heatmap_top genes.pdf')
    
    # 5. Top-likelihood genes analysis ----
    print("Analyzing top-likelihood genes...")
    
    # Plot top 15 genes by likelihood
    top_genes_all = adata.var['fit_likelihood'].sort_values(ascending=False).index
    scv.pl.scatter(adata, basis=top_genes_all[:15], ncols=5, frameon=False, 
                   save='Top-likelihood genes.svg', color='celltype', 
                   figsize=(4.5, 4))
    
    # Specific genes of interest
    var_names = ['Cd63', 'Itgae', 'Mcpt9', 'Cma2']
    if all(gene in adata.var.index for gene in var_names):
        scv.pl.scatter(adata, var_names, frameon=False, 
                       save='Top-likelihood genes choosen.svg', color='celltype')
        scv.pl.scatter(adata, x='latent_time', y=var_names, frameon=False, 
                       save='Top-likelihood genes latent_time.svg', color='celltype')
    
    # 6. Cluster-specific top genes ----
    print("Identifying cluster-specific dynamical genes...")
    
    scv.tl.rank_dynamical_genes(adata, groupby='celltype')
    df_ranks = scv.get_df(adata, 'rank_dynamical_genes/names')
    
    print("Top dynamical genes by cluster:")
    print(df_ranks.head(5))
    
    # Plot top genes for each cluster
    clusters = ['Cycling MC', 'Lrmda+ MC', 'Mcpt9 high MC', 'Mcpt9 medium MC', 'Nr4a1 high MC']
    
    for cluster in clusters:
        if cluster in df_ranks.columns:
            cluster_genes = df_ranks[cluster][:5].dropna()
            if len(cluster_genes) > 0:
                scv.pl.scatter(adata, cluster_genes, ylabel=cluster, frameon=False, 
                               save=f'_{cluster.replace(" ", "_")}_genes.svg', figsize=(4.5, 4))
    
    # 7. Advanced velocity analysis ----
    print("Performing advanced velocity analysis...")
    
    # Speed and coherence
    scv.tl.velocity_confidence(adata)
    scv.pl.scatter(adata, c=['velocity_length', 'velocity_confidence'], 
                   cmap='coolwarm', perc=[5, 95], save='speed and coherence.svg', 
                   figsize=(4.5, 4))
    
    # Velocity graph
    scv.pl.velocity_graph(adata, threshold=.1, save='velocitygraph.svg', 
                          figsize=(4.5, 4))
    
    # PAGA velocity graph
    if 'paga' in adata.uns:
        scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
                    min_edge_width=2, node_size_scale=1.5, save='paga.svg')
    
    # Pseudotime
    if 'dpt_pseudotime' in adata.obs:
        scv.pl.scatter(adata, color='dpt_pseudotime', save='pseudotime.svg',
                       figsize=(4.5, 4))
    
    print("Dynamical model analysis completed!")
    
    # Summary statistics
    n_velocity_genes = sum(adata.var['velocity_genes'])
    mean_likelihood = adata.var['fit_likelihood'].mean()
    
    print(f"Number of velocity genes: {n_velocity_genes}")
    print(f"Mean fit likelihood: {mean_likelihood:.3f}")
    print(f"Latent time range: {adata.obs['latent_time'].min():.2f} - {adata.obs['latent_time'].max():.2f}")
    
    return adata

def analyze_gene_dynamics(adata, genes):
    """
    Analyze dynamics of specific genes
    
    Parameters:
    adata: AnnData object with dynamical model results
    genes: List of gene names to analyze
    """
    
    print(f"Analyzing dynamics for genes: {genes}")
    
    for gene in genes:
        if gene in adata.var.index:
            # Phase portrait
            scv.pl.velocity(adata, gene, save=f'_{gene}_dynamics.svg', figsize=(4.5, 4))
            
            # Expression over latent time
            scv.pl.scatter(adata, x='latent_time', y=gene, save=f'_{gene}_latent_time.svg',
                           color='celltype', figsize=(4.5, 4))
        else:
            print(f"Gene {gene} not found in data")

if __name__ == "__main__":
    # Run dynamical model analysis
    try:
        adata = run_dynamical_model()
        
        # Analyze specific genes
        mast_cell_genes = ['Mcpt9', 'Gata2', 'Sox4', 'Lrmda', 'Nr4a1']
        analyze_gene_dynamics(adata, mast_cell_genes)
        
        print("Dynamical analysis completed successfully!")
        
    except Exception as e:
        print(f"Error during dynamical analysis: {e}")
        print("Please ensure the basic RNA velocity analysis has been completed first.")
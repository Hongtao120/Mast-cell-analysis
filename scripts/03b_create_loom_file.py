# Convert CSV to loom file for pySCENIC analysis
# Author: Based on analysis from README.md
# Date: 2024

import os
import sys
import loompy as lp
import numpy as np
import scanpy as sc

def create_loom_from_csv(csv_file, output_loom):
    """
    Convert CSV expression matrix to loom format for pySCENIC
    
    Parameters:
    csv_file: Path to CSV file with expression data (cells as rows, genes as columns)
    output_loom: Path for output loom file
    """
    
    print(f"Reading CSV file: {csv_file}")
    # Read the CSV file using scanpy
    adata = sc.read_csv(csv_file)
    
    print(f"Data shape: {adata.shape}")
    print(f"Genes: {adata.n_vars}, Cells: {adata.n_obs}")
    
    # Prepare attributes for loom file
    row_attrs = {"Gene": np.array(adata.var_names)}
    col_attrs = {"CellID": np.array(adata.obs_names)}
    
    print(f"Creating loom file: {output_loom}")
    # Create loom file (genes as rows, cells as columns - transpose the matrix)
    lp.create(output_loom, adata.X.transpose(), row_attrs, col_attrs)
    
    print("Loom file created successfully!")
    return output_loom

if __name__ == "__main__":
    # Set working directory and check files
    print(f"Current working directory: {os.getcwd()}")
    print(f"Files in directory: {os.listdir(os.getcwd())}")
    
    # Create loom file from CSV
    csv_file = "sce_exp.csv"
    output_loom = "sce.loom"
    
    if os.path.exists(csv_file):
        create_loom_from_csv(csv_file, output_loom)
    else:
        print(f"Error: {csv_file} not found. Please run the R preprocessing script first.")
        sys.exit(1)
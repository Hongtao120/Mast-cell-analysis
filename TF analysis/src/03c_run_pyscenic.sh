#!/bin/bash
# pySCENIC analysis pipeline
# This script runs the complete pySCENIC workflow for transcription factor analysis
# Author: Based on analysis from README.md
# Date: 2024

#SBATCH -J pyscenic_analysis
#SBATCH -p Acluster
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -N 1

# Set paths (modify these according to your system)
LOOM_FILE="sce.loom"
TF_LIST="allTFs_mm.txt"
MOTIF_RANKINGS="mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
MOTIF_ANNOTATIONS="motifs-v9-nr.mgi-m0.001-o0.0.tbl"

# Output files
GRN_OUTPUT="adj.sample.tsv"
CTX_OUTPUT="reg.csv"
FINAL_OUTPUT="sample_SCENIC.loom"

echo "Starting pySCENIC analysis pipeline..."
echo "Input loom file: $LOOM_FILE"

# Check if input files exist
if [ ! -f "$LOOM_FILE" ]; then
    echo "Error: $LOOM_FILE not found. Please run the preprocessing scripts first."
    exit 1
fi

if [ ! -f "$TF_LIST" ]; then
    echo "Error: $TF_LIST not found. Please download the transcription factor list."
    exit 1
fi

if [ ! -f "$MOTIF_RANKINGS" ]; then
    echo "Error: $MOTIF_RANKINGS not found. Please download the motif rankings database."
    exit 1
fi

if [ ! -f "$MOTIF_ANNOTATIONS" ]; then
    echo "Error: $MOTIF_ANNOTATIONS not found. Please download the motif annotations."
    exit 1
fi

# Step 1: Gene Regulatory Network (GRN) inference using GRNBoost2
echo "Step 1: Running GRN inference with GRNBoost2..."
pyscenic grn \
    --num_workers 10 \
    --output $GRN_OUTPUT \
    --method grnboost2 \
    $LOOM_FILE \
    $TF_LIST

if [ $? -ne 0 ]; then
    echo "Error: GRN inference failed"
    exit 1
fi

echo "GRN inference completed. Output: $GRN_OUTPUT"

# Step 2: Regulon prediction using cisTarget
echo "Step 2: Running regulon prediction with cisTarget..."
pyscenic ctx \
    $GRN_OUTPUT \
    $MOTIF_RANKINGS \
    --annotations_fname $MOTIF_ANNOTATIONS \
    --expression_mtx_fname $LOOM_FILE \
    --mode "dask_multiprocessing" \
    --output $CTX_OUTPUT \
    --num_workers 3 \
    --mask_dropouts

if [ $? -ne 0 ]; then
    echo "Error: Regulon prediction failed"
    exit 1
fi

echo "Regulon prediction completed. Output: $CTX_OUTPUT"

# Step 3: Regulon activity scoring using AUCell
echo "Step 3: Running regulon activity scoring with AUCell..."
pyscenic aucell \
    $LOOM_FILE \
    $CTX_OUTPUT \
    --output $FINAL_OUTPUT \
    --num_workers 3

if [ $? -ne 0 ]; then
    echo "Error: AUCell scoring failed"
    exit 1
fi

echo "pySCENIC analysis completed successfully!"
echo "Final output: $FINAL_OUTPUT"

# Print summary
echo "================================"
echo "pySCENIC Analysis Summary:"
echo "GRN inference output: $GRN_OUTPUT"
echo "Regulon prediction output: $CTX_OUTPUT"
echo "Final regulon activities: $FINAL_OUTPUT"
echo "================================"
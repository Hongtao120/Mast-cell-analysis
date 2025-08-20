# Mast cell analysis - Gene Set Scoring
# This script performs gene set enrichment analysis using irGSEA
# Author: Based on analysis from README.md
# Date: 2024

# Load required libraries
library(UCell)
library(irGSEA)
library(Seurat)
library(Nebulosa)
library(dplyr)
library(org.Mm.eg.db)

# Set working directory (modify as needed)
# setwd("path/to/your/analysis")

# Load Mast cell object
# MC <- readRDS("./MC_by_trace.rds")

# Define color palettes
allcolour_celltype <- c(
  "#D0AFC4", "#89558D", "#AFC2D9", "#435B95", "#79B99D"
)
names(allcolour_celltype) <- c('Mcpt9 high MC', 'Mcpt9 medium MC', 'Cycling MC', 'Nr4a1 high MC', 'Lrmda+ MC')

allcolour_group <- c(
  "#D55640", "#E69F84", "#6CB8D2", "#479D88", "#415284"
)
names(allcolour_group) <- c("IE-hpD07", "IE-hpD14", "LP-hpD0", "LP-hpD07", "LP-hpD14")

# 1. Gene set scoring with irGSEA ----

# Set identities
Idents(MC) <- "celltype"

# 1.1 Hallmark gene set scoring ----
print("Running Hallmark gene set scoring...")
MC <- irGSEA.score(object = MC, assay = "RNA",
                   slot = "data", seeds = 123, ncores = 4,
                   min.cells = 3, min.feature = 0,
                   custom = FALSE, geneset = NULL, msigdb = TRUE,
                   species = "Mus musculus", category = "H",
                   subcategory = NULL, geneid = "symbol",
                   method = c("AUCell", "UCell", "singscore", "ssgsea"),
                   aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                   kcdf = 'Gaussian')

# 1.2 KEGG gene set scoring ----
print("Running KEGG gene set scoring...")
MC <- irGSEA.score(object = MC, assay = "RNA",
                   slot = "data", seeds = 123, ncores = 4,
                   min.cells = 3, min.feature = 0,
                   custom = FALSE, geneset = NULL, msigdb = TRUE,
                   species = "Mus musculus", category = "C2",
                   subcategory = "CP:KEGG", geneid = "symbol",
                   method = c("AUCell", "UCell", "singscore", "ssgsea"),
                   aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                   kcdf = 'Gaussian')

# 1.3 GO-BP gene set scoring ----
print("Running GO-BP gene set scoring...")
MC <- irGSEA.score(object = MC, assay = "RNA",
                   slot = "data", seeds = 123, ncores = 4,
                   min.cells = 3, min.feature = 0,
                   custom = FALSE, geneset = NULL, msigdb = TRUE,
                   species = "Mus musculus", category = "C5",
                   subcategory = "GO:BP", geneid = "symbol",
                   method = c("AUCell", "UCell", "singscore", "ssgsea"),
                   aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                   kcdf = 'Gaussian')

# Check available assays
print("Available assays after gene set scoring:")
print(Seurat::Assays(MC))

# 2. Integrate differential expression gene sets ----

# Integration by celltype
print("Integrating results by celltype...")
result.dge <- irGSEA.integrate(object = MC,
                               group.by = "celltype",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell", "UCell", "singscore", "ssgsea"))

# Integration by group
print("Integrating results by group...")
result.dge_group <- irGSEA.integrate(object = MC,
                                    group.by = "group",
                                    metadata = NULL, col.name = NULL,
                                    method = c("AUCell", "UCell", "singscore", "ssgsea"))

# 3. Visualization ----

# 3.1 Heatmaps ----

# Select significant gene sets for celltype
geneset.show <- result.dge$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::pull(Name)

if (length(geneset.show) > 0) {
  irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                        method = "RRA",
                                        cluster.color = allcolour_celltype,
                                        show.geneset = geneset.show)
  print(irGSEA.heatmap.plot)
}

# Select significant gene sets for group
geneset.show2 <- result.dge_group$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::pull(Name)

if (length(geneset.show2) > 0) {
  irGSEA.heatmap.plot2 <- irGSEA.heatmap(object = result.dge_group,
                                         method = "RRA",
                                         cluster.color = allcolour_group,
                                         show.geneset = geneset.show2)
  print(irGSEA.heatmap.plot2)
}

# 3.2 Bubble plots ----

if (length(geneset.show) > 0) {
  irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge,
                                      method = "RRA",
                                      cluster.color = allcolour_celltype,
                                      show.geneset = geneset.show)
  print(irGSEA.bubble.plot)
}

if (length(geneset.show2) > 0) {
  irGSEA.bubble.plot2 <- irGSEA.bubble(object = result.dge_group,
                                       method = "RRA",
                                       cluster.color = allcolour_group,
                                       show.geneset = geneset.show2)
  print(irGSEA.bubble.plot2)
}

# 3.3 Upset plots ----

irGSEA.upset.plot <- irGSEA.upset(object = result.dge, 
                                  method = "RRA",
                                  mode = "intersect",
                                  cluster.color = allcolour_celltype,
                                  upset.width = 20,
                                  upset.height = 10,
                                  set.degree = 2,
                                  pt_size = grid::unit(2, "mm"))
print(irGSEA.upset.plot)

irGSEA.upset.plot2 <- irGSEA.upset(object = result.dge_group, 
                                   method = "RRA",
                                   mode = "intersect",
                                   cluster.color = allcolour_group,
                                   upset.width = 20,
                                   upset.height = 10,
                                   set.degree = 2,
                                   pt_size = grid::unit(2, "mm"))
print(irGSEA.upset.plot2)

# 3.4 Bar plots ----

irGSEA.barplot.plot <- irGSEA.barplot(object = result.dge,
                                      color.cluster = allcolour_celltype,
                                      method = c("AUCell", "UCell", "singscore",
                                                "ssgsea", "RRA"))
print(irGSEA.barplot.plot)

irGSEA.barplot.plot2 <- irGSEA.barplot(object = result.dge_group,
                                       color.cluster = allcolour_group,
                                       method = c("AUCell", "UCell", "singscore",
                                                 "ssgsea", "RRA"))
print(irGSEA.barplot.plot2)

# 3.5 Density scatter plots for specific pathways ----

# TNF-α signaling pathway
scatterplot_tnfa <- irGSEA.density.scatterplot(object = MC,
                                               method = "UCell",
                                               show.geneset = "HALLMARK-TNFA-SIGNALING-VIA-NFKB",
                                               reduction = "umap")
print(scatterplot_tnfa)

# Hypoxia pathway
scatterplot_hypoxia <- irGSEA.density.scatterplot(object = MC,
                                                  method = "UCell",
                                                  show.geneset = "HALLMARK-HYPOXIA",
                                                  reduction = "umap")
print(scatterplot_hypoxia)

# 3.6 Density heatmaps ----

densityheatmap_tnfa <- irGSEA.densityheatmap(object = MC,
                                            method = "UCell",
                                            show.geneset = "HALLMARK-TNFA-SIGNALING-VIA-NFKB")
print(densityheatmap_tnfa)

densityheatmap_hypoxia <- irGSEA.densityheatmap(object = MC,
                                               method = "UCell",
                                               show.geneset = "HALLMARK-HYPOXIA")
print(densityheatmap_hypoxia)

# 4. Traditional GO/KEGG analysis ----

# Find markers for traditional enrichment analysis
Idents(MC) <- "celltype"
sce.markers <- FindAllMarkers(MC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Convert gene symbols to ENTREZ IDs
ids <- mapIds(org.Mm.eg.db, 
              keys = sce.markers$gene, 
              column = "ENTREZID", 
              keytype = "SYMBOL", 
              multiVals = "first")

sce.markers$ENTREZID <- ids
sce.markers <- na.omit(sce.markers)
top20_MC <- sce.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# Split genes by cluster for pathway analysis
gene_lists <- split(sce.markers$ENTREZID, sce.markers$cluster)

# Save results
save(list = c("MC", "result.dge", "result.dge_group", "sce.markers", "gene_lists", "top20_MC"),
     file = "gene_set_scoring_results.RData")

# Summary statistics
print("Gene set scoring analysis completed!")
print(paste("Number of Hallmark gene sets with p < 0.05 (celltype):", length(geneset.show)))
print(paste("Number of Hallmark gene sets with p < 0.05 (group):", length(geneset.show2)))
print(paste("Number of marker genes found:", nrow(sce.markers)))
print(paste("Number of cell clusters:", length(unique(MC$celltype))))
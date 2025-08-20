# Mast cell analysis - Basic Analysis
# This script contains the basic clustering, annotation and visualization of mast cells
# Author: Based on analysis from README.md
# Date: 2024

# Load required libraries
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(clustree)
library(harmony)
library(tidydr)

# Set working directory (modify as needed)
# setwd("D:/data/Thymus")

# Load the initial Seurat object
# MC <- readRDS("./3cell.rds")

# 1) Recluster of Mast cells ----

# Subset to mast cells only
Mastcell <- subset(MC, celltype == "Mast cell")

# Cell cycle scoring
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(Mastcell))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(Mastcell))

# Score cell cycle phases
Mastcell <- CellCycleScoring(object = Mastcell,
                                     s.features = s_genes,
                                     g2m.features = g2m_genes,
                                     set.ident = TRUE)

# Plot cell cycle scores
Mastcell@meta.data %>% ggplot(aes(S.Score, G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()

# 1.1) Normalize and cluster from the start ----
Mastcell=CreateSeuratObject(counts = GetAssayData(Mastcell, assay = "RNA", slot = 'counts'),
                       meta.data = Mastcell@meta.data) 
Mastcell <- NormalizeData(Mastcell) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA

# Run harmony for batch correction
Mastcell <- RunHarmony(Mastcell, group.by.var = "group")
Mastcell <- FindNeighbors(Mastcell, reduction = "harmony", dims = 1:10)
Mastcell <- FindClusters(Mastcell, resolution = 0.2)
table(Mastcell@meta.data$seurat_clusters)
Mastcell <- RunUMAP(Mastcell, reduction = "harmony", dims = 1:10, return.model = TRUE)

# 1.2) Annotation for mast cell subtypes ----

# Find cluster markers
Idents(Mastcell) <- "seurat_clusters"
Mastcell_markers <-
  FindAllMarkers(
    object = Mastcell,
    test.use = "wilcox",
    only.pos = TRUE,
    logfc.threshold = 0.1,
    min.pct = 0.25,
    slot = "data"
  )
Mastcell.markers = Mastcell_markers %>% dplyr::select(gene, everything()) %>% subset(p_val < 0.05)
top20_Mastcell = Mastcell.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# Annotate cell types based on markers
Idents(Mastcell) <- "seurat_clusters"
celltype = data.frame(ClusterID = 0:4,
                      celltype = 'unknown')
celltype[celltype$ClusterID %in% c(0), 2] = 'Cycling MC'
celltype[celltype$ClusterID %in% c(1), 2] = "Mcpt9 medium MC"
celltype[celltype$ClusterID %in% c(2), 2] = "Mcpt9 high MC"
celltype[celltype$ClusterID %in% c(3), 2] = "Nr4a1 high MC"
celltype[celltype$ClusterID %in% c(4), 2] = "Lrmda+ MC"

sce.in = Mastcell
sce.in@meta.data$celltype = "NA"
for (i in 1:nrow(celltype)) {
  sce.in@meta.data[which(sce.in@meta.data$RNA_snn_res.0.2 == celltype$ClusterID[i]), 'celltype'] <-
    celltype$celltype[i]
}

Mastcell <- sce.in
rm(sce.in)

# 1.3) Visualization ----

# Define color palette
allcolour=  c(
"#D0AFC4",
"#89558D",
"#AFC2D9",
"#435B95",
"#79B99D",
"#D55640",
"#E69F84",
"#6CB8D2",
"#479D88",
"#415284",
"#C6367A",
"#ECDC52",
"#D1352B",
"#9B5B33"
  )
names(allcolour) <-
  c('Mcpt9 high MC', 'Mcpt9 medium MC', 'Cycling MC', 'Nr4a1 high MC', 'Lrmda+ MC', "IE-hpD07", "IE-hpD14", "LP-hpD0", "LP-hpD07", "LP-hpD14", "BMCP", "GMP", "DN(P)", "DN(Q)")

# UMAP plots
DimPlot(Mastcell, group.by = "celltype",pt.size = 2,cols = allcolour)+theme_dr()  + 
          guides(color = guide_legend(override.aes = list(size=5)))+ theme(panel.grid=element_blank(),
          legend.title = element_blank(), 
        legend.text = element_text(size=20), 
        legend.key.size=unit(1,'cm'))+ggtitle("Mast cell subclusters")+theme(plot.title = element_text(size = 20, hjust = 0.5))

# Tissue origin plot
tissuecolor <- c("#C6307C", "#4991C1")
names(tissuecolor) <- c("Epithelium", "LP")
Mastcell$tissue <- str_replace(Mastcell$tissue, "IE", "Epithelium")

DimPlot(Mastcell, group.by = "tissue",pt.size = 2,cols = tissuecolor)+theme_dr()  + 
          guides(color = guide_legend(override.aes = list(size=5)))+ theme(panel.grid=element_blank(),
          legend.title = element_blank(), 
        legend.text = element_text(size=20), 
        legend.key.size=unit(1,'cm'))+ggtitle("Mast cell origins")+theme(plot.title = element_text(size = 20, hjust = 0.5))

# Group plot
DimPlot(Mastcell, group.by = "group",pt.size = 2,cols = allcolour)+theme_dr()  + 
          guides(color = guide_legend(override.aes = list(size=5)))+ theme(panel.grid=element_blank(),
          legend.title = element_blank(), 
        legend.text = element_text(size=20), 
        legend.key.size=unit(1,'cm'))+ggtitle("Mast cell samples")+theme(plot.title = element_text(size = 20, hjust = 0.5))

# Feature and dot plots
FeaturePlot(MC, features = c("Mki67", "Lrmda", "Mcpt9", "Nr4a1"))
DotPlot(MC,features = c("Mki67", "Lrmda", "Mcpt9", "Nr4a1"), group.by = "group")
DotPlot(MC,features = c("Mki67", "Lrmda", "Mcpt9", "Nr4a1"), group.by = "celltype")

# 1.4) Cell proportion analysis ----

# Calculate cell proportion in different samples
Idents(Mastcell) <- "celltype" 
prop.table(table(Idents(Mastcell)))
table(Idents(Mastcell), Mastcell$group)

Cellratio <- prop.table(table(Idents(Mastcell), Mastcell$group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var1 <- factor(factor(Cellratio$Var1, levels = c("Cycling MC", "Lrmda+ MC", "Mcpt9 high MC", "Mcpt9 medium MC", "Nr4a1 high MC")))

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = allcolour)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"), legend.title= element_blank())+ggtitle("Mast cell subtype proportion in different samples")+theme (text = element_text(size =15),plot.title = element_text (hjust=0.5))

# Calculate samples proportion in different celltypes
Cellratio <- prop.table(table(Mastcell$group, Idents(Mastcell)), margin = 2)
Cellratio <- as.data.frame(Cellratio)

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = allcolour)+ 
  theme_classic() +
  labs(x='Celltype',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"), legend.title= element_blank())+ggtitle("Different sample proportion in Mast cell subtype")+theme (text = element_text(size =15),plot.title = element_text (hjust=0.5))

# Save the processed Mastcell object
# saveRDS(Mastcell, "processed_mastcells.rds")
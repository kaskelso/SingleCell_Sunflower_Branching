# Analysis of preliminary Sunflower single cell sequencing
# Kenny Askelson
# 08/13/2025

### PSC8 ###

library("Seurat")
library("ggplot2")
library("dplyr")
library("tidyr")

# Original file names looked like
# 'PSC8-DM-032725.scRNA.filtered.barcodes.tsv.gz', seurat doesn't like this so I need to remove PSC8-DM-032725.scRNA.filtered.
# from each file

data_dir <- 'PSC8_Dragen'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = expression_matrix, min.features = 300, min.cells = 3)

seurat_object

#features are genes and samples are cells

VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#hist of nCount_RNA

ggplot(seurat_object@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, fill = "lightgreen", color = "black") +
  geom_vline(aes(xintercept = 7000), colour="black") +
  labs(title = "Total Number of Counts (nCount_RNA)", x = "nCount_RNA", y = "Number of Cells")

#hist of nFeature_RNA

ggplot(seurat_object@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, fill = "red", color = "black") +
  geom_vline(aes(xintercept = 2500), colour="black") +
  labs(title = "Total Number of Counts (nFeature_RNA)", x = "nFeature_RNA", y = "Number of Cells")


seurat_object_filtered <- subset(seurat_object, subset = nFeature_RNA > 150 & nFeature_RNA < 2500 & nCount_RNA < 7000)

# Repeat Violin and Scatter plot

VlnPlot(seurat_object_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

FeatureScatter(seurat_object_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

seurat_object_filtered_normalized <- NormalizeData(seurat_object_filtered)

#Feature selection

seurat_object_filtered_normalized_selected <- FindVariableFeatures(seurat_object_filtered_normalized, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object_filtered_normalized_selected), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object_filtered_normalized_selected)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#PCA
all.genes <- rownames(seurat_object_filtered_normalized_selected)

seurat_object_filtered_normalized_selected_scaled <- ScaleData(seurat_object_filtered_normalized_selected, features = all.genes)

seurat_object_filtered_normalized_selected_scaled_PCA <- RunPCA(seurat_object_filtered_normalized_selected_scaled, features = VariableFeatures(object = seurat_object_filtered_normalized_selected_scaled))

print(seurat_object_filtered_normalized_selected_scaled_PCA[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat_object_filtered_normalized_selected_scaled_PCA, dims = 1:2, reduction = "pca")

DimPlot(seurat_object_filtered_normalized_selected_scaled_PCA, reduction = "pca") + NoLegend()

DimHeatmap(seurat_object_filtered_normalized_selected_scaled_PCA, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(seurat_object_filtered_normalized_selected_scaled_PCA)

#cut off at 6?

# Cell clustering

seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster <- FindNeighbors(seurat_object_filtered_normalized_selected_scaled_PCA, dims = 1:6)

seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster <- FindClusters(seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster, resolution = 0.5)

seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster <- RunUMAP(seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster, dims = 1:6)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster, reduction = "umap")

#Ha412HOChr10g0435621 was not found 

FeaturePlot(seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster, features = c("gene:Ha412HOChr10g0435441", "gene:Ha412HOChr10g0435451", "gene:Ha412HOChr10g0435491", "gene:Ha412HOChr10g0435511"))

VlnPlot(seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster, features = c("gene:Ha412HOChr10g0435441", "gene:Ha412HOChr10g0435451", "gene:Ha412HOChr10g0435491", "gene:Ha412HOChr10g0435511"))

seurat.markers <- FindAllMarkers(seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster, only.pos = TRUE)

seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(seurat_object_filtered_normalized_selected_scaled_PCA_cellcluster, features = top10$gene)

#looks like something is going in with gene:Ha412HOChr10g0435511 in cluster 4

### XRQ ###

data_dir_XRQ <- 'XRQ_Dragen'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix_XRQ <- Read10X(data.dir = data_dir_XRQ)
seurat_object_XRQ = CreateSeuratObject(counts = expression_matrix_XRQ, min.features = 300, min.cells = 3)

seurat_object_XRQ

#features are genes and samples are cells

VlnPlot(seurat_object_XRQ, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

FeatureScatter(seurat_object_XRQ, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#hist of nCount_RNA

ggplot(seurat_object_XRQ@meta.data, aes(x = nCount_RNA)) + xlim(0, 5000) +
  geom_histogram(bins = 100, fill = "lightgreen", color = "black") +
  geom_vline(aes(xintercept = 1000), colour="black") +
  labs(title = "XRQ Total Number of Counts (nCount_RNA)", x = "nCount_RNA", y = "Number of Cells")

#hist of nFeature_RNA

ggplot(seurat_object_XRQ@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(bins = 100, fill = "red", color = "black") +
  geom_vline(aes(xintercept = 1000), colour="black") +
  labs(title = "XRQ Total Number of Counts (nFeature_RNA)", x = "nFeature_RNA", y = "Number of Cells")

seurat_object_XRQ_filtered <- subset(seurat_object_XRQ, subset = nFeature_RNA > 145 & nFeature_RNA < 1000 & nCount_RNA > 160 & nCount_RNA < 1000)

VlnPlot(seurat_object_XRQ_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

FeatureScatter(seurat_object_XRQ_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Normalize 

seurat_object_XRQ_filtered_normalized <- NormalizeData(seurat_object_XRQ_filtered)

#Feature selection

seurat_object_XRQ_filtered_normalized_selected <- FindVariableFeatures(seurat_object_XRQ_filtered_normalized, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object_XRQ_filtered_normalized_selected), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object_XRQ_filtered_normalized_selected)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#PCA
all.genes_XRQ <- rownames(seurat_object_XRQ_filtered_normalized_selected)

seurat_object_XRQ_filtered_normalized_selected_scaled <- ScaleData(seurat_object_XRQ_filtered_normalized_selected, features = all.genes)

seurat_object_XRQ_filtered_normalized_selected_scaled_PCA <- RunPCA(seurat_object_XRQ_filtered_normalized_selected_scaled, features = VariableFeatures(object = seurat_object_XRQ_filtered_normalized_selected_scaled))

print(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA, dims = 1:2, reduction = "pca")

DimPlot(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA, reduction = "pca") + NoLegend()

DimHeatmap(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA)

# 10 seems like a decent cutoff? 

# Cell clustering

seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster <- FindNeighbors(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA, dims = 1:10)

seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster <- FindClusters(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster, resolution = 0.5)

seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster <- RunUMAP(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

DimPlot(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster, reduction = "umap")

FeaturePlot(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster, features = c("gene:Ha412HOChr10g0435441", "gene:Ha412HOChr10g0435451", "gene:Ha412HOChr10g0435491", "gene:Ha412HOChr10g0435511"))

VlnPlot(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster, features = c("gene:Ha412HOChr10g0435441", "gene:Ha412HOChr10g0435451", "gene:Ha412HOChr10g0435491", "gene:Ha412HOChr10g0435511"))

seurat_XRQ.markers <- FindAllMarkers(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster, only.pos = TRUE)

seurat_XRQ.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

seurat_XRQ.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(seurat_object_XRQ_filtered_normalized_selected_scaled_PCA_cellcluster, features = top10$gene)


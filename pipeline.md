# scRNA-seq analysis showcase (R Code)
**Author:** Shadi Mohammadabadi  
**Purpose:** This markdown file contains R code for an RNA-seq analysis pipeline originally written in Jupyter Notebook.

### 1. Load Libraries
```r
# Load required libraries
library(Seurat)
library(SeuratDisk)
library(scDblFinder)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(DropletUtils)
library(viridis)
library(plotly)
library(bluster)
library(data.table)
library(jsonlite)
library(Matrix)
library(RColorBrewer)
library(tibble)
library(ComplexHeatmap)
library(purrr)
library(scales)
library(circlize)
```
### 2. Helper Functions
```r
run_pca <- function(srat, npcs, nfeatures, normalization.method){
    srat <- NormalizeData(srat, normalization.method = normalization.method, verbose = F)
    srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = nfeatures, verbose = F)
    srat <- ScaleData(srat, verbose = F)
    srat <- RunPCA(srat, npcs = npcs, verbose = F)
    return(srat)
}
     
run_umap <- function(srat, npcs, resolution=0.8){
    srat <- FindNeighbors(srat, dims = 1:npcs)
    srat <- FindClusters(srat, resolution=resolution)
    srat <- RunUMAP(srat, reduction = "pca", dims = 1:npcs, verbose = F)
    return(srat)
}
     
run_doubletDetection <- function(srat_tr, nfeatures=1352, dims=20, iter=3, includePCs = 20, dbr=0.1, returnType = "sce"){
    sce <- as.SingleCellExperiment(srat_tr)
    sce <- scDblFinder(sce, clusters="ident", nfeatures=nfeatures, dims=dims, iter=iter, includePCs=includePCs, dbr=dbr, returnType=returnType)
    srat_tr <- as.Seurat(sce, data=NULL)
}
```
### 3. Add a list of gene markers
```r
# This is an example of gene markers
gene_markers <- list(
  "Excitatory neuron" = c("NRGN"),
  "Inhibitory neuron" = c("GAD1"),
  "Astrocytes" = c("GFAP"))
```

### 4. Read in the data
```r
# Read in the Seurat file
srat_tr <- LoadH5Seurat("path/to/public_data.h5Seurat")
```
### 5. Filter data
```r
#  filter for min.cells and min.features
m <- GetAssayData(srat_tr, layer = "counts")
cell_data <- srat_tr[[]][colnames(m),,drop=FALSE]
srat_tr <- CreateSeuratObject(counts = m, project = "RNAseq_Project", min.cells = 3, min.features = 200, meta.data = cell_data)

# Optional: filter for nFeature_RNA and nCount_RNA
srat_tr <- subset(srat_tr, subset = nFeature_RNA < 10000 & nCount_RNA < 5000)
```

### 6. Run clustering
```r
n_pcs <- 20 # Number of PCs for dimensionality reduction
var_features <- 2000 # Number of Variable Features
     

# Run clustering
srat_tr <- run_pca(srat_tr, npcs = n_pcs, nfeatures = var_features, normalization.method = "LogNormalize")
srat_tr <- run_umap(srat_tr, npcs = n_pcs)
```
### 7. Barcode Rank Plot
```r

bcrank <- barcodeRanks(counts(as.SingleCellExperiment(srat_tr)))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank
rank[uniq], bcrank$total[uniq], log="xy",
    xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
        col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
```

### 8. Run Doublet Detection
```r
# parameters
scDblFinder_nfeatures <- 2000
scDblFinder_dims <- 20
scDblFinder_iters <- 3
scDblFinder_includePCs <- 20
scDblFinder_dbr <- 0.1

# Run scDblFinder
srat_sce <- run_doubletDetection(srat_tr, nfeatures=scDblFinder_nfeatures, dims=scDblFinder_dims, iter=scDblFinder_iters, includePCs = scDblFinder_includePCs, dbr=scDblFinder_dbr)
```
### 9. Make a UMAP of doublets vs singlets
```r
DimPlot(srat_sce, reduction = "UMAP", group.by = "scDblFinder.class", label = FALSE) +
  ggtitle("scDblFinder.class before doublet removal")
```
### 10. Violin Plot
```r
VlnPlot(srat_sce, features = "nCount_RNA", group.by = "scDblFinder.class") +
  ggtitle("UMI Counts per Cell Grouped by scDblFinder Class")
```
### 11. Compare ground truth to doublet detection results
```r
# Simulated Ground Truth
# Optional: In case ground truth is available in metadata
confusion_matrix <- table(GroundTruth = srat_sce.groundtruth, Predicted = srat_sce.scDblFinder.class)

accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(paste("Accuracy:", accuracy))
precision <- confusion_matrix["doublet", "doublet"] / sum(confusion_matrix[, "doublet"])
print(paste("Precision:", precision))
recall <- confusion_matrix["doublet", "doublet"] / sum(confusion_matrix["doublet", ])
print(paste("Recall:", recall))
```
### 12. Dot plot of gene markers
```r
DotPlot(srat_sce, features = gene_markers, group.by="seurat_clusters")  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
ggtitle("Marker genes")
```
### 13. Top Markers
```r
# find top markers
markers <- FindAllMarkers(srat_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Feature plot of top markers for a selected cluster
FeaturePlot(subset(srat_sce, seurat_clusters == "5"), features = top_markers$gene[13:16], cols = viridis(10))
```
### 14. Complex Heatmap
```r
create_heatmap <- function(metadata, value, seurat_object, top_markers, width) {
  mat <- GetAssayData(seurat_object[, seurat_object[[metadata]] == value], layer = "scale.data")
  mat <- mat[top_markers$gene, ]
  
  # Calculate the min, max, and median values of the expression
  min_expr <- min(mat, na.rm = TRUE)
  max_expr <- max(mat, na.rm = TRUE)
  median_expr <- median(mat, na.rm = TRUE)
  # Create the heatmap
  Heatmap(
    mat, 
    name = value, 
    show_row_names = TRUE, 
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    column_title_rot = 45, 
    row_split = factor(rep(c("Excitatory", "Astrocytes", "Inhibitory"), split_sizes), 
                       levels = c("Excitatory", "Astrocytes", "Inhibitory")),
    cluster_row_slices = FALSE,
    col = colorRamp2(c(min_expr, median_expr, max_expr), c("purple", "black", "yellow")),
    column_title = value,
    width = width
  )
}
#Visualize a list of heatmaps
h2 <- create_heatmap("cell_type", "Excitatory neuron", data, top_markers, 2)
h3 <- create_heatmap("cell_type", "Astrocytes", data, top_markers, 2)
h4 <- create_heatmap("cell_type", "Inhibitory neuron", data, top_markers, 2)
ht_list <- h2 + h3 + h4
draw(ht_list, heatmap_legend_side = "right")
```
### 15. 3D UMAP
```r
umap_3d <- FindNeighbors(data, dims = 1:20)
umap_3d <- FindClusters(umap_3d, resolution=6)
umap_3d <- RunUMAP(umap_3d, reduction = "pca", dims = 1:20, verbose = F, n.components = 3)

umap_coordinates <- Embeddings(umap_3d, "umap")
umap_df <- as.data.frame(umap_coordinates)
umap_df.seurat_clusters <- umap_3d.seurat_clusters
     

fig <- plot_ly(umap_df, x = ~umap_1, y = ~umap_2, z = ~umap_3, 
               color = ~seurat_clusters, 
               colors = c('#636EFA', '#EF553B', '#00CC96')) %>% 
  add_markers(size=2) %>% 
  layout(scene = list(xaxis = list(title = 'UMAP1'), 
                      yaxis = list(title = 'UMAP2'), 
                      zaxis = list(title = 'UMAP3')))

# Show the plot
fig
```




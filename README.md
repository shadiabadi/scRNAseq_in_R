# **Single-Cell RNA-seq Analysis Showcase**  

**Author:** Shadi Mohammadabadi  
**Files:**
- **`pipeline.md`**: Main analysis pipeline.
- **`interactive_heatmap.md`**: Shiny app code for an interactive heatmap.

---

## **Overview**  
This repository demonstrates a single-cell RNA-seq analysis workflow, including:  
1. Quality control and filtering  
2. Dimensionality reduction and clustering (PCA, UMAP)  
3. Doublet detection using `scDblFinder`  
4. Visualization of marker gene expression (UMAP, violin plots, heatmaps)  

---

## **Requirements**  

### **Core Packages:**  
- `Seurat`, `scDblFinder`, `SingleCellExperiment`, `ggplot2`  
- For interactive plots: `plotly`, `shiny` (used in `interactive_heatmap.md`)  

---

## **Main File: `pipeline.md`**  

### **Key Steps:**  
- **Load Data:** Reads in `.h5Seurat` single-cell RNA-seq files.
- **Quality Control:** Filters cells based on feature counts.
- **Clustering:** Performs PCA and UMAP, followed by clustering.
- **Doublet Detection:** Uses `scDblFinder` to classify cells as singlets or doublets.
- **Marker Gene Analysis:** Identifies and visualizes cell types using known marker genes.

---

## **Interactive Heatmap App (`interactive_heatmap.md`)**  
The Shiny app generates interactive heatmaps of gene expression for selected clusters.  

---

## **How to Use This Repository:**  
1. Clone the repository:
   ```bash
   git clone https://github.com/shadiabadi/scRNAseq_in_R.git
   ```
2. Open `pipeline.md` in RStudio or a markdown viewer.
3. Edit data paths as needed for your dataset.

---

## **Sample Data Sources:**  
- [10X Genomics Public Datasets](https://www.10xgenomics.com/resources/datasets)  
- [GEO Datasets](https://www.ncbi.nlm.nih.gov/geo/)  

# Single-Cell RNA-Seq Analysis Workflow 🔬

This repository contains an R script for analyzing a single-cell RNA sequencing (scRNA-seq) dataset, specifically the **68k Peripheral Blood Mononuclear Cells (PBMCs) from Donor A**. The analysis pipeline uses the **Seurat** R package to perform key steps, from data loading to clustering and visualization.

---

## Prerequisites

You'll need the following R packages installed to run this analysis:

- `BiocManager`
- `rhdf5`
- `hdf5r`
- `Seurat`
- `patchwork`
- `ggplot2`

Install any missing packages using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("rhdf5")
install.packages("hdf5r")
install.packages("Seurat")

# Check and install additional packages
packages <- c("Seurat", "SeuratWrappers", "patchwork", "ggplot2")
installed <- packages %in% installed.packages()[,"Package"]
packages[!installed]
```

Step-by-Step Analysis
## 1. Set Up the Environment

Set the working directory and load the necessary libraries:

```r
# Set the working directory
setwd("C:/Users/cseng/Downloads/Final YR Proj")

# Load libraries
library(rhdf5)
library(hdf5r)
library(Seurat)
```

## 2. Load the Data

Load the single-cell expression matrix from a HDF5 file or Cell Ranger MEX format.
```r
# Option A: Load from HDF5 file
h5_path <- "filtered_feature_bc_matrix.h5"
mat <- Read10X_h5(h5_path)

# Option B: Load from Cell Ranger MEX format (recommended)
data_dir <- "C:/Users/cseng/Downloads/Final YR Proj/dataset/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices/filtered_matrices_mex/hg19"
seurat_data <- Read10X(data.dir = data_dir)
```

### 3. Create the Seurat Object
```r
# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = seurat_data)

# Summary of the object
seurat_obj
```

### 4. Quality Control (QC) and Filtering

Filter out low-quality cells based on:

Number of features (`nFeature_RNA`)
Total counts (`nCount_RNA`)
Mitochondrial percentage (`percent.mt`)

```r
# Calculate percent mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

### 5. Normalization and Dimensionality Reduction

Prepare the data for clustering:
```r
# Normalize data
seurat_obj <- NormalizeData(seurat_obj)

# Identify variable genes
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scale data
seurat_obj <- ScaleData(seurat_obj)

# Run PCA
seurat_obj <- RunPCA(seurat_obj)

# Visualize Elbow Plot
ElbowPlot(seurat_obj)
```

### 6. Clustering and Visualization

Cluster cells and visualize using UMAP.
```r
# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot UMAP
DimPlot(seurat_obj, reduction = "umap")
```

### 7. Find Marker Genes

Identify cluster-specific markers and visualize known cell-type markers.
```r
# Find marker genes
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers
head(markers)

# Visualize known markers
FeaturePlot(seurat_obj, features = c("CD3D", "MS4A1", "LYZ", "GNLY"))
```






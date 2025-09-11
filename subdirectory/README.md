# Single-Cell RNA-Seq Analysis with Seurat 🔬

This repository provides a complete **R workflow** for analyzing single-cell RNA sequencing (scRNA-seq) data using the **Seurat** package.  
It is demonstrated on the **68k Peripheral Blood Mononuclear Cells (PBMCs) from Donor A** dataset (10x Genomics).

The workflow covers:

- Reading raw **HDF5 (.h5)** files  
- Exploring the dataset structure  
- Creating a **Seurat object**  
- Performing **quality control (QC)**  
- **Normalization and dimensionality reduction (PCA)**  
- **Clustering and UMAP visualization**  
- Identifying **cluster-specific marker genes**  

---

## 📦 Prerequisites

Install required R packages:

```r
# Bioconductor manager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# HDF5 tools
BiocManager::install("rhdf5")   # recommended
install.packages("hdf5r")      # optional backup

# Seurat and visualization
install.packages("Seurat")
install.packages("patchwork")
install.packages("ggplot2")
```

## 🛠 Step-by-Step Workflow 
1. Set up environment
```r
getwd()
setwd("C:/Users/cseng/SecondONE_DRIVE/OneDrive/Documents")

# File picker to choose your .h5 file
h5_path <- file.choose()
h5_path
```

## 2. Explore HDF5 file
```r
library(rhdf5)

# List datasets inside
contents <- h5ls(h5_path, recursive = TRUE)
contents

# Create paths and find big datasets
contents$path <- file.path(contents$group, contents$name)
big <- subset(contents, otype == "H5I_DATASET" & grepl("x", dim))
big <- big[order(-as.numeric(sub("x.*","",big$dim))), ]
head(big, 10)

# Check file size in GB
file.info(h5_path)$size / 1024^3

# Free memory
gc()

```
### Optional inspection with `hdf5r`:
```r
library(hdf5r)
f <- H5File$new(h5_path, mode = "r")
f$ls(recursive = TRUE)
```


## 3. Inspect raw molecule information
```r
library(rhdf5)

# Read first 10,000 rows instead of full ~173M
n <- 10000
molecule_info <- data.frame(
  barcode = h5read(h5_path, "barcode", index = list(1:n)),
  gene    = h5read(h5_path, "gene", index = list(1:n)),
  umi     = h5read(h5_path, "umi", index = list(1:n))
)
head(molecule_info)

```

## 4. Create Seurat object
```r
library(Seurat)

# Load from Cell Ranger MEX format (recommended)
data_dir <- "C:/Users/cseng/Downloads/Final YR Proj/dataset/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices/filtered_matrices_mex/hg19"
seurat_data <- Read10X(data.dir = data_dir)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = seurat_data)

# View summary
seurat_obj

```

## 5. Quality control
```r
# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter low-quality cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

```

## 6. Normalization & dimensionality reduction
```r
# Normalize
seurat_obj <- NormalizeData(seurat_obj)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scale data
seurat_obj <- ScaleData(seurat_obj)

# PCA
seurat_obj <- RunPCA(seurat_obj)

# Elbow plot to choose PCs
ElbowPlot(seurat_obj)

```

## 7. Clustering & UMAP visualization
```r
# Build graph & cluster
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP projection
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot clusters
DimPlot(seurat_obj, reduction = "umap")

```

## 8. Find marker genes
```r
# Identify markers
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Show top markers
head(markers)

# Visualize known markers
FeaturePlot(seurat_obj, features = c("CD3D", "MS4A1", "LYZ", "GNLY"))

```

## 📊 Outputs

* UMAP plots showing clusters of single cells
* QC violin plots for RNA count, features, mitochondrial %
* Cluster-specific marker genes

## 📌 Notes

* This workflow is tailored for 10x Genomics PBMC 68k dataset, but works for any scRNA-seq dataset in HDF5 or MEX format.
* For large files (~173M rows), use partial loading (h5read() with index) to test before full runs.
* Adjust filtering thresholds (nFeature_RNA, percent.mt) depending on your dataset.

## 📚 References
[Seurat Documentation](https://satijalab.org/seurat/)<br>
[10x Genomics PBMC Dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets)


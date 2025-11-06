# -----------------------------------------------------------
# Project: ISS & scRNA-seq Data Exploration
# Author: Glenn Ross-Dolan
# Date: 2025-09-19
# -----------------------------------------------------------

# Load libraries
library(tidyverse) # CRAN
library(Matrix) # CRAN
library(ggplot2) # CRAN
library(data.table) # CRAN
library(Seurat) # CRAN
library(dplyr) # CRAN
library(viridis) # CRAN
library(sctransform) # CRAN
library(scry) # bioconductor
library(SingleCellExperiment) # bioconductor

# Set data directory
data_dir <- "/home/glennrdx/Documents/lab_projects/microenvironments/01_data/"

################################################################################
# scRNA-seq

# Read scRNA-seq expression data
scRNA_gene_expression <- fread(file.path(data_dir, "scRNA_gene_expression.txt")) %>% as.data.frame()
# Load cell type annotations
scRNA_cell_types <- fread(file.path(data_dir, "scRNA_cell_types.txt")) %>% as.data.frame()
  
  # transpose the data so genes are rows not columns (which coerses the data to a matrix). Then convert it back to a dataframe.
scRNA_gene_expression <- as.data.frame(t(scRNA_gene_expression))
  # CHECK THIS! - no barcodes in expression file so I assume they are in the same order as the meta file (and assign them to the expression matrix)
colnames(scRNA_gene_expression) <- scRNA_cell_types$Barcode

  # Clean the formatting of the annotation object
rownames(scRNA_cell_types) <- scRNA_cell_types$Barcode
scRNA_cell_types <- scRNA_cell_types[, "Annotation", drop = FALSE]

# Look at unique cell types
unique(scRNA_cell_types$Annotation)

# Create Seurat objects
scRNA_obj = CreateSeuratObject(counts = scRNA_gene_expression, meta.data = scRNA_cell_types)

# QC
# The [[ operator adds columns to object metadata. This is a great place to stash QC stats
scRNA_obj[["percent.mt"]] <- PercentageFeatureSet(scRNA_obj, pattern = "^MT-")

# 1. Histogram of nCount_RNA
ggplot(scRNA_obj@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Total Counts", x = "total counts", y = "Frequency")

# 2. Violin plot of percentage mitochondrial genes
VlnPlot(scRNA_obj, features = c('percent.mt'), pt.size = 0.01)

# 3. Scatter plot of three metrics together
ggplot(scRNA_obj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = Annotation)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10() +
  labs (title = 'seq depth vs gene count coloured by mt perc.',
        x = 'nCount_RNA',
        y = 'nFeature_RNA',
        color = 'Annotation')

# Detect low quality cells by MAD filtering.
# Define function to detect outliers using MAD:
is_outlier <- function(metric, nmads) {
  M <- metric
  med <- median(M, na.rm = T)
  mad_val <- mad(M, na.rm = T)
  
  outlier <- (M < (med - nmads * mad_val)) | (M > (med + nmads * mad_val))
  return(outlier)
}

# calculate log1p transformed metrics
scRNA_obj$log1p_nCount_RNA <- log1p(scRNA_obj$nCount_RNA)
scRNA_obj$log1p_nFeature_RNA <- log1p(scRNA_obj$nFeature_RNA)

# Calculate pct counts in top 20 genes for each cell
gene_sums <- rowSums(scRNA_obj@assays$RNA$counts) #if error try @Layers$counts
top20_genes <- names(sort(gene_sums, decreasing = T)[1:20])

top20_counts <- colSums(scRNA_obj@assays$RNA$counts[top20_genes, ]) # get a sum of the counts of the top 20 genes for each cell. 
scRNA_obj$pct_counts_top20genes <- (top20_counts / scRNA_obj$nCount_RNA) * 100

# Apply MAD-based filtering
scRNA_obj$outlier <- (
  is_outlier(scRNA_obj$log1p_nCount_RNA, 5) |
    is_outlier(scRNA_obj$log1p_nFeature_RNA, 5) |
    is_outlier(scRNA_obj$pct_counts_top20genes, 5) |
    is_outlier(scRNA_obj$percent.mt, 3) |
    (scRNA_obj$percent.mt > 8)
)

table(scRNA_obj$outlier)

# Filter out the outliers
cells_to_keep <- (!scRNA_obj$outlier)
scRNA_obj <- subset(scRNA_obj, cells = WhichCells(scRNA_obj)[cells_to_keep])

# REMOVE AMBIENT RNA
# REMOVE DOUBLETS

# Random subsampling to 100 cells per cell type
# extract metadata
meta <- scRNA_obj@meta.data %>%
  rownames_to_column(var = 'cell')

set.seed(123)
sampled_cells <- meta %>%
  group_by(Annotation) %>%
  slice_sample(n = 100) %>%
  pull(cell)

# Subset the Seurat object to those cell names
scRNA_obj <- subset(scRNA_obj, cells = sampled_cells)

# Normalisation
scRNA_obj <- SCTransform(scRNA_obj, verbose = F) # not sure if necessary to regress out pct mt variable
head(GetAssayData(scRNA_obj, slot = "data")) # confirm normalised

# Feature selection

# Convert Seurat object to a SCE object
sce <- Seurat::as.SingleCellExperiment(scRNA_obj, assay = "SCT")

# Calculate deviance feature selection
sce_deviance <- devianceFeatureSelection(sce, assay = 'counts')

# Extract binomial deviance values
binomial_deviance <- rowData(sce_deviance)$binomial_deviance

# Select top highly deviant genes
n_top_genes <- 4000
highly_deviant_genes <- names(sort(binomial_deviance, decreasing = TRUE)[1:n_top_genes])

# Add deviance information to the Seurat object
# add binomial deviance scores to meta.data 
scRNA_obj[["SCT"]]@meta.features$binomial_deviance <- binomial_deviance[rownames(scRNA_obj)]

# Create a logical vector for the highly deviant genes
highly_deviant_logical <- rownames(scRNA_obj) %in% highly_deviant_genes
scRNA_obj[["SCT"]]@meta.features$highly_deviant <- highly_deviant_logical

# Visualise Feature Selection results

#calculate mean expression / dispersion (using normalised data)
normalised_counts <- GetAssayData(scRNA_obj, slot = "data", assay = "SCT")
gene_means <- Matrix::rowMeans(normalised_counts)
gene_vars <- apply(normalised_counts, 1, var)
gene_dispersions <- gene_vars / gene_means

# create a df for plotting
plot_data <- data.frame(
  gene = rownames(scRNA_obj),
  means = gene_means,
  dispersions = gene_dispersions,
  highly_deviant = highly_deviant_logical,
  binomial_deviance = binomial_deviance[rownames(scRNA_obj)]
)

# create scatter plot
ggplot(plot_data, aes(x = means, y = dispersions, color = highly_deviant)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "lightgray", "TRUE" = "red")) +
  xlim(0, 1.5) +
  ylim(0, 3) +
  labs(
    title = "Gene Selection by Deviance",
    x = "Mean Expression",
    y = "Dispersion",
    color = "Highly Deviant"
  ) +
  theme_minimal()

# Set highly deviant genes as variable features
VariableFeatures(scRNA_obj) <- highly_deviant_genes

# Dimensionality Reduction

#PCA
scRNA_obj <- RunPCA(scRNA_obj, features = VariableFeatures(scRNA_obj), verbose = F)
pca_feature_plot <- FeaturePlot(scRNA_obj, reduction = "pca", features = "nCount_RNA") +
  ggtitle("PCA - Total Counts")

print(pca_feature_plot)

# t-SNE
scRNA_obj <- RunTSNE(scRNA_obj, dims = 1:50, reduction = "pca", verbose = F) #make elbow plot and judge

tsne_plot <- DimPlot(scRNA_obj, reduction = "tsne", group.by = "Annotation") +
  ggtitle("t-SNE - Cell Type Annotations")

print(tsne_plot)

# UMAP
#calc neighbours
scRNA_obj <- FindNeighbors(scRNA_obj, dims = 1:50, reduction = "pca", verbose = F)

# Run UMAP
scRNA_obj <- RunUMAP(scRNA_obj, dims = 1:50, reduction = "pca", verbose = F)

umap_plot <- DimPlot(scRNA_obj, reduction = "umap", group.by = "Annotation") +
  ggtitle("UMAP - Cell Type Annotations")

print(umap_plot)

###############################################################################

# ISS
ISS_gene_expression <- fread(file.path(data_dir, "ISS_gene_expression.txt")) %>% as.data.frame()
ISS_cell_types   <- fread(file.path(data_dir, "ISS_cell_types.txt")) %>% as.data.frame()

# Load ISS cell meta data with spatial coordinates
cells <- fread(file.path(data_dir, "cells.csv")) %>% as.data.frame()

# Look at unique cell types
ISS_cell_types <- as.data.frame(t(ISS_cell_types))
ISS_cell_types$V1 <- NULL
colnames(ISS_cell_types) <- 'Annotation'
unique(ISS_cell_types$Annotation)

# Create Seurat objects
ISS_obj = CreateSeuratObject(counts = ISS_gene_expression, meta.data = ISS_cell_types)

## QC
# Not sure what QC is required as there are only 314 genes

## Normalisation - log transformed values
# Not sure


# plot data
cells <- cbind(cells, ISS_cell_types)
# Basic spatial scatter plot
ggplot(cells, aes(x = x_centroid, y = y_centroid)) +
  geom_point(aes(color = Annotation), size = 0.3, alpha = 0.7) +
  scale_size_continuous(name = "Cell Area", range = c(0.5, 3)) +
  scale_x_reverse() + 
  theme_minimal() +
  labs(title = "Spatial Data Coloured By Cell Type",
       x = "X Coordinate", y = "Y Coordinate")

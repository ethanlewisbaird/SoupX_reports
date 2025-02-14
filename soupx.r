.libPaths(c("/data/ebaird/R-packages", .libPaths()))
library(SoupX)
library(ggplot2)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
sample_dir <- file.path("/data/ebaird/scRNAseq/CR/outputs", sample_name, "outs")
output_dir <- file.path("data/ebaird/scRNAseq/soupx/outs", sample_name)

# Load data and estimate soup profile
sc = load10X(file.path(sample_dir))
sc = autoEstCont(sc)
out = adjustCounts(sc)

# Create Seurat objects from the original and adjusted counts
seurat_obj_original <- CreateSeuratObject(counts = sc$toc)
seurat_obj_adjusted <- CreateSeuratObject(counts = out)

# Normalize, find variable features, scale data, and run PCA for original data
seurat_obj_original <- NormalizeData(seurat_obj_original, verbose = FALSE)
seurat_obj_original <- FindVariableFeatures(seurat_obj_original, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seurat_obj_original <- ScaleData(seurat_obj_original, verbose = FALSE)
seurat_obj_original <- RunPCA(seurat_obj_original, npcs = 20, verbose = FALSE)

# Normalize, find variable features, scale data, and run PCA for adjusted data
seurat_obj_adjusted <- NormalizeData(seurat_obj_adjusted, verbose = FALSE)
seurat_obj_adjusted <- FindVariableFeatures(seurat_obj_adjusted, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seurat_obj_adjusted <- ScaleData(seurat_obj_adjusted, verbose = FALSE)
seurat_obj_adjusted <- RunPCA(seurat_obj_adjusted, npcs = 20, verbose = FALSE)

# Extract PCA coordinates for original and adjusted data
pca_coords_original <- Embeddings(seurat_obj_original, "pca")[, 1:2]
pca_coords_adjusted <- Embeddings(seurat_obj_adjusted, "pca")[, 1:2]

# Create data frames for ggplot
data_original <- data.frame(
  PC1 = pca_coords_original[, 1],
  PC2 = pca_coords_original[, 2]
)

data_adjusted <- data.frame(
  PC1 = pca_coords_adjusted[, 1],
  PC2 = pca_coords_adjusted[, 2]
)

# Create scatter plots
p_original <- ggplot(data_original, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA Plot of Original Data", x = "PC1", y = "PC2") +
  theme_minimal()

p_adjusted <- ggplot(data_adjusted, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA Plot of Adjusted Data", x = "PC1", y = "PC2") +
  theme_minimal()

# Save the plots as PDFs
ggsave(file.path(output_dir, "PCA_Plot_of_Original_Data.pdf"), plot = p_original, device = "pdf")
ggsave(file.path(output_dir, "PCA_Plot_of_Adjusted_Data.pdf"), plot = p_adjusted, device = "pdf")

# Print the plots
print(p_original)
print(p_adjusted)
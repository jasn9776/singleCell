install.packages(c("Seurat", "tidyverse", "caret"))
install.packages(c("caret"))

library(Seurat)
library(tidyverse)
library(caret)

# Load your single-cell RNA-seq data
# Replace 'your_data_file.csv' with the actual file name/path
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# store the percent mitochondria for each cell
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalise (default is log normalise)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# identify top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)

VariableFeatures(pbmc)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

plot1
plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# # repeat for CLR
# normalized_seurat_clr <- FindVariableFeatures(normalized_seurat_clr, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(normalized_seurat_clr)
# normalized_seurat_clr <- ScaleData(normalized_seurat_clr, features = all.genes)
# normalized_seurat_clr <- RunPCA(normalized_seurat_clr, features = VariableFeatures(object = normalized_seurat_clr))
# normalized_seurat_clr <- FindNeighbors(normalized_seurat_clr, dims = 1:10)
# normalized_seurat_clr <- FindClusters(normalized_seurat_clr, resolution = 0.5)
# normalized_seurat_clr <- RunUMAP(normalized_seurat_clr, dims = 1:10)
# DimPlot(normalized_seurat_clr, reduction = "umap")

# assign canonical markers based on website: 
# https://satijalab.org/seurat/articles/pbmc3k_tutorial#finding-differentially-expressed-features-cluster-biomarkers

levels(pbmc)
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
canonical <- Idents(pbmc)
canonical[5]

library(sctree)
# Cross-validation for cell type prediction
# Split data into training and testing sets
set.seed(123)

# reduced_pbmc <- Reductions(pbmc)
reduced_pbmc <- pbmc[["pca"]]@cell.embeddings
counts.df <- as.data.frame(reduced_pbmc)


splitIndex <- createDataPartition(canonical, p = 0.8, list = FALSE)
train_data <- counts.df[splitIndex, ]
test_data <- counts.df[-splitIndex, ]

y_train = canonical[splitIndex]
y_test = canonical[-splitIndex ]

Your_Model <- train(x=train_data,y=y_train, data = train_data, method = "rf")
pbmc

predictions <- predict(Your_Model, newdata = test_data)

confusion_matrix <- confusionMatrix(predictions, y_test)
print(confusion_matrix)

#repeat with UMAP
reduced_umap <- pbmc[["umap"]]@cell.embeddings
umap_counts.df = reduced_umap
train_data <- counts.df[splitIndex, ]
test_data <- counts.df[-splitIndex, ]
Your_Model <- train(x=train_data,y=y_train, data = train_data, method = "rf")
pbmc

predictions <- predict(Your_Model, newdata = test_data)
confusion_matrix <- confusionMatrix(predictions, y_test)
print(confusion_matrix)


# # Repeat for CLR normalisation XDD
# 
# # reduced_pbmc <- Reductions(pbmc)
# reduced_pbmc_clr <- normalized_seurat_clr[["pca"]]@cell.embeddings
# reduced_umap_clr <- normalized_seurat_clr[["umap"]]@cell.embeddings
# counts.df <- as.data.frame(reduced_pbmc_clr)
# 
# canonical <- normalized_seurat_clr@meta.data[["seurat_clusters"]]
# 
# splitIndex <- createDataPartition(canonical, p = 0.8, list = FALSE)
# train_data <- counts.df[splitIndex, ]
# test_data <- counts.df[-splitIndex, ]
# 
# pbmc$groups
# 
# train_data@assays$RNA$counts
# 
# train_data
# 
# y_train = canonical[splitIndex]
# y_test = canonical[-splitIndex ]
# 
# Your_Model <- train(x=train_data,y=y_train, data = train_data, method = "rf")
# pbmc
# 
# predictions <- predict(Your_Model, newdata = test_data)
# 
# confusion_matrix <- confusionMatrix(predictions, y_test)
# print(confusion_matrix)
# 
# # Repeat for RC normalisation XDD
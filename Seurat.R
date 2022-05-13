#install.packages('Seurat')
library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)

# Read the matrix and create Seurat object
mtx <- Matrix::readMM("C:/Users/Александр/Desktop/Data/R/filtered_feature_bc_matrix/matrix.mtx")
barcodes <- readr::read_tsv("C:/Users/Александр/Desktop/Data/R/filtered_feature_bc_matrix/barcodes.tsv", col_names = FALSE)
genes <- readr::read_tsv("C:/Users/Александр/Desktop/Data/R/filtered_feature_bc_matrix/features.tsv", col_names = FALSE)
colnames(mtx) <- barcodes$X1
rownames(mtx) <- genes$X2

# Check that the column contains nothing but gene expression data
table(genes$X3)
mtx
genes$X2

# Create seurat object. sc means single cell
sc <- CreateSeuratObject(mtx, min.cells = 3, min.features = 200)
sc
table(genes$X3)

# QC and selecting cells for further analysis
# QC - Quality Control

# mitochondrial genes analysis
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
head(sc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

# Cell selection
# We select cells that have unique feature counts over 2000 or less than 10000
# and less than 5% of mitochondrial genes expressed
# This is done to discard empty droplets or double-celled droplets or dead cells
sc2 <- subset(sc, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 5)
sc2
head(sc2)

#Normalisation
sc3 <- NormalizeData(sc2, normalization.method = "LogNormalize", scale.factor = 10000)

# Identifying 2000 most variable features

sc4 <- FindVariableFeatures(sc3, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc4), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc4)
plot2 <- LabelPoints(plot = plot1, points = "ACTA2", repel = TRUE)
plot1
top10


# Scaling the data
all.genes <- rownames(sc4)
sc5 <- ScaleData(sc4, features = all.genes)

# Linear dimensional reduction
sc6 <- RunPCA(sc5, features = VariableFeatures(object = sc5))
# Examine and visualize PCA results a few different ways
print(sc6[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(sc6, reduction = "pca")

#Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
sc7 <- JackStraw(sc6, num.replicate = 100)
sc7 <- ScoreJackStraw(sc7, dims = 1:20)
JackStrawPlot(sc7, dims = 1:20)
ElbowPlot(sc7)

#Clusterisation. I suggest first 10 PCAs enough sufficient according to graphs
sc8 <- FindNeighbors(sc7, dims = 1:10)
sc9 <- FindClusters(sc8, resolution = 0.5)
head(Idents(sc9), 5)


#Umap 
#Lets install umap first if needed:
#reticulate::py_install(packages = 'umap-learn')

sc10 <- RunUMAP(sc9, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(sc10, reduction = "umap")

FeaturePlot(sc10, features = c("ACTA2"))
VlnPlot(sc10, features = c("ACTA2"))

#save rds-file
saveRDS(sc10, file = "C:/Users/Александр/Desktop/Data/R/final/seurat_result.rds")
sc10 <- readRDS(file = "C:/Users/Александр/Desktop/Data/R/final/seurat_result.rds")

#Count global Foldchange
sc10.markers <- FindAllMarkers(sc10, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0, return.thresh = 1)
sc10.markers
sc10.markers %>%
  group_by(cluster)
  #slice_max(n = 2, order_by = avg_log2FC)
head(sc10)
sc10.markers
write.csv(sc10.markers,'C:/Users/Александр/Desktop/Data/R/global_fold.csv', row.names = FALSE)

# find all markers of cluster 4 (suspicious)
cluster4.markers <- FindMarkers(sc10, ident.1 = 4, min.pct = 0, only.pos = FALSE, logfc.threshold = 0, return.thresh = 1)
head(cluster4.markers, n = 5)
write.csv(cluster4.markers,'C:/Users/Александр/Desktop/Data/R/final/cluster4_markers.csv', row.names = TRUE)

# find all markers of cluster 0 (non-responding)
cluster0.markers <- FindMarkers(sc10, ident.1 = 0, min.pct = 0, only.pos = FALSE, logfc.threshold = 0, return.thresh = 1)
head(cluster0.markers, n = 5)
write.csv(cluster0.markers,'C:/Users/Александр/Desktop/Data/R/final/cluster0_markers.csv', row.names = TRUE)

# find all markers of cluster 3 (responding)
cluster3.markers <- FindMarkers(sc10, ident.1 = 3, min.pct = 0, only.pos = FALSE, logfc.threshold = 0, return.thresh = 1)
head(cluster3.markers, n = 5)
write.csv(cluster3.markers,'C:/Users/Александр/Desktop/Data/R/final/cluster3_markers.csv', row.names = TRUE)


#Compare clusters 0 and 3
cluster0.markers <- FindMarkers(sc10, ident.1 = 0, ident.2 = 3, min.pct = 0, only.pos = FALSE, logfc.threshold = 0, return.thresh = 1)
write.csv(cluster0.markers,'C:/Users/Александр/Desktop/Data/R/final/local_fold_0vs3.csv', row.names = TRUE)


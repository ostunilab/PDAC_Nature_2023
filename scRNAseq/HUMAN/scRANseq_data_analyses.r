####### scRNA-seq analyses template
library(Seurat)
library(SeuratWrappers)
library(scDblFinder)


# Read 10x data for each sample
raw_counts <- Read10X("filtered_feature_bc_matrix/")
Seurat.obj<- CreateSeuratObject(raw_counts, min.cells = 3,  min.features = 200)

# Doublets prediction for each sample
sceDblF <- scDblFinder(Seurat.obj@assays$RNA@counts,dbr =0.05)
score.scdblfinder <- sceDblF@colData@listData[["scDblFinder.score"]]
names(score.scdblfinder) <- rownames(sceDblF@colData)

# Merge data with merge() function from Seurat
# Seurat.obj.unfiltered.rds object contains merged dataset with Neutrophils and doublets annotated

human.obj <- readRDS("Seurat.obj.unfiltered.rds")


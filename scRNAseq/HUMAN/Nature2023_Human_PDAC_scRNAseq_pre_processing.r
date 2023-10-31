#### HUMAN TUMOUR scRNA-seq DATA ####

library(Seurat)
library(ggplot2)
library(scDblFinder)
library(copykat)
library(SeuratWrappers)
library(harmony)
library(foreach)
library(parallel)
library(dplyr)
library(tidyr)
library(magrittr)
library(nichenetr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(slingshot)
library(viridis)
library(scales)
library(msigdbr)

#### PRE-PROCESSING ####

#load datasets
Sample.data <- Read10X("GSM6727545/filtered_feature_bc_matrix/")
Sample_30_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "LPDAC_30_tumor")

Sample.data <- Read10X("GSM6727548/filtered_feature_bc_matrix")
Sample_50_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "PDAC_50_tumor")

Sample.data <- Read10X("GSM6727549/filtered_feature_bc_matrix")
Sample_51_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "PDAC_51_tumor")

Sample.data <- Read10X("GSM6727550/filtered_feature_bc_matrix")
Sample_55_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "PDAC_55_tumor")

Sample.data <- Read10X("GSM6727546/filtered_feature_bc_matrix")
Sample_47_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "PDAC_47_tumor")

Sample.data <- Read10X("GSM6727547/filtered_feature_bc_matrix")
Sample_48_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "PDAC_48_tumor")

Sample.data <- Read10X("GSM6727543/filtered_feature_bc_matrix/")
Sample_25_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "LPDAC_25_Tumor")

Sample.data <- Read10X("GSM6727551/filtered_feature_bc_matrix/")
Sample_60_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "PDAC_60_Tumor")

Sample.data <- Read10X("GSM6727544/filtered_feature_bc_matrix/")
Sample_26_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "LPDAC_26_Tumor")

Sample.data <- Read10X("GSM6727542/filtered_feature_bc_matrix/")
Sample_15_tumor <- CreateSeuratObject(Sample.data, min.cells = 3,  project = "LPDAC_15_Tumor")

#merge samples
Sample.merge<- merge(Sample_30_tumor, y = c(Sample_50_tumor,  Sample_51_tumor, Sample_55_tumor, Sample_47_tumor,Sample_60_tumor, Sample_48_tumor, Sample_25_tumor, Sample_26_tumor, Sample_15_tumor), add.cell.ids = c("LPDAC_30_tumor", "PDAC_50_tumor", "PDAC_51_tumor", "PDAC_55_tumor", "PDAC_47_tumor","PDAC_60_tumor", "PDAC_48_tumor", "LPDAC_25_tumor", "LPDAC_26_tumor", "LPDAC_15_tumor"), project = "humanPDAC")

# cn prediction with copykat
copykat.PDAC48 <- copykat(rawmat=as.matrix(Sample_48_tumor@assays$RNA@counts), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FALSE")
copykat.PDAC60 <- copykat(rawmat=as.matrix(Sample_60_tumor@assays$RNA@counts), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FLASE")
copykat.LPDAC25 <- copykat(rawmat=as.matrix(Sample_25_tumor@assays$RNA@counts), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FLASE")
copykat.LPDAC26 <- copykat(rawmat=as.matrix(Sample_26_tumor@assays$RNA@counts), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FALSE", cell.line=T)
copykat.LPDAC30 <- copykat(rawmat=as.matrix(Sample_30_tumor@assays$RNA@counts),, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FLASE")
copykat.PDAC47 <- copykat(rawmat=as.matrix(Sample_47_tumor@assays$RNA@counts), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FALSE")
copykat.PDAC55 <- copykat(rawmat=as.matrix(Sample_55_tumor@assays$RNA@counts), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FALSE")
copykat.PDAC51 <- copykat(rawmat=as.matrix(Sample_51_tumor@assays$RNA@counts), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FALSE")
copykat.PDAC50 <- copykat(rawmat=as.matrix(Sample_50_tumor@assays$RNA@counts), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FALSE")

# add percentage of expression of mitochondrial genes and ribosomal protein genes 
Sample.merge <- PercentageFeatureSet(Sample.merge, pattern = "^MT-", col.name = "percent.mito")
Sample.merge <- PercentageFeatureSet(Sample.merge, pattern = "^RPL", col.name = "percent.ribo")

#filtering
Sample.merge <- subset(Sample.merge, subset = percent.mt < 40 & nCount_RNA > 1000 & nFeature_RNA > 500)

#filtering for Neutrophils annotation
#Sample.merge <- subset(Sample.merge, subset = percent.mt < 40 & nFeature_RNA > 200)

#cell-cycle prediction
s.genes <- readLines('genes_Sphase.txt')
g2m.genes <- readLines('genes_G2Mphase.txt')
Sample.merge <- CellCycleScoring(Sample.merge, g2m.features=g2m.genes[g2m.genes %in% rownames(Sample.merge@assays$RNA@data)], s.features=s.genes[s.genes %in% rownames(Sample.merge@assays$RNA@data)], set.ident = FALSE)
Sample.merge@meta.data$CC.Difference <- Sample.merge@meta.data$S.Score - Sample.merge@meta.data$G2M.Score

#doublet calling with scDblFinder
doublets.scdblfinder <- unlist(lapply(unique(Sample.merge$orig.ident), function(x) {
sel_cells <- rownames(Sample.merge@meta.data[which(Sample.merge@meta.data$orig.ident == x),])
sceDblF <- scDblFinder(Sample.merge@assays$RNA@counts[,sel_cells],dbr =0.05)
doublets_anno <- as.vector(sceDblF@colData$scDblFinder.class)
names(doublets_anno) <- row.names(sceDblF@colData)
return(doublets_anno)
}))

#filtering doublets
Sample.merge <- AddMetaData(Sample.merge, doublets.scdblfinder[colnames(Sample.merge)], "is.doublet")
Sample.merge <- subset(Sample.merge, subset = is.doublet == 'singlet')

#filtering mitochondrial genes and ribosomal protein genes
mito.genes.expr <- grep("^MT-", rownames(Sample.merge@assays$RNA@counts), value = T)
ribo.genes.expr <- grep("^RPL", rownames(Sample.merge@assays$RNA@counts), value = T)
keep_genes = rownames(Sample.merge@assays$RNA@counts)
keep_genes = keep_genes[!(keep_genes %in% c(mito.genes.expr,ribo.genes.expr))]
Sample.merge <- subset(Sample.merge, features = keep_genes)

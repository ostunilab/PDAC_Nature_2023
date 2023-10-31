library(Seurat)
library(SeuratWrappers)
library(harmony)
library(dplyr)
set.seed(123)

### ALL CELLS ###

Sample_expr <- NormalizeData(Sample_expr, normalization.method = "LogNormalize", scale.factor = 1e4, assay='RNA')
Sample_expr <- RunFastMNN(object.list = SplitObject(Sample_expr, split.by = "orig.ident"))
Sample_expr <- RunUMAP(Sample_expr, reduction='mnn', dims = 1:30)
Sample_expr <- FindNeighbors(Sample_expr, reduction = 'mnn', dims = 1:30)
Sample_expr <- FindClusters(Sample_expr, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))

for (i in c(0,15,4,10,12,9,11,18,2,13,17)){
  eval(parse(text=(paste("sub_cl <- subset(Sample_expr, subset = RNA_snn_res.0.8 == ",i,")",sep=""))))
  Idents(sub_cl) <- sub_cl$orig.ident
  eval(parse(text=(paste("Cluster_",i,"_WT_vs_KO <- FindMarkers(sub_cl, ident.1 ='WT', ident.2='KO', min.pct=0.1, only.pos = FALSE, pseudocount.use = 0.1, logfc.threshold = 0.5, assay = 'RNA')",sep=""))))
  #eval(parse(text=(paste("write.table(Cluster_",i,"_WT_vs_KO[Cluster_",i,"_WT_vs_KO$p_val_adj<0.01,], 'ALLCELLS_Cluster",i,"_DEG_WT_vs_COX2_KO_res0.8.txt', sep='\t', quote=F, col.names=T, row.names=T)",sep=""))))
}

sub_cl <- subset(Sample_expr, subset = RNA_snn_res.0.8 %in% c(1,3,5,8))
Idents(sub_cl) <- sub_cl$orig.ident
Cluster_Macro_WT_vs_KO <- FindMarkers(sub_cl, ident.1 ='WT', ident.2='KO', min.pct=0.1, only.pos = FALSE, pseudocount.use = 0.1, logfc.threshold = 0.5, assay = 'RNA')
#write.table(Cluster_Macro_WT_vs_KO[Cluster_Macro_WT_vs_KO$p_val_adj<0.01,], 'ALLCELLS_ClusterMacro_DEG_WT_vs_COX2_KO_res0.8.txt', sep='\t', quote=F, col.names=T, row.names=T)

sub_cl <- subset(Sample_expr, subset = RNA_snn_res.0.8 %in% c(6,16))
Idents(sub_cl) <- sub_cl$orig.ident
Cluster_Fibroblasts_WT_vs_KO <- FindMarkers(sub_cl, ident.1 ='WT', ident.2='KO', min.pct=0.1, only.pos = FALSE, pseudocount.use = 0.1, logfc.threshold = 0.5, assay = 'RNA')
#write.table(Cluster_Fibroblasts_WT_vs_KO[Cluster_Fibroblasts_WT_vs_KO$p_val_adj<0.01,], 'ALLCELLS_ClusterFibroblasts_DEG_WT_vs_COX2_KO_res0.8.txt', sep='\t', quote=F, col.names=T, row.names=T)

sub_cl <- subset(Sample_expr, subset = RNA_snn_res.0.8 %in% c(7,14))
Idents(sub_cl) <- sub_cl$orig.ident
Cluster_DCs_WT_vs_KO <- FindMarkers(sub_cl, ident.1 ='WT', ident.2='KO', min.pct=0.1, only.pos = FALSE, pseudocount.use = 0.1, logfc.threshold = 0.5, assay = 'RNA')
#write.table(Cluster_DCs_WT_vs_KO[Cluster_DCs_WT_vs_KO$p_val_adj<0.01,], 'ALLCELLS_ClusterDCs_DEG_WT_vs_COX2_KO_res0.8.txt', sep='\t', quote=F, col.names=T, row.names=T)

### TUMOR-ASSOCIATED MACROPHAGES ###

Sample_expr_TAM <- subset(Sample_expr, subset = Annotation_2 == 'TAMs')

Sample_expr_TAM <- NormalizeData(Sample_expr_TAM, normalization.method = "LogNormalize", scale.factor = 1e4, assay='RNA')
Sample_expr_TAM <- FindVariableFeatures(Sample_expr_TAM,selection.method = "vst", nfeatures = 3000)
Sample_expr_TAM <- ScaleData(Sample_expr_TAM, vars.to.regress = c("CC.Difference"))
Sample_expr_TAM <- RunPCA(Sample_expr_TAM)
Sample_expr_TAM <- RunHarmony(Sample_expr_TAM, group.by.vars = c('orig.ident'), dims.use = 1:30, theta=2, reduction.save = 'harmony')
Sample_expr_TAM <- RunUMAP(Sample_expr_TAM, reduction='harmony', dims = 1:20)
Sample_expr_TAM <- FindNeighbors(Sample_expr_TAM, reduction = 'harmony', dims = 1:20)
Sample_expr_TAM <- FindClusters(Sample_expr_TAM, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))

sub_Il1bTAM <- subset(Sample_expr_TAM, subset = Annotation_TAMs == 'Il1b_TAMs')
Idents(sub_Il1bTAM) <- sub_Il1bTAM$orig.ident
Il1bTAM_WT_vs_KO <- FindMarkers(sub_Il1bTAM, ident.1 ='WT', ident.2='KO', min.pct=0.1, only.pos = FALSE, pseudocount.use = 0.1, logfc.threshold = 0, assay = 'RNA')

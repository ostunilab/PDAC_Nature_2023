library(Seurat)
library(SeuratWrappers)
library(harmony)
library(SeuratExtend)
library(parallel)
library(foreach)
library(dplyr)
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)
library(msigdbr)
set.seed(123)


#### ALL CELLS ####

Sample_expr <- NormalizeData(Sample_expr, normalization.method = "LogNormalize", scale.factor = 1e4, assay='RNA')
Sample_expr_FastMNN <- RunFastMNN(object.list = SplitObject(Sample_expr, split.by = "orig.ident"))
Sample_expr_FastMNN <- RunUMAP(Sample_expr_FastMNN, reduction='mnn', dims = 1:20)
Sample_expr_FastMNN <- FindNeighbors(Sample_expr_FastMNN, reduction = 'mnn', dims = 1:20)
Sample_expr_FastMNN <- FindClusters(Sample_expr_FastMNN, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))

Idents(Sample_expr_FastMNN) <- 'RNA_snn_res.0.5'
i <- 0
while(i<=23){
  eval(parse(text=(paste("cluster",i,".markers0.5 <- FindMarkers(Sample_expr_FastMNN, ident.1 =",i,", min.pct=0.25, only.pos = TRUE, pseudocount.use = 0.1, logfc_threshold = 1, assay = 'RNA')", sep=""))))
  eval(parse(text=(paste("cluster",i,".markers0.5 <- cluster",i,".markers0.5[order(cluster",i,".markers0.5$avg_log2FC, decreasing = TRUE),]", sep=""))))
  #eval(parse(text=(paste("write.table(cluster",i,".markers0.5, 'MarkerGenes_in_Cluster",i,"_res0.5.txt', sep='\t', quote=F, col.names=T, row.names=T)", sep=""))))
  print(paste("Evaluated the markers' significance of cluster n.",i))
  i<-i+1}

#### MONONUCLEAR PHAGOCYTES ####

Sample_expr_MP <- subset(Sample_expr, subset = Annotation_2 == 'MNPs')

Sample_expr_MP <- NormalizeData(Sample_expr_MP, normalization.method = "LogNormalize", scale.factor = 1e4, assay='RNA')
Sample_expr_MP <- FindVariableFeatures(Sample_expr_MP,selection.method = "vst", nfeatures = 3000)
Sample_expr_MP <- ScaleData(Sample_expr_MP, vars.to.regress = c("CC.Difference"))
Sample_expr_MP <- RunPCA(Sample_expr_MP)
Sample_expr_MP <- RunHarmony(Sample_expr_MP, group.by.vars = c('orig.ident'), dims.use = 1:30, theta=2, reduction.save = 'harmony')
Sample_expr_MP <- RunUMAP(Sample_expr_MP, reduction='harmony', dims = 1:20)
Sample_expr_MP <- FindNeighbors(Sample_expr_MP, reduction = 'harmony', dims = 1:20)
Sample_expr_MP <- FindClusters(Sample_expr_MP, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))

Idents(Sample_expr_MP) <- 'RNA_snn_res.1'
i <- 0
while(i<=15){
  eval(parse(text=(paste("cluster",i,".markers1 <- FindMarkers(Sample_expr_MP, ident.1 =",i,", min.pct=0.1, only.pos = TRUE, pseudocount.use = 0.1, logfc_threshold = 1, assay = 'RNA')", sep=""))))
  eval(parse(text=(paste("cluster",i,".markers1 <- cluster",i,".markers1[order(cluster",i,".markers1$avg_log2FC, decreasing = TRUE),]", sep=""))))
  eval(parse(text=(paste("write.table(cluster",i,".markers1, 'MarkerGenes_in_Cluster",i,"_res1.txt', sep='\t', quote=F, col.names=T, row.names=T)", sep=""))))
  print(paste("Evaluated the markers' significance of cluster n.",i))
  i<-i+1}
  
#### TUMOR-ASSOCIATED MACROPHAGES ####

Sample_expr_TAM <- subset(Sample_expr, subset = Annotation_3 == 'TAMs')

Sample_expr_TAM <- NormalizeData(Sample_expr_TAM, normalization.method = "LogNormalize", scale.factor = 1e4, assay='RNA')
Sample_expr_TAM <- FindVariableFeatures(Sample_expr_TAM,selection.method = "vst", nfeatures = 3000)
Sample_expr_TAM <- ScaleData(Sample_expr_TAM, vars.to.regress = c("CC.Difference"))
Sample_expr_TAM <- RunPCA(Sample_expr_TAM)
Sample_expr_TAM <- RunHarmony(Sample_expr_TAM, group.by.vars = c('orig.ident'), dims.use = 1:30, theta=2, reduction.save = 'harmony')
Sample_expr_TAM <- RunUMAP(Sample_expr_TAM, reduction='harmony', dims = 1:20)
Sample_expr_TAM <- FindNeighbors(Sample_expr_TAM, reduction = 'harmony', dims = 1:20)
Sample_expr_TAM <- FindClusters(Sample_expr_TAM, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))

Idents(Sample_expr_MP) <- 'RNA_snn_res.1'
i <- 0
while(i<=15){
  eval(parse(text=(paste("cluster",i,".markers1 <- FindMarkers(Sample_expr_MP, ident.1 =",i,", min.pct=0.1, only.pos = TRUE, pseudocount.use = 0.1, logfc_threshold = 1, assay = 'RNA')", sep=""))))
  eval(parse(text=(paste("cluster",i,".markers1 <- cluster",i,".markers1[order(cluster",i,".markers1$avg_log2FC, decreasing = TRUE),]", sep=""))))
  eval(parse(text=(paste("write.table(cluster",i,".markers1, 'MarkerGenes_in_Cluster",i,"_res1.txt', sep='\t', quote=F, col.names=T, row.names=T)", sep=""))))
  print(paste("Evaluated the markers' significance of cluster n.",i))
  i<-i+1}
  
TAM_annotation <- ifelse(Sample_expr_TAM$RNA_snn_res.0.4 == 0, 'Il1b_TAMs',
                        ifelse(Sample_expr_TAM$RNA_snn_res.0.4 == 1, 'Cxcl9_TAMs', 
                        ifelse(Sample_expr_TAM$RNA_snn_res.0.4 == 2, 'Spp1_TAMs',
                        ifelse(Sample_expr_TAM$RNA_snn_res.0.4 == 3, 'Folr2_TAMs',
                        ifelse(Sample_expr_TAM$RNA_snn_res.0.4 == 4, 'Clps_TAMs',
                        ifelse(Sample_expr_TAM$RNA_snn_res.0.4 == 5, 'Proliferating_TAMs', 'Marco_TAMs'))))))

Sample_expr_TAM$TAM_Annotation <- TAM_annotation

#### MONOCYTES AND MACROPHAGES ####

Sample_merge_MonoMacro <- NormalizeData(Sample_merge_MonoMacro, normalization.method = "LogNormalize", scale.factor = 1e4, assay='RNA')
Sample_merge_MonoMacro <- RunFastMNN(object.list = SplitObject(Sample_merge_MonoMacro, split.by = "orig.ident"))
palantir_so<-RunPalantirDiffusionMap(Sample_merge_MonoMacro_subset, reduction = "mnn", n_components = 20)
Sample_merge_MonoMacro_subset[["tsne_mnn"]] <-
    read.csv("tmp/tsne.csv", row.names = 1) %>%
    set_colnames(paste0("TSNE_FASTMNN_", 1:ncol(.))) %>%
    as.matrix() %>%
    CreateDimReducObject(key = "TSNEFASTMNN_", assay = DefaultAssay(Sample_merge_MonoMacro_subset))

## prepare annotations for velocity and Cellrank analysis

annotated_clusters <- as.data.frame(Sample_merge_MonoMacro_subset$Annotation)
colnames(annotated_clusters) <- 'clusters_refined'
write.csv(annotated_clusters, file='annotated_clusters.csv')

### run python notebook scripts for velocity analysis + Cellrank

## prepare data for optimal transport analysis

### run WOT scripts for optimal transport analysis

#### EPITHELIAL AND TUMOR CELLS ####

Sample_expr_Epithelial <- subset(Sample_expr, subset = Annotation_2 == 'Epithelial_cells')

Sample_expr_Epithelial <- NormalizeData(Sample_expr_Epithelial, normalization.method = "LogNormalize", scale.factor = 1e4, assay='RNA')
Sample_expr_Epithelial <- FindVariableFeatures(Sample_expr_Epithelial,selection.method = "vst", nfeatures = 3000)
Sample_expr_Epithelial <- ScaleData(Sample_expr_Epithelial, vars.to.regress = c("CC.Difference"))
Sample_expr_Epithelial <- RunPCA(Sample_expr_Epithelial)
Sample_expr_Epithelial <- RunHarmony(Sample_expr_Epithelial, group.by.vars = c('orig.ident'), dims.use = 1:20, theta=1, reduction.save = 'harmony')
Sample_expr_Epithelial <- RunUMAP(Sample_expr_Epithelial, reduction='harmony', dims = 1:20)
Sample_expr_Epithelial <- FindNeighbors(Sample_expr_Epithelial, reduction = 'harmony', dims = 1:20)
Sample_expr_Epithelial <- FindClusters(Sample_expr_Epithelial, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))

Idents(Sample_expr_Epithelial)<-Sample_expr_Epithelial$orig.ident

Epithelial_Day10_vs_Healthy <- FindMarkers(Sample_expr_Epithelial, ident.1 ='Tumor_d10',ident.2 ='Healthy', min.pct=0.1, only.pos = FALSE, pseudocount.use = 0.1, logfc_threshold = 0, assay = 'RNA')
#write.table(Epithelial_Day10_vs_Healthy[which(Epithelial_Day10_vs_Healthy$p_val_adj < 0.01 & abs(Epithelial_Day10_vs_Healthy$avg_log2FC) >= 1),], 'Epithelial_Day10_vs_Healthy.txt', sep='\t', quote=F, col.names=T, row.names=T)

Epithelial_Day20_vs_Healthy <- FindMarkers(Sample_expr_Epithelial, ident.1 ='Tumor_d20',ident.2 ='Healthy', min.pct=0.1, only.pos = FALSE, pseudocount.use = 0.1, logfc_threshold = 0, assay = 'RNA')
#write.table(Epithelial_Day20_vs_Healthy[which(Epithelial_Day20_vs_Healthy$p_val_adj < 0.01 & abs(Epithelial_Day20_vs_Healthy$avg_log2FC) >= 1),], 'Epithelial_Day20_vs_Healthy.txt', sep='\t', quote=F, col.names=T, row.names=T)

Epithelial_Day30_vs_Healthy <- FindMarkers(Sample_expr_Epithelial, ident.1 ='Tumor_d30',ident.2 ='Healthy', min.pct=0.1, only.pos = FALSE, pseudocount.use = 0.1, logfc_threshold = 0, assay = 'RNA')
#write.table(Epithelial_Day30_vs_Healthy[which(Epithelial_Day30_vs_Healthy$p_val_adj < 0.01 & abs(Epithelial_Day30_vs_Healthy$avg_log2FC) >= 1),], 'Epithelial_Day30_vs_Healthy.txt', sep='\t', quote=F, col.names=T, row.names=T)

hallmark_gene_sets = msigdbr(species = "mouse", category = "H")
mouse = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2018.archive.ensembl.org")

my_term_mouse=data.frame(hallmark_gene_sets$gs_name,hallmark_gene_sets$entrez_gene)

clusters_ordered_mouse = c(10,20,30)
for (i in clusters_ordered_mouse){
  eval(parse(text=(paste("tmp <- Epithelial_Day",i,"_vs_Healthy[,'avg_log2FC']", sep=""))))
  eval(parse(text=(paste("names(tmp) <- rownames(Epithelial_Day",i,"_vs_Healthy)", sep=""))))
  tmp <- tmp[which(tmp != "NA")]
  tmp <- sort(tmp, decreasing=TRUE)
  bioM_mouse=getBM(filters="mgi_symbol",values=names(tmp), attributes=c("entrezgene","mgi_symbol","description"),mart = mouse)
  gene_id<-as.character(unlist(mclapply(names(tmp), function(x) ifelse(x%in%bioM_mouse$mgi_symbol,bioM_mouse[which(bioM_mouse$mgi_symbol==x),1],"NA"),mc.cores = 4)))
  geneList <- tmp
  names(geneList) <- as.character(gene_id)
  geneList=geneList[which(names(geneList) != "NA")]
  eval(parse(text=(paste("GSEA_Epithelial_Day",i,"_vs_Healthy.mouse_HALLMARK <- GSEA(geneList, TERM2GENE = my_term_mouse, nPerm=100000, minGSSize= 15, maxGSSize=500, pvalueCutoff = 1,verbose = FALSE)", sep=""))))
  #eval(parse(text=(paste("write.table(GSEA_Epithelial_Day",i,"_vs_Healthy.mouse_HALLMARK@result,'Epithelial_Day",i,"_vs_Healthy_GSEA_HALLMARK.txt', sep='\t', quote=F, col.names=T, row.names=F)", sep=""))))
}


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


#### FULL DATASET ANALYSIS ####

Sample.merge <- NormalizeData(Sample.merge, normalization.method = "LogNormalize", scale.factor = 1e4, assay='RNA')
Sample.merge <- ScaleData(Sample.merge, vars.to.regress = c("CC.Difference"))
Sample.merge <- FindVariableFeatures(object = Sample.merge)
Sample.merge <- RunPCA(Sample.merge, pcs.compute=50)
Sample.merge <- RunFastMNN(object.list = SplitObject(Sample.merge, split.by = "orig.ident"))
Sample.merge <- RunUMAP(Sample.merge, reduction = "mnn", dims = 1:30)
Sample.merge <- FindNeighbors(Sample.merge, reduction = "mnn", dims = 1:30)
Sample.merge <- FindClusters(Sample.merge,  resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))


#### MONONUCLEAR PHAGOCYTES ####

Sample.merge_MP <- subset(Sample.merge, subset = RNA_snn_res.0.5 %in% c(0,9,12))
Sample.merge_MP <- ScaleData(Sample.merge_MP, vars.to.regress = c("CC.Difference"), features=rownames(Sample.merge_MP))
Sample.merge_MP <- FindVariableFeatures(Sample.merge_MP)
Sample.merge_MP <- RunPCA(Sample.merge_MP, pcs.compute=50)
Sample.merge_MP <- RunHarmony(Sample.merge_MP, "orig.ident", dims.use = 1:30, max.iter.harmony = 30)
Sample.merge_MP <- RunUMAP(Sample.merge_MP, reduction="harmony", dims = 1:ncol(Embeddings(Sample.merge_MP, "harmony")), reduction.name="umap", reduction.key="UMAPHARMONY_")
Sample.merge_MP <- FindNeighbors(Sample.merge_MP, reduction = "harmony", dims = 1:30)
Sample.merge_MP <- FindClusters(Sample.merge_MP,  resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))

#### TUMOR-ASSOCIATED MACROPHAGES ####

Sample.merge_TAM <- subset(Sample.merge_MP, subset = RNA_snn_res.1 %in% c(1,2,5,6,7,9,11))
Sample.merge_TAM <- ScaleData(Sample.merge_TAM, vars.to.regress = c("CC.Difference"), features=rownames(Sample.merge_TAM))
Sample.merge_TAM <- FindVariableFeatures(Sample.merge_TAM)
Sample.merge_TAM <- RunPCA(Sample.merge_TAM, pcs.compute=50)
Sample.merge_TAM <- RunHarmony(Sample.merge_TAM, "orig.ident", dims.use = 1:30, max.iter.harmony = 30, theta=3)
Sample.merge_TAM <- RunUMAP(Sample.merge_TAM, reduction="harmony", dims = 1:ncol(Embeddings(Sample.merge_MP, "harmony")), reduction.name="umap", reduction.key="UMAPHARMONY_")
Sample.merge_TAM <- FindNeighbors(Sample.merge_TAM, reduction = "harmony", dims = 1:30)
Sample.merge_TAM <- FindClusters(Sample.merge_TAM,  resolution = c(0.2,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))

# find markers for TAM subsets
Idents(Sample.merge_TAM) <- 'RNA_snn_res.0.36'
DEGs_TAMsubsets <- Reduce("rbind",lapply(unique(Sample.merge_TAM$RNA_snn_res.0.36), function(x) {
    Markers <- FindMarkers(Sample.merge_TAM, ident.1 = x, ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 1, pseudocount.use = 0.1)
    Markers <- Markers[which(Markers$p_val_adj < 0.01),]
    Markers$gene <- rownames(Markers)
    Markers$Cluster <- rep(paste("Cluster",x),nrow(Markers))
    return(Markers)
}))

# GSEA on GO BP
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2018.archive.ensembl.org")
bioM=getBM(filters="hgnc_symbol", values=rownames(Sample.merge_TAM), attributes=c("entrezgene","hgnc_symbol"), mart = mart)

GO_BP_GSEA_TAMsubsets <- Reduce("rbind",lapply(unique(Sample.merge_TAM$RNA_snn_res.0.36), function(x) {
    AllMarkers <- FindMarkers(Sample.merge_TAM, ident.1 = x, ident.2 = NULL, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0, pseudocount.use = 0.1)
    logFC = AllMarkers$avg_log2FC
    gene_id<-as.character(unlist(mclapply(rownames(AllMarkers), function(x) ifelse(x%in%bioM$hgnc_symbol,bioM[which(bioM$hgnc_symbol==x),1],"NA"),mc.cores = 4)))
    logFC <- logFC[!(is.na(gene_id))]
    names(logFC)= gene_id[!(is.na(gene_id))]
    logFC = sort(logFC, decreasing = TRUE)
    GSEA_bp <- gseGO(geneList = logFC, OrgDb = org.Hs.eg.db, ont= "BP",  minGSSize= 10, maxGSSize=500, pvalueCutoff = 1,verbose = FALSE)
    GSEA_bp<-data.frame(ID=GSEA_bp@result$ID, Description=GSEA_bp@result$Description, setSize=GSEA_bp@result$setSize, NES=GSEA_bp@result$NES,pvalue=GSEA_bp@result$pvalue, qvalues=GSEA_bp@result$qvalues)
    GSEA_bp$Cluster <- rep(paste("Cluster",x),nrow(GSEA_bp))
    GSEA_bp <- GSEA_bp[which(GSEA_bp$qvalues < 0.01),]
    return(GSEA_bp)
}))

#### HUMAN-MOUSE TAMs OVERLAP - GSEA 

# compute orthologous of expressed genes in TAMs subset
Idents(Sample_Macro_hg38) <- 'RNA_snn_res.0.36'
Idents(Sample_Macro_mm10) <- 'RNA_snn_res.0.4'
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "jul2018.archive.ensembl.org")
mouse = useMart("ensembl", dataset= "mmusculus_gene_ensembl", host = "jul2018.archive.ensembl.org")
genes_human = rownames(Sample_Macro_hg38)
genes_mouse = rownames(Sample_Macro_mm10)
genes_human_converted = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genes_human , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=F)
genes_mouse_converted = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes_mouse , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F)
unambiguous_mouse_genes = genes_mouse_converted %>% group_by(MGI.symbol) %>% count() %>% filter(n<2) %>% .$MGI.symbol
ambiguous_mouse_genes = genes_mouse_converted  %>% group_by(MGI.symbol) %>% count() %>% filter(n>=2) %>% .$MGI.symbol
geneinfo_ambiguous_solved = genes_mouse_converted %>% filter(MGI.symbol %in% ambiguous_mouse_genes) %>% filter(HGNC.symbol==toupper(MGI.symbol))
genes_mouse_converted = genes_mouse_converted %>% filter(MGI.symbol %in% unambiguous_mouse_genes) %>% bind_rows(geneinfo_ambiguous_solved)
rownames(genes_mouse_converted) =genes_mouse_converted[,1]	
genes_mouse_converted=genes_mouse_converted[!(duplicated(genes_mouse_converted[,1])),]
rownames(genes_mouse_converted) =genes_mouse_converted[,1]
expressed_genes_TAM_hg38 <- unique(unlist(lapply(unique(Sample_Macro_hg38$RNA_snn_res.0.36), function(x){
    cells <- rownames(Sample_Macro_hg38@meta.data[which(Sample_Macro_hg38@meta.data$RNA_snn_res.0.36 == x),])
    pct <- rowSums(Sample_Macro_hg38@assays$RNA@data[,cells]>0)/length(cells)
    return(names(pct[which(pct > 0.1)]))
})))
expressed_genes_TAM_mm10 <- unique(unlist(lapply(unique(Sample_Macro_mm10$RNA_snn_res.0.4), function(x){
    cells <- rownames(Sample_Macro_mm10@meta.data[which(Sample_Macro_mm10@meta.data$RNA_snn_res.0.4 == x),])
    pct <- rowSums(Sample_Macro_mm10@assays$RNA@data[,cells]>0)/length(cells)
    return(names(pct[which(pct > 0.1)]))
})))
tmp=genes_mouse_converted[which(genes_mouse_converted$MGI.symbol %in% expressed_genes_TAM_mm10),]
gene_to_mouse_common_expressed_TAM = tmp[which(tmp$HGNC.symbol %in% expressed_genes_TAM_hg38),]

term2gene <- Reduce("rbind", lapply(unique(Sample_Macro_hg38$RNA_snn_res.0.36), function(x){
    Markers <- FindMarkers(Sample_Macro_hg38, ident.1 = x, ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.8, pseudocount.use = 0.1)
	Markers <- Markers[which(Markers$p_val_adj < 0.01),]
	Markers_to_mm10 <- unlist(lapply(rownames(Markers), function(x) ifelse(x %in% gene_to_mouse_common_expressed_TAM$HGNC.symbol, gene_to_mouse_common_expressed_TAM[which(gene_to_mouse_common_expressed_TAM$HGNC.symbol == x),1], NA)))
	Markers_to_mm10 <- Markers_to_mm10[!(is.na(Markers_to_mm10))]
    term2gene = data.frame(id=rep(paste("Cluster_",x,sep=""),length(Markers_to_mm10)),gene=Markers_to_mm10)
    return(term2gene)
}))
term2name = data.frame(id=unique(term2gene[,1])[order(unique(term2gene[,1]))],Description=c("Hu_SPP1+","Hu_IL1B+","Hu_FOLR2+","Hu_HSP+","Hu_MT+","Hu_MKI67+"))

GSEA_human_to_Mouse_TAMs<-Reduce("rbind",lapply(unique(Sample_Macro_mm10$RNA_snn_res.0.4), function(x){
	AllMarkers <- FindMarkers(Sample_Macro_mm10, ident.1 = x, ident.2 = NULL, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0, pseudocount.use = 0.1)
    ranks=AllMarkers[order(AllMarkers$avg_log2FC,decreasing=T),"avg_log2FC"]
    names(ranks)=rownames(AllMarkers[order(AllMarkers$avg_log2FC,decreasing=T),])
    fgseaplot=GSEA(ranks, minGSSize = 10, maxGSSize = 500, eps = 1e-50, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = term2gene,TERM2NAME = term2name)
    fgseaplot@result[,7]=-log10(fgseaplot@result[,7])
    gsea <- data.frame(fgseaplot@result[,c(2,5,7)],rep(paste("Cluster",x,sep=""),nrow(fgseaplot@result)))
    return(gsea)
}))
colnames(GSEA_human_to_Mouse_TAMs) = c("Hu_TAMs","NES","log_padj","mouseTAMs_Cluster")

#### TUMOR CELLS IN NAIVE SAMPLES ####

cells_Naive<-rownames(Sample.merge@meta.data[which(Sample.merge@meta.data$orig.ident %in% c("LPDAC_30_tumor","PDAC_50_tumor","PDAC_55_tumor","PDAC_60_Tumor")),])
Sample.merge_Naive <- subset(Sample.merge, cells = cells_Naive)
Sample.merge_Naive <- NormalizeData(Sample.merge_Naive, normalization.method = "LogNormalize", scale.factor = 1e4)
#Sample.merge_Naive <- ScaleData(Sample.merge_Naive, vars.to.regress = c("CC.Difference"), features=rownames(Sample.merge_Naive))
#Sample.merge_Naive <- FindVariableFeatures(object = Sample.merge_Naive)
#Sample.merge_Naive <- RunPCA(Sample.merge_Naive, pcs.compute=50)
Sample.merge_Naive <- RunFastMNN(object.list = SplitObject(Sample.merge_Naive, split.by = "orig.ident"))
Sample.merge_Naive <- RunUMAP(Sample.merge_Naive, reduction = "mnn", dims = 1:30)
Sample.merge_Naive <- FindNeighbors(Sample.merge_Naive, reduction = "mnn", dims = 1:30)
Sample.merge_Naive <- FindClusters(Sample.merge_Naive,  resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))

cells_Tumor=rownames(Sample.merge_Naive@meta.data[which(Sample.merge_Naive@meta.data$RNA_snn_res.1 %in% c(1,3,13,15,5,8,17)),])
Sample_Tumor<-subset(Sample.merge_Naive, cells=cells_Tumor)
#Sample_Tumor <- ScaleData(Sample_Tumor, vars.to.regress = c("CC.Difference"), features=rownames(Sample_Tumor))
#Sample_Tumor <- FindVariableFeatures(object = Sample_Tumor)
#Sample_Tumor <- RunPCA(Sample_Tumor, pcs.compute=50)
Sample_Tumor <- RunFastMNN(object.list = SplitObject(Sample_Tumor, split.by = "orig.ident"))
Sample_Tumor <- RunUMAP(Sample_Tumor, reduction="mnn", dims = 1:30)
Sample_Tumor <- FindNeighbors(Sample_Tumor, reduction = "mnn", dims = 1:30)
Sample_Tumor <- FindClusters(Sample_Tumor,  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))

# find markers for Tumor cells subsets
Idents(Sample_Tumor) <- 'RNA_snn_res.0.3'
DEGs_TAMsubsets <- Reduce("rbind",lapply(unique(Sample_Tumor$RNA_snn_res.0.3), function(x) {
    Markers <- FindMarkers(Sample_Tumor, ident.1 = x, ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 1, pseudocount.use = 0.1)
    Markers <- Markers[which(Markers$p_val_adj < 0.01),]
    Markers$gene <- rownames(Markers)
    Markers$Cluster <- rep(paste("Cluster",x),nrow(Markers))
    return(Markers)
}))

# re-analysis of clusters enriched in T1RS+ cells
cells_T1RS <- rownames(Sample_Tumor@meta.data[which(Sample_Tumor@meta.data$RNA_snn_res.0.1 == 1),])

Sample_Tumor_T1RS <- subset(Sample_Tumor, cells=cells_T1RS)
#Sample_Tumor_T1RS <- ScaleData(Sample_Tumor_T1RS, vars.to.regress = c("CC.Difference"), features=rownames(Sample_Tumor_T1RS))
#Sample_Tumor_T1RS <- FindVariableFeatures(object = Sample_Tumor_T1RS)
#Sample_Tumor_T1RS_Naive <- RunPCA(Sample_Tumor_T1RS_Naive, pcs.compute=50)
Sample_Tumor_T1RS <- RunFastMNN(object.list = SplitObject(Sample_Tumor_T1RS, split.by = "orig.ident"))
Sample_Tumor_T1RS <- RunUMAP(Sample_Tumor_T1RS, reduction="mnn", dims = 1:20)
Sample_Tumor_T1RS <- FindNeighbors(Sample_Tumor_T1RS, reduction = "mnn", dims = 1:20)
Sample_Tumor_T1RS <- FindClusters(Sample_Tumor_T1RS,  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))

# pseudotime analysis with slingshot
Tumor_sl<- slingshot(Embeddings(Sample_Tumor_T1RS, "mnn")[,c(1:10)], clusterLabels = Sample_Tumor_T1RS$RNA_snn_res.0.2)
pt <- slingPseudotime(Tumor_sl)

pct <- rowSums(Sample_Tumor_T1RS@assays$RNA@counts > 0)/ncol(Sample_Tumor_T1RS@assays$RNA@counts)
expressedGenes <- names(pct[which(pct > 0.1)])
t <- na.omit(pt[,1])
y<- Sample_Tumor_T1RS@assays$RNA@scale.data[expressedGenes,names(t)]
corr <- apply(y,1,function(z){
    cor(t,z, method = "pearson")
})

## NICHENET analysis

# load nichenet networks
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

Idents(Sample.merge_Naive) <- 'Annotation_nichnet'
receiver="PDAC_cluster_T1RS_enriched"  # cells subjected to pseudotime analysis
sender="IL1B_TAM"  # IL1B+ TAMs

# reciver (PDAC cells subjected to pseudotime analysis) expressed genes 
DEG_PDAC<-FindMarkers(Sample.merge_Naive, ident.1="PDAC_cluster_T1RS_enriched", ident.2="other_PDAC_clusters", only.pos=TRUE, logfc.threshold=0.5, pseudocount.use=0.1)
expressed_genes_receiver = get_expressed_genes(receiver, Sample.merge_Naive, pct = 0.15) %>% .[. %in% rownames(DEG_PDAC)] 
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] 
expressed_receptors = intersect(receptors,expressed_genes_receiver)

# sender (IL1B+ TAMs) expressed genes 
DEG_IL1B<-FindMarkers(Sample.merge_Naive, ident.1="IL1B_TAM", ident.2="other_TAM", only.pos=TRUE, logfc.threshold=0.2, pseudocount.use=0.1)
list_expressed_genes_sender = sender %>% unique() %>% lapply(get_expressed_genes, Sample.merge_Naive, 0.15) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique() %>% .[. %in% rownames(DEG_IL1B)]
expressed_ligands = intersect(ligands,expressed_genes_sender)

# target genes: genes upregulated in T1RS+ PDAC cells (cluster at the end-point of pseudotime curve)
Idents(Sample_Tumor_T1RS) <- 'RNA_snn_res.0.4'
markers_Cluster2 <- FindMarkers(Sample_Tumor_T1RS, ident.1=2, logfc.threshold=1, min.pct=0.3, pseudocount.use=0.1, only.pos=TRUE)
geneset_oi = rownames(markers_Cluster2) 
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)] 

# MODEL
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
# targets
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
# receptors
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
# ligand pearson
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

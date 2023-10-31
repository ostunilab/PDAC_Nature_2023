library(data.table)
library(Matrix)
library(ggplot2)
library(future)
library(dplyr)
library(grid)
library(Seurat)
library(tidyr)
library(dendextend)
library(Giotto)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)


load('Spatial.filt.Robj')
images <- Images(Spatial.filt, assay = DefaultAssay(object = Spatial.filt))
image.use <- Spatial.filt[[images]]
coordinates <- GetTissueCoordinates(object = image.use)

# import proportions predicted by DestVI
proportions<-read.csv('CellProp_DestVI.csv', row.names=1)


#### Clustering and analysis of MonoMacro erniched spots 

SpatialPCs_MonoMacro<-read.csv('SpatialPCs_MonoMacro.csv', row.names=1)
Spatial.filt_MonoMacro<-subset(Spatial.filt, cells=rownames(SpatialPCs_MonoMacro))

List_simulations <- lapply(c(1:6), function(i) {
    sim<-read.csv(paste('Simulation_',i,'_MonoMacro.csv',sep=""), row.names=1)
    sim=sim[rownames(Spatial.filt_MonoMacro@meta.data),]
    return(sim)
})
simulation_mean<-Reduce("+",List_simulations)/length(List_simulations)

simulationMean <- CreateSeuratObject(t(simulation_mean), min.cells = 0,  project = "MonoMacro", min.features = 0)
simulationMean <- AddMetaData(simulationMean, Spatial.filt_MonoMacro@meta.data)
simulationMean <- ScaleData(simulationMean)
simulationMean[['SpatialPCA']] <- CreateDimReducObject(embeddings = as.matrix(SpatialPCs_MonoMacro[rownames(Spatial.filt_MonoMacro@meta.data),]), key="SpatialPCA_")
#simulationMean <- RunUMAP(simulationMean, reduction = "SpatialPCA", dims = 1:5)
#simulationMean <- FindNeighbors(simulationMean, reduction = "SpatialPCA", dims = 1:5)
#simulationMean <- FindClusters(simulationMean,  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))

# load TAMs markers from TABLE 5 (sheet TAM_markers_day30_MonoMacro)
Il1b_markers=Il1b_markers[Il1b_markers %in% rownames(simulationMean)]
Folr2_markers=Folr2_markers[Folr2_markers %in% rownames(simulationMean)]
Spp1_markers=Spp1_markers[Spp1_markers %in% rownames(simulationMean)]

# correlation with signatures gene expression ans Spatial PCs coordinates
pca_spatial<-Spatial.filt_MonoMacro@reductions$SpatialPCs@cell.embeddings[,1:5]
mean.exp_Il1b <- log(colMeans(as.matrix(expm1(simulationMean@assays$RNA@data[Il1b_markers,rownames(pca_spatial)])))+1)
mean.exp_Folr2 <- log(colMeans(as.matrix(expm1(simulationMean@assays$RNA@data[Folr2_markers,rownames(pca_spatial)])))+1)
mean.exp_Spp1 <- log(colMeans(as.matrix(expm1(simulationMean@assays$RNA@data[Spp1_markers,rownames(pca_spatial)])))+1)

dotplot_correlations=matrix(nrow=15, ncol=4)
dotplot_correlations=as.data.frame(dotplot_correlations)
colnames(dotplot_correlations) = c("TAM_subset","corr","p_value","PC")
for (i in 0:4) {
j=i+1
c<-cor.test(pca_spatial[order(pca_spatial[,j]),j],  mean.exp_Il1b[rownames(pca_spatial[order(pca_spatial[,j]),])], method=c("pearson"))
dotplot_correlations[i*3+1,1]="Il1b"
dotplot_correlations[i*3+1,2]=c$estimate
dotplot_correlations[i*3+1,3]=c$p.value
eval(parse(text=paste("dotplot_correlations[",i*3+1,",4]=\'PC_",j,"\'",sep="")))

c<-cor.test(pca_spatial[order(pca_spatial[,j]),j],  mean.exp_Folr2[rownames(pca_spatial[order(pca_spatial[,j]),])], method=c("pearson"))
dotplot_correlations[i*3+2,1]="Folr2"
dotplot_correlations[i*3+2,2]=c$estimate
dotplot_correlations[i*3+2,3]=c$p.value
eval(parse(text=paste("dotplot_correlations[",i*3+2,",4]=\'PC_",j,"\'",sep="")))

c<-cor.test(pca_spatial[order(pca_spatial[,j]),j],  mean.exp_Spp1[rownames(pca_spatial[order(pca_spatial[,j]),])], method=c("pearson"))
dotplot_correlations[i*3+3,1]="Spp1"
dotplot_correlations[i*3+3,2]=c$estimate
dotplot_correlations[i*3+3,3]=c$p.value
eval(parse(text=paste("dotplot_correlations[",i*3+3,",4]=\'PC_",j,"\'",sep="")))
}
dotplot_correlations[,3]=-log10(dotplot_correlations[,3])
dotplot_correlations[,3]= MinMax(dotplot_correlations[,3], min = 0, max = 30)
dotplot_correlations[,1]=factor(dotplot_correlations[,1], levels=c("Spp1","Folr2","Il1b"))

# load TAMs markers from TABLE 5 (sheet TAM_markers_day30_MonoMacro)
TAMs<-list(Il1b_markers,Folr2_markers,Spp1_markers)
TAMsSign<-makeSignMatrixPAGE(sign_list = TAMs, sign_names=c('Il1b_markers','Folr2_markers','Spp1_markers'))

# Giotto signature enrichemnt analysis (PAGE) for TAMs markers (same analysis for other lists of genes)
giotto.obj = createGiottoObject(raw_exprs = Spatial.filt@assays$Spatial@counts, spatial_locs = coordinates)
giotto.obj <- normalizeGiotto(gobject = giotto.obj, scalefactor = 6000, verbose = T)
TAMs<-list(Il1b_markers,Folr2_markers,Spp1_markers)
TAMsSign<-makeSignMatrixPAGE(sign_list = TAMs, sign_names=c('Il1b_markers','Folr2_markers','Spp1_markers'))
giotto.obj <- runPAGEEnrich(gobject = giotto.obj, p_value = TRUE, sign_matrix = TAMsSign, output_enrichment='original', min_overlap_genes=5, include_depletion=F, expression_values='normalized') # with pvalues; -log10(p) is returned
# highlight only MonoMacro enriched spots with p < 0.001 
enrichment_PAGE=as.data.frame(giotto.obj@spatial_enrichment$PAGE) 
rownames(enrichment_PAGE)=enrichment_PAGE[,1]
for (i in 2:ncol(enrichment_PAGE)){         
enrichment_PAGE[,i] = unlist(lapply(enrichment_PAGE[,i], function(x) ifelse(x >= 3, 1, 0)))
}
enrichment_PAGE=as.data.frame(enrichment_PAGE) 
for (i in 1:nrow(enrichment_PAGE)){
	if(! (enrichment_PAGE$cell_ID[i] %in% colnames(Spatial.filt_MonoMacro))) {
		enrichment_PAGE[i,2:ncol(enrichment_PAGE)] = rep(0,ncol(enrichment_PAGE)-1)
	}
}

## GO_BP enrichment analysis 
Idents(SpatialA1FilteredFil) <- 'SCT_snn_res.0.3'
Cluster_Il1bvsAll_spatial_seurat <- FindMarkers(SpatialA1FilteredFil, ident.1 = 4, ident.2 = NULL, only.pos = FALSE, min.pct = 0.1, pseudocount.use = 0.1,logfc.threshold = 0)

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2018.archive.ensembl.org")
genes <- rownames(Cluster_Il1bvsAll_spatial_seurat)
logFC <- data.frame(Cluster_Il1bvsAll_spatial_seurat$avg_log2FC)
bioM <- getBM(filters="mgi_symbol",values=genes, attributes=c("entrezgene","mgi_symbol","description"),mart = mart)
gene_id <- as.character(unlist(mclapply(genes, function(x) ifelse(x%in%bioM$mgi_symbol,bioM[which(bioM$mgi_symbol==x),1],"NA"),mc.cores = 4)))
logFC <- logFC[!(is.na(gene_id))]
names(logFC) <- gene_id[!(is.na(gene_id))]
logFC <- sort(logFC, decreasing = TRUE)
GSEA_Il1b <- gseGO(geneList = logFC, OrgDb = org.Mm.eg.db, ont= "BP", minGSSize= 15, maxGSSize=500, pvalueCutoff = 1,verbose = FALSE)
geneSets_list <- GSEA_Il1b@geneSets
GSEA_Il1b <- data.frame(ID=GSEA_Il1b@result$ID, Description=GSEA_Il1b@result$Description, setSize=GSEA_Il1b@result$setSize, NES=GSEA_Il1b@result$NES,pvalue=GSEA_Il1b@result$pvalue, qvalues=GSEA_Il1b@result$qvalues)
GSEA_Il1b <- GSEA_Il1b[which(GSEA_Il1b$qvalues < 0.01),]
# extract gene names for Bio processes terms
geneName_list_go_il1b <- sapply(GSEA_Il1b$ID, function(i){
bioM=getBM(filters="entrezgene",values=geneSets_list[[i]], attributes=c("entrezgene","mgi_symbol","description"),mart = mart)
gene_symbol<-as.character(unlist(mclapply(geneSets_list[[i]], function(x) ifelse(x%in%bioM$entrezgene,bioM[which(bioM$entrezgene==x),2],NA),mc.cores = 4)))
}, simplify = FALSE, USE.NAMES = TRUE)

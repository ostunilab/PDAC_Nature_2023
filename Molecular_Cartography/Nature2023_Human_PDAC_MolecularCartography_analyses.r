library(Matrix)
library(graphics)
library(parallel)
library(Seurat)
library(dplyr)
library(spdep)
library(dbscan)
library(ggplot2)
library(scales)

 LoadResolve <- function (data, fov, assay = "Resolve") {
    segs <- CreateSegmentation(data$segmentations)
    cents <- CreateCentroids(data$centroids)
    segmentations.data <- list(centroids = cents, segmentation = segs)
    coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", "centroids"), molecules = data$microns, assay = assay)
    obj <- CreateSeuratObject(counts = data$transcripts, assay = assay)
    coords <- subset(x = coords, cells = intersect(x = Cells(x = coords[["segmentation"]]), y = Cells(x = obj)))
    obj[[fov]] <- coords
    return(obj)
}

##### UPLOAD baysor output and create a Seurat object
setwd('GSM7655264/outs')
#setwd('GSM7655265/outs')
#setwd('GSM7655264/outs')


data <- vector(mode = 'list', length = 4)
names(x = data) = c("transcripts","centroids","segmentations","microns")

###### centroids
cell_centroids<-read.csv('segmentation_cell_stats.csv')
cell_centroids$cell = paste("Cell_",cell_centroids$cell,sep="")
data[['centroids']] = cell_centroids[,c("x","y","cell")]

###### microns
baysor_out<-read.csv('segmentation.csv')
table(baysor_out$is_noise)
baysor_out = baysor_out[ baysor_out$is_noise == "false",]
baysor_out$cell = paste("Cell_",baysor_out$cell,sep="")
data[['microns']] = baysor_out[,c("x","y","gene")]

###### transcripts
genes<-unique(baysor_out$gene)
genes<-genes[order(genes)]
cells_id <- cell_centroids$cell
counts<-mclapply(cells_id, function(x) {
	as.vector(table(factor(baysor_out[baysor_out[,"cell"] == x,"gene"], levels=genes)))
})
counts <- do.call(cbind.data.frame, counts)
rownames(counts) = genes
colnames(counts) = cells_id
counts <- Matrix(as.matrix(counts), sparse = TRUE)
data[['transcripts']] = counts


####### Compute segmentation 
cells_id <- cell_centroids$cell
id_edges_segmentation <- unlist(mclapply(cells_id, function(x) {
	test<-chull(baysor_out[baysor_out[,"cell"] == x,c("x","y")])
	c(rownames(baysor_out[baysor_out[,"cell"] == x,][test,]),rownames(baysor_out[baysor_out[,"cell"] == x,][test,])[1])
	})
	)
segmentation = baysor_out[id_edges_segmentation,c("x","y","cell")]
rownames(segmentation) = c(1:nrow(segmentation))
data[['segmentations']] = segmentation

resolve_B2_1<- LoadResolve(data, "GSM7655264")
resolve_B2_1<-RenameCells(resolve_B2_1, add.cell.id= "B2_1")
DefaultBoundary(resolve_B2_1[["B2_1"]]) <- "segmentation"
Sample_B2_1 <- CreateSeuratObject(resolve_B2_1@assays$Resolve@counts, min.cells = 10,  project = "B2_1", min.features = 4)

#resolve_C2_1<- LoadResolve(data, "GSM7655265")
#resolve_C2_1<-RenameCells(resolve_C2_1, add.cell.id= "C2_1")
#DefaultBoundary(resolve_C2_1[["C2_1"]]) <- "segmentation"
#Sample_C2_1 <- CreateSeuratObject(resolve_C2_1@assays$Resolve@counts, min.cells = 10,  project = "C2_1", min.features = 4)

#resolve_D2_1<- LoadResolve(data, "GSM7655266")
#resolve_D2_1<-RenameCells(resolve_D2_1, add.cell.id= "D2_1")
#DefaultBoundary(resolve_D2_1[["D2_1"]]) <- "segmentation"
#Sample_D2_1 <- CreateSeuratObject(resolve_D2_1@assays$Resolve@counts, min.cells = 10,  project = "D2_1", min.features = 4)

############


# Data analysis
Sample.merge<- merge(Sample_B2_1, y = c(Sample_C2_1,Sample_D2_1))
Sample.merge <- subset(Sample.merge, subset =  nCount_RNA >= 10 & nFeature_RNA  <= 25)
Sample.merge <- SCTransform(Sample.merge, assay = "RNA", clip.range = c(-10, 10), )
Sample.merge <- RunPCA(Sample.merge, npcs = 30, features = rownames(Sample.merge))
Sample.merge <- RunUMAP(Sample.merge, dims = 1:20)
Sample.merge <- FindNeighbors(Sample.merge, reduction = "pca", dims = 1:20)
Sample.merge <- FindClusters(Sample.merge, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

Idents(Sample.merge) <- 'SCT_snn_res.0.5'
DEGs_res0.5 <- Reduce("rbind",lapply(unique(Sample.merge$SCT_snn_res.0.5), function(x) {
    Markers <- FindMarkers(Sample.merge, ident.1 = x, ident.2 = NULL, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 1, pseudocount.use = 0.3)
    Markers$gene <- rownames(Markers)
    Markers$Cluster <- rep(paste("Cluster",x),nrow(Markers))
    return(Markers)
}))

resolve_B2_1 <- subset(Sample.merge, cells= rownames(Sample.merge@meta.data[which(Sample.merge@meta.data$orig.ident == "B2"),]))
resolve_C2_1 <- subset(Sample.merge, cells= rownames(Sample.merge@meta.data[which(Sample.merge@meta.data$orig.ident == "C2"),]))
resolve_D2_1 <- subset(Sample.merge, cells= rownames(Sample.merge@meta.data[which(Sample.merge@meta.data$orig.ident == "D2"),]))

# Neighbourhood analysis
resolve.obj=list(B2_1=resolve_B2_1,C2_1=resolve_C2_1, D2_1=resolve_D2_1)

fraction_NN=Reduce("+",lapply(resolve.obj, function(x) {
    xy_cells<-GetTissueCoordinates(x[[names(x)[2]]][["centroids"]])
    rownames(xy_cells)=xy_cells[,"cell"]

    cells_sel<-rownames(x@meta.data[which(x@meta.data$clusters == 16),])

    dim_NN<-c(length=40)    
    for (nNeighbours in 1:40){   
    knn_spatial <- dbscan::kNN(x = xy_cells[, c("x", "y")] %>% as.matrix(), k = nNeighbours)
    knn_spatial.norm <- data.frame(from = rep(1:nrow(knn_spatial$id), nNeighbours),
                                 to = as.vector(knn_spatial$id),
                                 weight = 1/(1 + as.vector(knn_spatial$dist)),
                                 distance = as.vector(knn_spatial$dist))
    knn_spatial.norm$from = rownames(xy_cells)[knn_spatial.norm$from]                                 
    knn_spatial.norm$to= rownames(xy_cells)[knn_spatial.norm$to]         

    nn <- unique(knn_spatial.norm[which(knn_spatial.norm$from %in% cells_sel & knn_spatial.norm$distance < 400),"to"])
    dim_NN<-rbind(dim_NN,table(factor(x@meta.data[nn,"clusters"], levels=c(0:18)) ))
    }
    dim_NN = dim_NN[-1,]
    return(dim_NN)
}))

for (i in 1:nrow(fraction_NN)){
    fraction_NN[i,] = fraction_NN[i,]/table(factor(Sample.merge@meta.data[,"SCT_snn_res.0.5"], levels=c(0:18)))
}


# Neighbourhood enrichment
resolve.obj=list(B2_1=resolve_B2_1,C2_1=resolve_C2_1, D2_1=resolve_D2_1)

rand_out_lists <- mclapply(resolve.obj, function(x) {
    xy_cells<-GetTissueCoordinates(x[[names(x)[2]]][["centroids"]])
    rownames(xy_cells)=xy_cells[,"cell"]
    nNeighbours =40
    knn_spatial <- dbscan::kNN(x = xy_cells[, c("x", "y")] %>% as.matrix(), k = nNeighbours)
    knn_spatial.norm <- data.frame(from = rep(1:nrow(knn_spatial$id), nNeighbours), to = as.vector(knn_spatial$id),distance = as.vector(knn_spatial$dist))
    knn_spatial.norm$from = rownames(xy_cells)[knn_spatial.norm$from]                                 
    knn_spatial.norm$to= rownames(xy_cells)[knn_spatial.norm$to]   
    randomization<-matrix(nrow=ncol(x), ncol=1000)
    randomization<-apply(randomization,2, function(y) {
        y = sample(as.vector(x@meta.data$clusters))
    })
    rownames(randomization) = rownames(x@meta.data)

    rand_out=apply(randomization, 2, function(y) {
        y=factor(y,levels=levels(x@meta.data$clusters))
        cells_sel<-list()
        for (i in levels(y)){
            cells_sel[[i]] = names(y[y == i])
        }

        dim_NN=lapply(cells_sel, function(z) {
            nn <- unique(knn_spatial.norm[which(knn_spatial.norm$from %in% z & knn_spatial.norm$distance < 400),"to"])
            table(factor(y[nn], levels=levels(x@meta.data$clusters)))
        })
        return(dim_NN)
    })

    rand_out_list<-list()
    for (i in levels(x@meta.data$clusters)){
        rand_out_list[[i]] = Reduce('rbind', lapply(rand_out, function(z) z[[i]]))
    }
    return(rand_out_list)
},mc.cores = 16)

randomized_NN<-list()
for (i in levels(resolve.obj[[1]]@meta.data$clusters)){
    randomized_NN[[i]] = Reduce('+', lapply(rand_out_lists, function(z) z[[i]]))
}

real_NN=Reduce("+",lapply(resolve.obj, function(x) {
    xy_cells<-GetTissueCoordinates(x[[names(x)[2]]][["centroids"]])
    rownames(xy_cells)=xy_cells[,"cell"]
    nNeighbours =40
    knn_spatial <- dbscan::kNN(x = xy_cells[, c("x", "y")] %>% as.matrix(), k = nNeighbours)
    knn_spatial.norm <- data.frame(from = rep(1:nrow(knn_spatial$id), nNeighbours),
                                 to = as.vector(knn_spatial$id),
                                 distance = as.vector(knn_spatial$dist))
    knn_spatial.norm$from = rownames(xy_cells)[knn_spatial.norm$from]                                 
    knn_spatial.norm$to= rownames(xy_cells)[knn_spatial.norm$to]   

    cells_sel<-list()
    for (i in levels(x@meta.data$clusters)){
        cells_sel[[i]] = rownames(x@meta.data[which(x@meta.data$clusters == i),])
    }
        
    dim_NN=Reduce('rbind',lapply(cells_sel, function(z) {
        nn <- unique(knn_spatial.norm[which(knn_spatial.norm$from %in% z & knn_spatial.norm$distance < 400),"to"])
        dim_nn <- table(factor(x@meta.data[nn,"clusters"], levels=levels(x@meta.data$clusters)))
        return(dim_nn)
    }))
    return(dim_NN)
}))
rownames(real_NN) = levels(resolve.obj[[1]]@meta.data$clusters)

z_scores <- Reduce('rbind',lapply(names(randomized_NN), function(x){
    z_scores<-c()
    for (i in 1:ncol(randomized_NN[[x]])){
        z_scores<-c(z_scores,(real_NN[x,i] - mean(randomized_NN[[x]][,i]))/ sd(randomized_NN[[x]][,i]))
    }
    return(z_scores)
}))



library(data.table)
library(RColorBrewer)
library(anndata)
library(Matrix)
library(ggplot2)
library(future)
library(dplyr)
library(scales)
library(grid)
library(Seurat)
library(SeuratData)
library(tidyr)
library(dendextend)
library(cowplot)
library(patchwork)


####### HOW TO UPLOAD HIGH RES IMAGE (tissue_hires_image.png)
library(jsonlite)
library(png)
Read10X_Image <- function(image.dir, image.name = "tissue_hires_image.png", filter.matrix = TRUE, ...) {
  image <- readPNG(source = file.path(image.dir, image.name))
  scale.factors <- fromJSON(txt = file.path(image.dir, 'scalefactors_json.json'))
  tissue.positions.path <- Sys.glob(paths = file.path(image.dir, 'tissue_positions*'))
  tissue.positions <- read.csv(
    file = tissue.positions.path[1],
    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
    header = ifelse(
      test = basename(tissue.positions.path[1]) == "tissue_positions.csv",
      yes = TRUE,
      no = FALSE
    ),
    as.is = TRUE,
    row.names = 1
  )
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_hires_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  return(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_hires_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ))
}
#############

image<-Read10X_Image('GSM6727528/outs/spatial/')
Spatial<-Load10X_Spatial('GSM6727528/outs/', image=image)

#image<-Read10X_Image('GSM6727529/outs/spatial/')
#Spatial<-Load10X_Spatial('GSM6727529/outs/', image=image)


Spatial.filt <- subset(Spatial, subset =  nFeature_Spatial  > 100)
Spatial.filt <- SCTransform(Spatial.filt, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE )
Spatial.filt <- RunPCA(Spatial.filt, assay = "SCT", verbose = FALSE, pcs.compute=50)
Spatial.filt <- FindNeighbors(Spatial.filt, reduction = "pca", dims = 1:20)
Spatial.filt <- FindClusters(Spatial.filt, verbose = FALSE, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))
Spatial.filt <- RunUMAP(Spatial.filt, reduction = "pca", dims = 1:20)

write.csv(rownames(Spatial.filt@meta.data),"SelectedSpots.csv")
save(Spatial.filt,file="Spatial.filt.Robj")
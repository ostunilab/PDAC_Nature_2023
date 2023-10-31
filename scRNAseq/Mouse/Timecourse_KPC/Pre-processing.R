library(Seurat)
library(scDblFinder)

#### PRE-PROCESSING ####

Sample.Healthy <- Read10X('GSM6727561/filtered_feature_bc_matrix/')
Sample.Healthy <- CreateSeuratObject(Sample.Healthy, min.cells = 3,  project ="Healthy")

Sample.d10.Tumor <- Read10X('GSM6727558/filtered_feature_bc_matrix/')
Sample.d10.Tumor <- CreateSeuratObject(Sample.d10.Tumor, min.cells = 3,  project ="Tumor_d10")

Sample.d20.Tumor <- Read10X('GSM6727559/filtered_feature_bc_matrix/')
Sample.d20.Tumor <- CreateSeuratObject(Sample.d20.Tumor, min.cells = 3,  project ="Tumor_d20")

Sample.d30.Tumor <- Read10X('GSM6727560/filtered_feature_bc_matrix/')
Sample.d30.Tumor <- CreateSeuratObject(Sample.d30.Tumor, min.cells = 3,  project ="Tumor_d30")

Sample_expr <- merge(Sample.d10.Tumor, y = c(Sample.d20.Tumor, Sample.Healthy, Sample.d30.Tumor), 
                         add.cell.ids = c('Tumor_d10','Tumor_d20','Healthy','Tumor_d30'))
 
Sample_expr[['percent.mt']] <- PercentageFeatureSet(Sample_expr, pattern = '^mt-')
Sample_expr[['percent.ribo']] <- PercentageFeatureSet(Sample_expr, pattern = '^Rp[sl]')
s.genes <- readLines('ccgenes_mm_Sphase.txt')
g2m.genes <- readLines('ccgenes_mm_G2Mphase.txt')
Sample_expr <- CellCycleScoring(Sample_expr, g2m.features=g2m.genes[g2m.genes %in% rownames(Sample_expr@assays$RNA@data)], s.features=s.genes[s.genes %in% rownames(Sample_expr@assays$RNA@data)], set.ident = FALSE)
Sample_expr@meta.data$CC.Difference <- Sample_expr@meta.data$S.Score - Sample_expr@meta.data$G2M.Score

Sample_expr <- subset(Sample_expr, subset = percent.mt < 25 & nCount_RNA > 1000 & nFeature_RNA > 200)

for (i in c('Healthy','Tumor_d10','Tumor_d20','Tumor_d30')){
	sub <- subset(Sample_expr, subset = orig.ident == i)
	eval(parse(text=paste("sceDblF_",i," <- scDblFinder(sub@assays$RNA@counts, dbr = 0.07)",sep="")))
	eval(parse(text=paste("score.",i," <- sceDblF_",i,"@colData@listData[['scDblFinder.score']]",sep="")))
	eval(parse(text=paste("names(score.",i,") <- rownames(sceDblF_",i,"@colData)",sep="")))
}

doublets.info <- rbind(sceDblF_Tumor_d10@colData,sceDblF_Tumor_d20@colData,sceDblF_Healthy@colData,sceDblF_Tumor_d30@colData)
Sample_expr$is.doublet <- doublets.info$scDblFinder.class

Sample_expr <- subset(Sample_expr, subset = is.doublet == 'singlet')
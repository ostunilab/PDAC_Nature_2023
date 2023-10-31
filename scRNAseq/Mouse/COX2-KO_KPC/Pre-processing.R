library(Seurat)
library(scDblFinder)
set.seed(123)

#### PRE-PROCESSING ####

Sample.d7.WT <- Read10X('GSM6727566/filtered_feature_bc_matrix/')
Sample.d7.WT <- CreateSeuratObject(Sample.d7.WT, min.cells = 3,  project ="WT")

Sample.d7.KO <- Read10X('GSM6727567/filtered_feature_bc_matrix/')
Sample.d7.KO <- CreateSeuratObject(Sample.d7.KO, min.cells = 3,  project ="KO")

Sample_expr <- merge(Sample.d7.WT, y = c(Sample.d7.KO), add.cell.ids = c('WT','KO'))
 
Sample_expr[['percent.mt']] <- PercentageFeatureSet(Sample_expr, pattern = '^mt-')
Sample_expr[['percent.ribo']] <- PercentageFeatureSet(Sample_expr, pattern = '^Rp[sl]')
s.genes <- readLines('ccgenes_mm_Sphase.txt')
g2m.genes <- readLines('ccgenes_mm_G2Mphase.txt')
Sample_expr <- CellCycleScoring(Sample_expr, g2m.features=g2m.genes[g2m.genes %in% rownames(Sample_expr@assays$RNA@data)], s.features=s.genes[s.genes %in% rownames(Sample_expr@assays$RNA@data)], set.ident = FALSE)
Sample_expr@meta.data$CC.Difference <- Sample_expr@meta.data$S.Score - Sample_expr@meta.data$G2M.Score

Sample_expr <- subset(Sample_expr, subset = percent.mt < 25 & nFeature_RNA > 200)

for (i in c('WT','KO')){
	sub <- subset(Sample_expr, subset = orig.ident == i)
	eval(parse(text=paste("sceDblF_",i," <- scDblFinder(sub@assays$RNA@counts, dbr = 0.05)",sep="")))
	eval(parse(text=paste("score.",i," <- sceDblF_",i,"@colData@listData[['scDblFinder.score']]",sep="")))
	eval(parse(text=paste("names(score.",i,") <- rownames(sceDblF_",i,"@colData)",sep="")))
}

doublets.info <- rbind(sceDblF_WT@colData,sceDblF_KO@colData)
Sample_expr$is.doublet <- doublets.info$scDblFinder.class

Sample_expr <- subset(Sample_expr, subset = is.doublet == 'singlet')
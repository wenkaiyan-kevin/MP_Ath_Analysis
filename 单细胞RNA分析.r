#R version 3.6.1
#Seurat version 3.1.5
library(Seurat)
library(cowplot)

#First Read the three replicate of Ath root from WT(Ryu et al., 2019)
#############
WT1.data <-Read10X(data.dir = "../2_RNA_cellranger_previous/scRNA-seq_rep1/outs/filtered_feature_bc_matrix")
WT1 <- CreateSeuratObject(counts = WT1.data, project = "WT1", min.cells = 5)
WT1 <- RenameCells(object = WT1, new.names = gsub("-1","-WT1",colnames(WT1)))
WT1[["percent.mt"]] <- PercentageFeatureSet(WT1, pattern = "^ATM")

doublets_WT1 <- read.table("../2_RNA_cellranger_previous/filter_doublet/scRNA-seq_rep1.doublet.txt")
WT1 <- subset(WT1, cells=WhichCells(WT1, doublets_WT1[,1]), invert=TRUE)
WT1 <- subset(WT1, cells = Cells(WT1), features = rownames(WT1)[grep('AT[1-5]G', rownames(WT1))])
WT1 <- subset(WT1, subset = nFeature_RNA  > 200)

WT1$protocol <- "WT1"
WT1$tech <- "sc-RNA"
WT1 <- NormalizeData(WT1)
WT1 <- FindVariableFeatures(WT1, selection.method = "vst", nfeatures = 2000)
################
WT2.data <-Read10X(data.dir = "../2_RNA_cellranger_previous/scRNA-seq_rep2/outs/filtered_feature_bc_matrix")
WT2 <- CreateSeuratObject(counts = WT2.data, project = "WT2", min.cells = 5)
WT2 <- RenameCells(object = WT2, new.names = gsub("-1","-WT2",colnames(WT2)))
WT2[["percent.mt"]] <- PercentageFeatureSet(WT2, pattern = "^ATM")

WT2 <- subset(WT2, cells = Cells(WT2), features = rownames(WT2)[grep('AT[1-5]G', rownames(WT2))])
WT2 <- subset(WT2, subset = nFeature_RNA  > 200)

WT2$protocol <- "WT2"
WT2$tech <- "sc-RNA"
WT2 <- NormalizeData(WT2)
WT2 <- FindVariableFeatures(WT2, selection.method = "vst", nfeatures = 2000)
################
WT3.data <-Read10X(data.dir = "../2_RNA_cellranger_previous/scRNA-seq_rep3/outs/filtered_feature_bc_matrix")
WT3 <- CreateSeuratObject(counts = WT3.data, project = "WT3", min.cells = 5)
WT3 <- RenameCells(object = WT3, new.names = gsub("-1","-WT3",colnames(WT3)))
WT3[["percent.mt"]] <- PercentageFeatureSet(WT3, pattern = "^ATM")

WT3 <- subset(WT3, cells = Cells(WT3), features = rownames(WT3)[grep('AT[1-5]G', rownames(WT3))])
WT3 <- subset(WT3, subset = nFeature_RNA  > 200)

WT3$protocol <- "WT3"
WT3$tech <- "sc-RNA"
WT3 <- NormalizeData(WT3)
WT3 <- FindVariableFeatures(WT3, selection.method = "vst", nfeatures = 2000)

#Second Read the Five replicate of Ath root (Andrew Farmer et al., 2021)
#############
rep1.data <-Read10X(data.dir = "../1_RNA_cellranger_current/sNucRNA-seq_rep1/outs/filtered_feature_bc_matrix")
rep1 <- CreateSeuratObject(counts = rep1.data, project = "rep1", min.cells = 5)
rep1 <- RenameCells(object = rep1, new.names = gsub("-1","-rep1",colnames(rep1)))
rep1[["percent.mt"]] <- PercentageFeatureSet(rep1, pattern = "^ATM")

rep1 <- subset(rep1, cells = Cells(rep1), features = rownames(rep1)[grep('AT[1-5]G', rownames(rep1))])
rep1 <- subset(rep1, subset = nFeature_RNA  > 200)

rep1$protocol <- "rep1"
rep1$tech <- "sn-RNA"
rep1 <- NormalizeData(rep1)
rep1 <- FindVariableFeatures(rep1, selection.method = "vst", nfeatures = 2000)
#############
rep2.data <-Read10X(data.dir = "../1_RNA_cellranger_current/sNucRNA-seq_rep2/outs/filtered_feature_bc_matrix")
rep2 <- CreateSeuratObject(counts = rep2.data, project = "rep2", min.cells = 5)
rep2 <- RenameCells(object = rep2, new.names = gsub("-1","-rep2",colnames(rep2)))
rep2[["percent.mt"]] <- PercentageFeatureSet(rep2, pattern = "^ATM")

rep2 <- subset(rep2, cells = Cells(rep2), features = rownames(rep2)[grep('AT[1-5]G', rownames(rep2))])
rep2 <- subset(rep2, subset = nFeature_RNA  > 200)

rep2$protocol <- "rep2"
rep2$tech <- "sn-RNA"
rep2 <- NormalizeData(rep2)
rep2 <- FindVariableFeatures(rep2, selection.method = "vst", nfeatures = 2000)
#############
rep3.data <-Read10X(data.dir = "../1_RNA_cellranger_current/sNucRNA-seq_rep3/outs/filtered_feature_bc_matrix")
rep3 <- CreateSeuratObject(counts = rep3.data, project = "rep3", min.cells = 5)
rep3 <- RenameCells(object = rep3, new.names = gsub("-1","-rep3",colnames(rep3)))
rep3[["percent.mt"]] <- PercentageFeatureSet(rep3, pattern = "^ATM")

rep3 <- subset(rep3, cells = Cells(rep3), features = rownames(rep3)[grep('AT[1-5]G', rownames(rep3))])
rep3 <- subset(rep3, subset = nFeature_RNA  > 200)

rep3$protocol <- "rep3"
rep3$tech <- "sn-RNA"
rep3 <- NormalizeData(rep3)
rep3 <- FindVariableFeatures(rep3, selection.method = "vst", nfeatures = 2000)
#############
rep4.data <-Read10X(data.dir = "../1_RNA_cellranger_current/sNucRNA-seq_rep4/outs/filtered_feature_bc_matrix")
rep4 <- CreateSeuratObject(counts = rep4.data, project = "rep4", min.cells = 5)
rep4 <- RenameCells(object = rep4, new.names = gsub("-1","-rep4",colnames(rep4)))
rep4[["percent.mt"]] <- PercentageFeatureSet(rep4, pattern = "^ATM")

rep4 <- subset(rep4, cells = Cells(rep4), features = rownames(rep4)[grep('AT[1-5]G', rownames(rep4))])
rep4 <- subset(rep4, subset = nFeature_RNA  > 200)

rep4$protocol <- "rep4"
rep4$tech <- "sn-RNA"
rep4 <- NormalizeData(rep4)
rep4 <- FindVariableFeatures(rep4, selection.method = "vst", nfeatures = 2000)
#############
rep5.data <-Read10X(data.dir = "../1_RNA_cellranger_current/sNucRNA-seq_rep5/outs/filtered_feature_bc_matrix")
rep5 <- CreateSeuratObject(counts = rep5.data, project = "rep5", min.cells = 5)
rep5 <- RenameCells(object = rep5, new.names = gsub("-1","-rep5",colnames(rep5)))
rep5[["percent.mt"]] <- PercentageFeatureSet(rep5, pattern = "^ATM")

doublets_rep5 <- read.table("../1_RNA_cellranger_current/filter_doublet/sNucRNA-seq_rep5.doublet.txt")
rep5 <- subset(rep5, cells=WhichCells(rep5, doublets_rep5[,1]), invert=TRUE)
rep5 <- subset(rep5, cells = Cells(rep5), features = rownames(rep5)[grep('AT[1-5]G', rownames(rep5))])
rep5 <- subset(rep5, subset = nFeature_RNA  > 200)

rep5$protocol <- "rep5"
rep5$tech <- "sn-RNA"
rep5 <- NormalizeData(rep5)
rep5 <- FindVariableFeatures(rep5, selection.method = "vst", nfeatures = 2000)

#####################################################################################
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(WT1,WT2,WT3,rep1,rep2,rep3,rep4,rep5))

protocol.anchors <- FindIntegrationAnchors(object.list = list(WT1,WT2,WT3,rep1,rep2,rep3,rep4,rep5),
                                           dims = 1:50,
                                           anchor.features = features,
                                           reduction = "cca")

protocol.combined <- IntegrateData(anchorset = protocol.anchors, dims = 1:50)

DefaultAssay(protocol.combined) <- "integrated"

protocol.combined <- ScaleData(protocol.combined, verbose = FALSE)
protocol.combined <- RunPCA(protocol.combined, npcs = 50, verbose = FALSE)
protocol.combined <- RunUMAP(protocol.combined, reduction = "pca", dims = 1:50)
protocol.combined <- FindNeighbors(protocol.combined, reduction = "pca", dims = 1:50)
protocol.combined <- FindClusters(protocol.combined, resolution = 0.5)

pdf("splitdimplot-umap.pdf", width =12, height = 8)
DimPlot(protocol.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

pdf("splitdimplot-tech-umap.pdf", width =12, height = 8)
DimPlot(protocol.combined, reduction = "umap", group.by = "tech")
dev.off()

pdf("vlnPlot-Genes.pdf")
VlnPlot(protocol.combined, features=c("nFeature_RNA"), pt.size=0)
dev.off()
pdf("vlnPlot-UMI.pdf")
VlnPlot(protocol.combined, features=c("nCount_RNA"), pt.size=0)
dev.off()
pdf("vlnPlot-Mito.pdf")
VlnPlot(protocol.combined, features=c("percent.mt"), pt.size=0)
dev.off()



write.table(table(protocol.combined@meta.data$seurat_clusters, protocol.combined@meta.data$orig.ident), "cluster_counts.txt", sep="\t", quote=F)
write.table(prop.table(table(protocol.combined@meta.data$seurat_clusters, protocol.combined@meta.data$orig.ident),2), "cluster_proportions.by_sample.txt", sep="\t", quote=F)
saveRDS(protocol.combined, "protocol.combined.rds")
write.csv(FetchData(protocol.combined, "ident"), "seurat_bc_clustermap.csv", quote=F)
write.csv(protocol.combined@reductions$umap@cell.embeddings, file = "umap.csv", quote=F)






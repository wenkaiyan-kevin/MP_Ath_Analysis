
library(Signac)
library(Seurat)
library(ggplot2)
set.seed(1234)

counts <- Read10X_h5(filename = "/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep1/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep1/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)
n5 <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)
n5$lib <- 'n5'
fragment.path <- '/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep1/outs/fragments.tsv.gz'
n5 <- SetFragments(
  object = n5,
  file = fragment.path
)
n5 <- NucleosomeSignal(object = n5, region='1-1-30427671')
n5$pct_reads_in_peaks <- n5$peak_region_fragments / n5$total * 100
n5 <- subset(n5, subset = peak_region_fragments > 1000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & nucleosome_signal < 10)
n5 <- RunTFIDF(n5)
n5 <- FindTopFeatures(n5, min.cutoff = 'q0')
n5 <- RunSVD(
  object = n5,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
n5 <- RunUMAP(object = n5, reduction = 'lsi', dims = 1:30)


counts <- Read10X_h5(filename = "/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep2/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep2/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)
g10 <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)
g10$lib <- 'g10'
fragment.path <- '/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep2/outs/fragments.tsv.gz'
g10 <- SetFragments(
  object = g10,
  file = fragment.path
)
g10 <- NucleosomeSignal(object = g10, region='1-1-30427671')
g10$pct_reads_in_peaks <- g10$peak_region_fragments / g10$total * 100
g10 <- subset(g10, subset = peak_region_fragments > 1000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & nucleosome_signal < 10)
g10 <- RenameCells(object = g10, add.cell.id = "B")
g10 <- RunTFIDF(g10)
g10 <- FindTopFeatures(g10, min.cutoff = 'q0')
g10 <- RunSVD(
  object = g10,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
g10 <- RunUMAP(object = g10, reduction = 'lsi', dims = 1:30)

intersecting.regions <- GetIntersectingFeatures(
  object.1 = n5,
  object.2 = g10,
  sep.1 = c(":", "-"),
  sep.2 = c(":", "-")
)



#temporary punt- try Harmony integration
atac_combined <- MergeWithRegions(
  object.1 = n5,
  object.2 = g10,
  sep.1 = c(":", "-"),
  sep.2 = c(":", "-")
)
atac_combined <- FindTopFeatures(atac_combined, min.cutoff = 0)
atac_combined <- RunTFIDF(atac_combined)
atac_combined <- RunSVD(atac_combined, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
atac_combined <- RunUMAP(atac_combined, dims = 1:20, reduction = 'lsi',umap.method = "umap-learn")
p1 <- DimPlot(atac_combined, group.by = 'lib', pt.size = 0.1) + ggplot2::ggtitle("Unintegrated")

ggsave('Unintegrated.pdf',p1)

library("harmony")
atac_combined <- RunHarmony(
  object = atac_combined,
  group.by.vars = 'lib',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)
atac_combined <- RunUMAP(atac_combined, dims = 1:20, reduction = 'harmony', umap.method = "umap-learn")
#atac_combined <- RunUMAP(atac_combined, dims = 1:20, reduction = 'harmony')

p2 <- DimPlot(atac_combined, group.by = 'lib', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")

ggsave('refilter.harmony_integration.pdf',p2)

pdf("refilter.harmony_integration.pdf")
CombinePlots(plots = list(p1, p2), ncol = 2)
dev.off()

#for clustering only with ATAC, per request of reviewer
atac_combined <- FindNeighbors(object = atac_combined, reduction = 'lsi', dims = 1:20)
atac_combined <- FindClusters(object = atac_combined, verbose = FALSE, algorithm = 3, resolution = 0.5)
saveRDS(atac_combined, "refilter.atac_combined.harmonized_no_RNAseq_integration.RDS")
write.table(FetchData(atac_combined, "ident"), "refilter.atac_harmonized_no_RNAseq_integration.clusters.csv", sep=",", quote=F)
pdf("refilter.atac_harmonized_no_RNAseq_integration.clusters.pdf")
DimPlot(
object=atac_combined,
 group.by = 'ident',
  label = TRUE,
  repel = TRUE) + ggplot2::ggtitle('scATAC-seq self-clustering')
dev.off()

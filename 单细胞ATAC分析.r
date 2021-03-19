# R version 3.6.1
# Seurat version 3.1.5
# Signac version 0.2.5

library(Signac)
library(Seurat)
library(ggplot2)
set.seed(1234)

####################
rep1_counts <- Read10X_h5(filename = "/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep1/outs/filtered_peak_bc_matrix.h5")

rep1_metadata <- read.csv(
  file = "/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep1/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

rep1 <- CreateSeuratObject(
  counts = rep1_counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = rep1_metadata
)

rep1$lib <- 'rep1'

rep1_fragment.path <- '/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep1/outs/fragments.tsv.gz'
rep1 <- SetFragments(object = rep1,file = rep1_fragment.path)

rep1 <- NucleosomeSignal(object = rep1, region='1-1-30427671')
rep1$pct_reads_in_peaks <- rep1$peak_region_fragments / rep1$total * 100
rep1 <- subset(rep1, subset = peak_region_fragments > 1000 & 
               peak_region_fragments < 20000 & 
               pct_reads_in_peaks > 15 & 
               nucleosome_signal < 10)
rep1 <- RunTFIDF(rep1)
rep1 <- FindTopFeatures(rep1, min.cutoff = 'q0')
rep1 <- RunSVD(
  object = rep1,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
rep1 <- RunUMAP(object = rep1, reduction = 'lsi', dims = 1:20)

###########################
rep2_counts <- Read10X_h5(filename = "/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep2/outs/filtered_peak_bc_matrix.h5")

rep2_metadata <- read.csv(
  file = "/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep2/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

rep2 <- CreateSeuratObject(
  counts = rep2_counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = rep2_metadata
)
rep2 <- RenameCells(object = rep2, new.names = gsub("-1","-2",colnames(rep2)))
rep2$lib <- 'rep2'

rep2_fragment.path <- '/gss1/home/yanwk/ywk_work/5_other_Analysis/3_MP_Ath_analysis/3_ATAC_cellranger/sNucATAC-seq_rep2/outs/fragments.tsv.gz'
rep2 <- SetFragments(object = rep2,file = rep2_fragment.path)

rep2 <- NucleosomeSignal(object = rep2, region='1-1-30427671')
rep2$pct_reads_in_peaks <- rep2$peak_region_fragments / rep2$total * 100
rep2 <- subset(rep2, subset = peak_region_fragments > 1000 & 
                              peak_region_fragments < 20000 & 
                              pct_reads_in_peaks > 15 & 
                              nucleosome_signal < 10)

rep2 <- RunTFIDF(rep2)
rep2 <- FindTopFeatures(rep2, min.cutoff = 'q0')
rep2 <- RunSVD(object = rep2,
               assay = 'peaks',
               reduction.key = 'LSI_',
               reduction.name = 'lsi')
rep2 <- RunUMAP(object = rep2, reduction = 'lsi', dims = 1:20)

############################
intersecting.regions <- GetIntersectingFeatures(object.1 = rep1,
                                                object.2 = rep2,
                                                sep.1 = c(":", "-"),
                                                sep.2 = c(":", "-"))

#temporary punt- try Harmony integration
atac_combined <- MergeWithRegions(object.1 = rep1,
                                  object.2 = rep2,
                                  sep.1 = c(":", "-"),
                                  sep.2 = c(":", "-"))

atac_combined <- FindTopFeatures(atac_combined, min.cutoff = 0)
atac_combined <- RunTFIDF(atac_combined)
atac_combined <- RunSVD(atac_combined, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
atac_combined <- RunUMAP(atac_combined, dims = 1:30, reduction = 'lsi')
#atac_combined <- RunUMAP(atac_combined, dims = 1:30, reduction = 'lsi',umap.method = "umap-learn")

p1 <- DimPlot(atac_combined, group.by = 'lib', pt.size = 0.1) + ggplot2::ggtitle("Unintegrated")
ggsave('Unintegrated.pdf',p1)

#######################################
library(harmony)

atac_combined <- RunHarmony(object = atac_combined,
                            group.by.vars = 'lib',
                            reduction = 'lsi',
                            assay.use = 'peaks',
                            project.dim = FALSE)
#atac_combined <- RunUMAP(atac_combined, dims = 1:30, reduction = 'harmony', umap.method = "umap-learn")
atac_combined <- RunUMAP(atac_combined, dims = 1:30, reduction = 'harmony')

p2 <- DimPlot(atac_combined, group.by = 'lib', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
ggsave('refilter.harmony_integration.pdf',p2)

#pdf("refilter.harmony_integration.pdf")
#CombinePlots(plots = list(p1, p2), ncol = 2)
#dev.off()

#for clustering only with ATAC, per request of reviewer
atac_combined <- FindNeighbors(object = atac_combined, reduction = 'harmony', dims = 1:30)
atac_combined <- FindClusters(object = atac_combined,  verbose = FALSE, algorithm = 3)
#atac_combined <- FindClusters(object = atac_combined, verbose = FALSE, algorithm = 3, resolution = 0.5)

p3 <- DimPlot( object=atac_combined, label = TRUE, repel = TRUE) + ggplot2::ggtitle('scATAC-seq self-clustering')
ggsave("refilter.atac_harmonized_no_RNAseq_integration.clusters.pdf", p3)

saveRDS(atac_combined, "refilter.atac_combined.harmonized_no_RNAseq_integration.RDS")
write.table(FetchData(atac_combined, "ident"), "refilter.atac_harmonized_no_RNAseq_integration.clusters.csv", sep=",", quote=F)
umap <- cbind("Barcode" = rownames(Embeddings(object = atac_combined, reduction = "umap")), 
              Embeddings(object = atac_combined, reduction = "umap"))
write.table(umap, file="refilter.atac_harmonized_no_RNAseq_integration.umap.csv", 
            sep = ",", quote = F, row.names = F, col.names = T)







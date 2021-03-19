
library(Signac)
library(Seurat)
set.seed(1234)
counts <- Read10X_h5(filename = "N5_cellranger-1_2_0/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "N5_cellranger-1_2_0/singlecell.csv",
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

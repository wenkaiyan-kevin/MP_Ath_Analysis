#!/usr/bin/env python3
import sys
import numpy as np
import doubletdetection
import tarfile
import matplotlib.pyplot as plt

matrix_path = sys.argv[1]
raw_counts = doubletdetection.load_mtx(matrix_path)
zero_genes = np.sum(raw_counts, axis=0) == 0
raw_counts = raw_counts[:, ~zero_genes]
clf = doubletdetection.BoostClassifier(n_iters=50)
doublets = clf.fit(raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

db = doublets.nonzero()
f = open(sys.argv[2], "w")
for d in db[0] :
    f.write(str(d+1) + "\n")
f.close()



import numpy as np
import doubletdetection
import scanpy as sc
import matplotlib.pyplot as plt

matrix_path = sys.argv[1]
adata = sc.read_10x_h5(
    "pbmc_10k_v3_filtered_feature_bc_matrix.h5", 
    backup_url="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
)
adata.var_names_make_unique()

# remove "empty" genes
sc.pp.filter_genes(adata, min_cells=1)

clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True)
doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.8)
#doublet_score = clf.doublet_score()

db = doublets.nonzero()
f = open(sys.argv[2], "w")
for d in db[0] :
    f.write(str(d+1) + "\n")
f.close()






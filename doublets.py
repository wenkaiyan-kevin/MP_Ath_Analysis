import numpy as np
import doubletdetection
import scanpy as sc
import matplotlib.pyplot as plt

matrix_path = sys.argv[1]
adata = sc.read_10x_h5(matrix_path)
adata.var_names_make_unique()

# remove "empty" genes
sc.pp.filter_genes(adata, min_cells=1)

clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True)
doublets = clf.fit(adata.X).predict(p_thresh=1e-7, voter_thresh=0.8)
#doublet_score = clf.doublet_score()

db = doublets.nonzero()
f = open(sys.argv[2], "w")
for d in db[0] :
    f.write(str(d+1) + "\n")
f.close()






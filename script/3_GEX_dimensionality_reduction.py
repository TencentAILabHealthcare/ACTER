import os

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler, StandardScaler


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80)

samples = [
    "D10-Ln-Na",
    "D10-Sp-Na",
    "D10-Tm-Na",
    "D20-Ln-Na",
    "D20-Ln-Neo",
    "D20-Ln-NeoPD",
    "D20-Ln-PD",
    "D20-Sp-Na",
    "D20-Tm-Na",
    "D20-Tm-Neo",
    "D20-Tm-NeoPD",
    "D20-Tm-PD",
]


adata = sc.read_h5ad('../Result/GSE178881/inter_file/gex_scanpy_qc_filtered.h5ad')
print(adata)

# normalize to depth 10 000
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
# logaritmize
sc.pp.log1p(adata)

# store normalized counts in the raw slot, 
# we will subset adata.X for variable genes, but want to keep all genes matrix as well.
adata.raw = adata

# compute variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Highly variable genes: %d"%sum(adata.var.highly_variable))

#plot variable genes
sc.pl.highly_variable_genes(adata)
plt.savefig("../Result/GSE178881/figure/gex_highly_variable_genes.pdf")

# subset for variable genes in the dataset
adata = adata[:, adata.var['highly_variable']]

# regress out unwanted variables
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# scale data, clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50)
plt.savefig("../Result/GSE178881/figure/gex_PCA_ranking.pdf")

# UMAP
sc.pp.neighbors(adata, n_pcs = 50, n_neighbors = 20)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['sample', 'day', 'tissue', 'treat'], ncols=1)
plt.savefig("../Result/GSE178881/figure/gex_UMAP.pdf")


# Save data
print(adata.X.shape)
print(adata.raw.X.shape)
adata.write_h5ad('../Result/GSE178881/inter_file/gex_scanpy_dr.h5ad')

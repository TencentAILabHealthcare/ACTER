# -*- coding: utf-8 -*-
import os
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler, StandardScaler

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80)


adata = sc.read_h5ad("../Result/GSE178881/inter_file/gexvdj_scanpy_mnn_corrected.h5ad")
print(adata)

sc.tl.leiden(adata, key_added="leiden_1.0")  # default resolution in 1.0
sc.tl.leiden(adata, resolution=0.6, key_added="leiden_0.6")
sc.tl.leiden(adata, resolution=0.4, key_added="leiden_0.4")
sc.tl.leiden(adata, resolution=1.4, key_added="leiden_1.4")

sc.pl.umap(adata, color=["leiden_0.4", "leiden_0.6", "leiden_1.0", "leiden_1.4"])
plt.savefig("../Result/GSE178881/figure/gexvdj_clustering.pdf")
plt.close()

tmp04 = pd.crosstab(adata.obs["leiden_0.4"], adata.obs["treat"], normalize="index")
fig=plt.figure(figsize=(18, 10))
axis=fig.add_axes([0.1,0.1,0.8,0.8])
tmp04.plot.bar(ax=axis, stacked=True).legend(bbox_to_anchor=(1.1, 1.0), loc="upper right")
plt.savefig("../Result/GSE178881/figure/gexvdj_clustering_stat04.pdf")
plt.close()

tmp06 = pd.crosstab(adata.obs["leiden_0.6"], adata.obs["treat"], normalize="index")
tmp06.plot.bar(stacked=True, figsize=(18, 10)).legend(bbox_to_anchor=(1.1, 1.0), loc="upper right")
plt.savefig("../Result/GSE178881/figure/gexvdj_clustering_stat06.pdf")
plt.close()

tmp10 = pd.crosstab(adata.obs["leiden_1.0"], adata.obs["treat"], normalize="index")
tmp10.plot.bar(stacked=True, figsize=(18, 10)).legend(bbox_to_anchor=(1.1, 1.0), loc="upper right")
plt.savefig("../Result/GSE178881/figure/gexvdj_clustering_stat10.pdf")
plt.close()

tmp14 = pd.crosstab(adata.obs["leiden_1.4"], adata.obs["treat"], normalize="index")
tmp14.plot.bar(stacked=True, figsize=(18, 10)).legend(bbox_to_anchor=(1.1, 1.0), loc="upper right")
plt.savefig("../Result/GSE178881/figure/gexvdj_clustering_stat14.pdf")
plt.close()

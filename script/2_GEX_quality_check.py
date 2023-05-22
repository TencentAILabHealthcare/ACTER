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


adata = sc.read_h5ad('../Result/GSE178881/merged_gex.h5ad')
print(adata.obs)

adata.obs["sample"] = adata.obs["batch"].cat.rename_categories({str(i[0]): i[1] for i in zip(range(len(samples)), samples)})
print(adata.obs)

adata.obs["day"] = [i.split("-")[0] for i in adata.obs["sample"].tolist()]
adata.obs["day"] = adata.obs["day"].astype("category")

adata.obs["tissue"] = [i.split("-")[1] for i in adata.obs["sample"].tolist()]
adata.obs["tissue"] = adata.obs["tissue"].astype("category")

adata.obs["treat"] = [i.split("-")[2] for i in adata.obs["sample"].tolist()]
adata.obs["treat"] = adata.obs["treat"].astype("category")

print(adata.obs)
print(adata.obs.dtypes)
print(adata.obs['sample'].value_counts())
print(adata.var)

# mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('mt-') 
# ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(("Rps","Rpl"))

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
print(adata)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
    jitter=0.4,
    groupby="sample",
    rotation=45,
)
plt.savefig("../Result/GSE178881/figure/n_genes_by_counts.pdf")

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(adata.n_obs, adata.n_vars)

sc.pl.highest_expr_genes(adata, n_top=20)
plt.savefig("../Result/GSE178881/figure/highest_expr_genes.pdf")

# filter for percent mito
adata = adata[adata.obs['pct_counts_mt'] < 20, :]
# filter for percent ribo > 0.05
adata = adata[adata.obs['pct_counts_ribo'] > 5, :]
print("Remaining cells %d"%adata.n_obs)


# Predict doublets
adata.raw = adata
import scrublet as scr
scrub = scr.Scrublet(adata.raw.X)
adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
scrub.plot_histogram()
plt.savefig("../Result/GSE178881/figure/predict_doublets.pdf")
sum(adata.obs['predicted_doublets'])

# add in column with singlet/doublet instead of True/False
adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)

# also revert back to the raw counts as the main matrix in adata
adata = adata.raw.to_adata() 
adata = adata[adata.obs['doublet_info'] == 'False',:]
print(adata.shape)

# Save QC-filtered data
save_file = '../Result/GSE178881/inter_file/scanpy_qc_filtered.h5ad'
adata.write_h5ad(save_file)

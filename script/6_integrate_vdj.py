import sys
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

adata = sc.read_h5ad("../Result/GSE178881/inter_file/gex_scanpy_dr.h5ad")
print(adata.var)

print(adata.X.shape)
adata3 = adata.raw.to_adata()
print(adata3.X.shape)

var_genes_all = adata.var.highly_variable
print("Highly variable genes: %d"%sum(var_genes_all))

sc.pp.highly_variable_genes(adata3, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'sample')
print("Highly variable genes intersection: %d"%sum(adata3.var.highly_variable_intersection))
print("Number of batches where gene is variable:")
print(adata3.var.highly_variable_nbatches.value_counts())
var_genes_batch = adata3.var.highly_variable_nbatches > 0
print("Any batch var genes: %d"%sum(var_genes_batch))
print("All data var genes: %d"%sum(var_genes_all))
print("Overlap: %d"%sum(var_genes_batch & var_genes_all))
print("Variable genes in all batches: %d"%sum(adata3.var.highly_variable_nbatches == 12))
print("Overlap batch instersection and all: %d"%sum(var_genes_all & adata3.var.highly_variable_intersection))
var_select = adata3.var.highly_variable_nbatches > 6
var_genes = var_select.index[var_select]
print("var_genes:", len(var_genes))


# Read TCR
clones_file = "../Result/GSE178881/merged_clones.tsv"
tmpdata = open(clones_file, "r")
clones_file_header = tmpdata.readline()
tmpdata.close()
clones_df = pd.read_csv(clones_file, sep="\t")

# Check format
for vj in "vj":
    for ab in "ab":
        tag = f"{vj}{ab}_gene"
        alt_tag = f"{vj}{ab}"
        if tag not in clones_df.columns and alt_tag in clones_df.columns:
            clones_df.rename(columns={alt_tag: tag}, inplace=True)

CLONES_FILE_REQUIRED_COLUMNS = "clone_id va_gene ja_gene cdr3a cdr3a_nucseq vb_gene jb_gene cdr3b cdr3b_nucseq".split()
for colname in CLONES_FILE_REQUIRED_COLUMNS:
    if colname not in clones_df.columns:
        print("ERROR clones_file is missing required column:", colname)
        print("Required columns:", CLONES_FILE_REQUIRED_COLUMNS)
        print("clones_file:", clones_file)
        sys.exit()



id2tcr = {}
tcr2id = {}

include_tcr_nucseq = True

for l in clones_df.itertuples(index=False):
    if include_tcr_nucseq:
        atcr = (l.va_gene, l.ja_gene, l.cdr3a, l.cdr3a_nucseq)
        btcr = (l.vb_gene, l.jb_gene, l.cdr3b, l.cdr3b_nucseq)
    else:
        atcr = (l.va_gene, l.ja_gene, l.cdr3a)  # , l.cdr3a_nucseq )
        btcr = (l.vb_gene, l.jb_gene, l.cdr3b)  # , l.cdr3b_nucseq )
    tcr = (atcr, btcr)
    # tcr = ( atcr, btcr, l.clone_id ) # hack to guarantee uniqueness!!!
    id2tcr[l.clone_id] = tcr
    tcr2id[tcr] = l.clone_id


# read the barcode/clonotype mapping info
barcode2tcr = {}
for line in open('../Result/GSE178881/merged_clones.tsv.barcode_mapping.tsv', "r"):
    l = line[:-1].split("\t")
    if l[0] == "clone_id":
        continue  # header line
    if not l[1]:
        continue
    barcodes = l[1].split(",")
    clone_id = l[0]
    if clone_id not in id2tcr: continue # maybe short cdr3?
    tcr = id2tcr[clone_id]
    for bc in barcodes:
        barcode2tcr[bc] = tcr
        assert bc in barcodes  # the barcodes list before preprocessing...


# Reducing
def store_tcrs_in_adata(adata, tcrs):
    """ returns NOTHING, modifies adata
    tcrs is a list of (atcr,btcr) tuples, where
     atcr = (va,ja,cdr3a,cdr3a_nucseq) ...
    """
    assert len(tcrs) == adata.shape[0]

    tcr_keys = "va ja cdr3a cdr3a_nucseq vb jb cdr3b cdr3b_nucseq".split()
    tcr_indices = ((x, y) for x in range(2) for y in range(4))

    for tag, (i, j) in zip(tcr_keys, tcr_indices):
        adata.obs[tag] = [x[i][j] for x in tcrs]

    # ensure lower case
    adata.obs["cdr3a_nucseq"] = adata.obs.cdr3a_nucseq.str.lower()
    adata.obs["cdr3b_nucseq"] = adata.obs.cdr3b_nucseq.str.lower()

    return


print(adata3.obs.index)

mask = [x in barcode2tcr for x in adata3.obs.index]
print(f"Reducing to the {np.sum(mask)} barcodes (out of {adata3.shape[0]}) with paired TCR sequence data")
# adata = adata[mask,:]
adata_w_TCR = adata3[mask, :].copy()
tcrs = [barcode2tcr[x] for x in adata_w_TCR.obs.index]
store_tcrs_in_adata(adata_w_TCR, tcrs)
print(adata_w_TCR)
print(adata_w_TCR.obs)

adata = adata_w_TCR

print(adata.obs['cdr3a'].str.len().describe())
print(adata.obs['cdr3b'].str.len().describe())

af_table = pd.read_csv('../Result/GSE178881/Atchley_factors.csv', index_col='aa')
print(af_table)


def cal_af_full(peptides):
    af = []
    for peptide in peptides:
        tmp = []
        for aa in peptide:
            tmp.extend(list(af_table.loc[aa]))
        af.append(tmp)
    af = pd.DataFrame(af)
    af.fillna(0, inplace=True)
    return af


cdr3a_af_full = cal_af_full(adata.obs["cdr3a"])
cdr3b_af_full = cal_af_full(adata.obs["cdr3b"])

print(adata)

X_new = np.concatenate((adata.X.toarray(), np.array(cdr3a_af_full), np.array(cdr3b_af_full)), axis=1)
print(X_new)

l = cdr3a_af_full.shape[1] + cdr3b_af_full.shape[1]
d = {"gene_ids": ["AF" + str(i) for i in range(l)], "feature_types": ["Atchley Factors"] * l, "highly_variable": [True] * l}
var_af = pd.DataFrame(d)
var_af.index = var_af["gene_ids"]
var_af['feature_types'] = var_af['feature_types'].astype('category')
print(var_af.dtypes)

var_old = adata.var[["gene_ids", "feature_types", "highly_variable"]].copy()
var_old['highly_variable'] = False
var_old.loc[var_genes, 'highly_variable'] = True
var_new = pd.concat([var_old, var_af])
print(var_new)

adata_gex_vdj = ad.AnnData(X_new, obs=adata.obs, var=var_new)
print(adata_gex_vdj)
print(adata_gex_vdj.var)


# Integration
# split per batch into new objects.
batches = adata_gex_vdj.obs['sample'].cat.categories.tolist()
alldata = {}
for batch in batches:
    alldata[batch] = adata_gex_vdj[adata_gex_vdj.obs['sample'] == batch,]
print(alldata)


var_select = adata_gex_vdj.var.highly_variable == True
var_select_gex_vdj = var_select.index[var_select]
print(len(var_select_gex_vdj))

cdata = sc.external.pp.mnn_correct(*alldata.values(), svd_dim=50, batch_key="sample", save_raw=True, var_subset=var_select_gex_vdj)

corr_data = cdata[0][:,var_select_gex_vdj]
print(corr_data.X.shape)

corr_data.obs["sample"] = corr_data.obs["sample"].cat.rename_categories({str(i[0]): i[1] for i in zip(range(len(samples)), samples)})
corr_data.obs["batch"] = corr_data.obs["batch"].astype("category")
corr_data.obs["day"] = corr_data.obs["day"].astype("category")
corr_data.obs["tissue"] = corr_data.obs["tissue"].astype("category")
corr_data.obs["treat"] = corr_data.obs["treat"].astype("category")

print(corr_data.obs)

# the variable genes defined are used by default by the pca function, 
# now we want to run on all the genes in the dataset
sc.tl.pca(corr_data, svd_solver = 'arpack', use_highly_variable = False)
sc.pl.pca(corr_data, components = ['1,2','3,4','5,6','7,8'], ncols=1, color='sample')
plt.savefig("../Result/GSE178881/figure/gexvdj_PCAs.pdf")
plt.close()

# tSNE
sc.tl.tsne(corr_data, n_pcs = 50)
# UMAP, first with neighbor calculation 
sc.pp.neighbors(corr_data, n_pcs = 50, n_neighbors = 20)
sc.tl.umap(corr_data)

sc.tl.pca(adata_gex_vdj, svd_solver = 'arpack', use_highly_variable = False)
# tSNE
sc.tl.tsne(adata_gex_vdj, n_pcs = 50)
# UMAP, first with neighbor calculation 
sc.pp.neighbors(adata_gex_vdj, n_pcs = 50, n_neighbors = 20)
sc.tl.umap(adata_gex_vdj)

fig, axs = plt.subplots(2, 2, figsize=(12,10),constrained_layout=True)
sc.pl.tsne(corr_data, color="sample", title="MNN Corrected tsne", ax=axs[0,0], show=False)
sc.pl.tsne(adata_gex_vdj, color="sample", title="Uncorrected tsne", ax=axs[0,1], show=False)
sc.pl.umap(corr_data, color="sample", title="MNN Corrected umap", ax=axs[1,0], show=False)
sc.pl.umap(adata_gex_vdj, color="sample", title="Uncorrected umap", ax=axs[1,1], show=False)
plt.savefig("../Result/GSE178881/figure/gexvdj_MNN.pdf")
plt.close()

corr_data.write_h5ad('../Result/GSE178881/inter_file/gexvdj_scanpy_mnn_corrected.h5ad')

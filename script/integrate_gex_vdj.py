# -*- coding: utf-8 -*-
import argparse
import os
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr

# from sklearn.decomposition import PCA
# from sklearn.preprocessing import MinMaxScaler, StandardScaler

parser = argparse.ArgumentParser(description="Integrate GEX and VDJ information into AnnData object.")
parser.add_argument(
    "--samples",
    help="tsvfile with 4 columns: sample name, clones file, gex data, gex data type",
    default=None,
    required=True,
)
parser.add_argument("--merged_clones_file", default=None, required=True)
parser.add_argument("--merged_gex_data", default=None, required=True)
parser.add_argument("--organism", choices=["human", "mouse"], required=True)
args = parser.parse_args()

samples_df = pd.read_csv(args.samples, sep="\t")
samples = list(samples_df.sample_name)
clones_df = pd.read_csv(args.merged_clones_file, sep="\t")
bcmap_df = pd.read_csv(args.merged_clones_file + ".barcode_mapping.tsv", sep="\t")
adata = sc.read_h5ad(args.merged_gex_data)

id2tcr = {}
tcr2id = {}
for line in clones_df.itertuples(index=False):
    tcra = (line.va_gene, line.ja_gene, line.cdr3a, line.cdr3a_nucseq)
    tcrb = (line.vb_gene, line.jb_gene, line.cdr3b, line.cdr3b_nucseq)
    tcr = (tcra, tcrb)
    id2tcr[line.clone_id] = tcr
    tcr2id[tcr] = line.clone_id

barcode2tcr = {}
for line in bcmap_df.itertuples(index=False):
    barcodes = line.barcodes.split(",")
    clone_id = line.clone_id
    # if clone_id not in id2tcr: continue # maybe short cdr3?
    tcr = id2tcr[clone_id]
    for bc in barcodes:
        barcode2tcr[bc] = tcr

mask = [x in barcode2tcr for x in adata.obs.index]
print(f"Reducing to the {np.sum(mask)} barcodes (out of {adata.shape[0]}) with paired TCR sequence data")
adata = adata[mask, :].copy()

# add tcr/vdj info to adata.obs
tcrs = [barcode2tcr[x] for x in adata.obs.index]
tcr_keys = "va ja cdr3a cdr3a_nucseq vb jb cdr3b cdr3b_nucseq".split()
tcr_indices = ((x, y) for x in range(2) for y in range(4))
for tag, (i, j) in zip(tcr_keys, tcr_indices):
    adata.obs[tag] = [x[i][j] for x in tcrs]
adata.obs["cdr3a_nucseq"] = adata.obs.cdr3a_nucseq.str.lower()
adata.obs["cdr3b_nucseq"] = adata.obs.cdr3b_nucseq.str.lower()

# sample such as D10-Ln-Na, D10-Sp-Na
adata.obs["sample"] = adata.obs["batch"].cat.rename_categories(
    {str(i[0]): i[1] for i in zip(range(len(samples)), samples)}
)

adata.obs["day"] = [i.split("-")[0] for i in adata.obs["sample"].tolist()]
adata.obs["day"] = adata.obs["day"].astype("category")

adata.obs["tissue"] = [i.split("-")[1] for i in adata.obs["sample"].tolist()]
adata.obs["tissue"] = adata.obs["tissue"].astype("category")

adata.obs["treatment"] = [i.split("-")[2] for i in adata.obs["sample"].tolist()]
adata.obs["treatment"] = adata.obs["treatment"].astype("category")

# filtered out cells that have less than 200 genes expressed
sc.pp.filter_cells(adata, min_genes=200)
# filtered out 14115 genes that are detected in less than 3 cells
sc.pp.filter_genes(adata, min_cells=3)

if args.organism == "mouse":
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
elif args.organism == "human":
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
else:
    sys.exit()

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)
# filter for percent mito
adata = adata[adata.obs["pct_counts_mt"] < 20, :]
# filter for percent ribo > 0.05
adata = adata[adata.obs["pct_counts_ribo"] > 5, :]

print(f"Remaining {adata.n_obs} cells %d")

# Predict doublets
adata.raw = adata
import scrublet as scr
scrub = scr.Scrublet(adata.raw.X)
adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()

# add column with singlet/doublet instead of True/False
adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)

# also revert back to the raw counts as the main matrix in adata
adata = adata.raw.to_adata() 
adata = adata[adata.obs['doublet_info'] == 'False',:]
print(adata.shape)

# Save QC-filtered data
save_file = '../Result/GSE178881/inter_file/scanpy_qc_filtered.h5ad'
adata.write_h5ad(save_file)

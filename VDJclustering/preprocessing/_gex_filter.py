# -*- coding: utf-8 -*-
import scanpy as sc


def gex_filter(adata, min_genes, min_cells, pct_counts_mt, pct_counts_ribo):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    # filter for percent mito
    adata = adata[adata.obs['pct_counts_mt'] < 20, :]
    # filter for percent ribo
    adata = adata[adata.obs['pct_counts_ribo'] > 5, :]
    print(f"Remaining cells {adata.n_obs}, genes {adata.n_vars}")
    return adata

# -*- coding: utf-8 -*-
import scanpy as sc


def pca_neighbors(adata: sc.AnnData, svd_solver="arpack", use_highly_variable=True, n_comps=50, n_neighbors=20):
    sc.tl.pca(adata, svd_solver=svd_solver, use_highly_variable=use_highly_variable, n_comps=n_comps)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps)

    return adata

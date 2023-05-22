# -*- coding: utf-8 -*-
import scanpy as sc


def gex_norm_logtrans(adata: sc.AnnData, counts_per_cell_after: int = 1e4, save_raw: bool = True):
    if save_raw:
        adata.raw = adata
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=counts_per_cell_after)
    sc.pp.log1p(adata)

    return adata

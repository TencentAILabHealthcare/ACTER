# -*- coding: utf-8 -*-
import scanpy as sc


def zscore_transform(adata: sc.AnnData, save_raw: bool = True):
    if save_raw:
        adata.raw = adata
    sc.pp.scale(adata, max_value=10)

    return adata

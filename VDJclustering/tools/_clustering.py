# -*- coding: utf-8 -*-
import scanpy as sc


def leiden(adata: sc.AnnData, resolution: str, key_added: str):
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)

# -*- coding: utf-8 -*-
import scanpy as sc


def umap(adata: sc.AnnData):
    sc.pl.umap(adata)

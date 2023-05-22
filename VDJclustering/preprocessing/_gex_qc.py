# -*- coding: utf-8 -*-
import sys

import scanpy as sc


def _mt_ribo_names(organism):
    """Based on organism, return a dict to specify mitochondrial and ribosomal gene names."""
    mt_ribo_gene_names = {}
    if organism == "human":
        mt_ribo_gene_names["mt"] = "MT-"
        mt_ribo_gene_names["ribo"] = ("RPS", "RPL")
    elif organism == "mouse":
        mt_ribo_gene_names["mt"] = "mt-"
        mt_ribo_gene_names["ribo"] = ("Rps", "Rpl")
    else:
        print(f"Cannot recognize {organism}")
        print(f"organism must be one of ['human','mouse']")
        sys.exit()
    return mt_ribo_gene_names


def gex_calculate_qc_metrics(adata, organism, inplace=True):
    mt_ribo_gene_names = _mt_ribo_names(organism)
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith(mt_ribo_gene_names["mt"])
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(mt_ribo_gene_names["ribo"])
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=inplace)

    return adata

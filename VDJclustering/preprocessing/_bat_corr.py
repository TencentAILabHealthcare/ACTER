# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd


def batch_correct_mnn(adata: sc.AnnData, sample_config: str, use_highly_variable: bool = True, save_raw: bool = True):
    if use_highly_variable:
        adata = adata[:, adata.var["highly_variable"] == True]

    adata.raw = adata

    # split per batch into new objects.
    batches = adata.obs["batch"].cat.categories.tolist()
    alldata = {}
    for batch in batches:
        alldata[batch] = adata[
            adata.obs["batch"] == batch,
        ]

    var_genes = None
    # if use_highly_variable:
    #     var_select = adata.var.highly_variable == True
    #     var_genes = var_select.index[var_select]

    cdata = sc.external.pp.mnn_correct(
        *alldata.values(), svd_dim=50, batch_key="batch", var_subset=var_genes, save_raw=save_raw
    )
    corr_data = cdata[0]

    corr_data.obs["sample"] = corr_data.obs["sample"].astype("category")
    samples_df = pd.read_csv(sample_config, sep="\t")
    for cat in samples_df.columns[4:]:
        corr_data.obs[cat] = corr_data.obs[cat].astype("category")

    return corr_data


def batch_correct_mnn2(adata: sc.AnnData, sample_config: str, use_highly_variable: bool = True, save_raw: bool = True):
    """For test only, not use in practice"""
    adata.raw = adata
    # if use_highly_variable:
    #     adata = adata[:, adata.var["highly_variable"] == True]

    # split per batch into new objects.
    batches = adata.obs["batch"].cat.categories.tolist()
    alldata = {}
    for batch in batches:
        alldata[batch] = adata[
            adata.obs["batch"] == batch,
        ]

    var_genes = None
    if use_highly_variable:
        var_select = adata.var["highly_variable"] == True
        var_genes = var_select.index[var_select]

    cdata = sc.external.pp.mnn_correct(
        *alldata.values(), svd_dim=50, batch_key="batch", var_subset=var_genes, save_raw=save_raw
    )
    corr_data = cdata[0]

    corr_data.obs["sample"] = corr_data.obs["sample"].astype("category")
    samples_df = pd.read_csv(sample_config, sep="\t")
    for cat in samples_df.columns[4:]:
        corr_data.obs[cat] = corr_data.obs[cat].astype("category")

    return corr_data

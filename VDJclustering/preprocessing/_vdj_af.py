# -*- coding: utf-8 -*-

import os

import numpy as np
import pandas as pd
import scanpy as sc


def _cal_af(peptides):
    """Calculate Atchley factors"""
    wd = os.path.dirname(__file__)
    af_df = pd.read_csv(os.path.join(wd, "data/Atchley_factors.csv"), index_col="aa")
    crd3_af = []
    for peptide in peptides:
        line = []
        for aa in peptide:
            line.extend(list(af_df.loc[aa]))
        crd3_af.append(line)
    crd3_af = pd.DataFrame(crd3_af)

    crd3_af.fillna(0, inplace=True)

    # fill NaN with the forward value in row
    # crd3_af.fillna(method="ffill",axis=1, inplace=True)
    
    # fill NaN with column mean
    # for column in list(crd3_af.columns[crd3_af.isnull().sum() > 0]):
    #     mean_val = crd3_af[column].mean()
    #     crd3_af[column].fillna(mean_val, inplace=True)

    # Min(af_df) == -4.7596375
    # Max(af_df) == 3.0973596
    # crd3_af = crd3_af + 5

    return crd3_af


def generate_gex_vdj(adata, var_genes):
    cdr3a_af = _cal_af(adata.obs["cdr3a"])
    cdr3b_af = _cal_af(adata.obs["cdr3b"])
    cdr3ab_len = cdr3a_af.shape[1] + cdr3b_af.shape[1]
    var_cdr3ab = pd.DataFrame(
        {
            "gene_ids": ["AF" + str(i) for i in range(cdr3ab_len)],
            "feature_types": ["Atchley Factor"] * cdr3ab_len,
            "highly_variable": [True] * cdr3ab_len,
        }
    )
    var_cdr3ab.index = var_cdr3ab["gene_ids"]
    var_cdr3ab['feature_types'] = var_cdr3ab['feature_types'].astype('category')

    var_gex = adata.var[["gene_ids", "feature_types", "highly_variable"]].copy()
    var_gex['highly_variable'] = False
    var_gex.loc[var_genes, 'highly_variable'] = True
    var_gex_cdr3ab = pd.concat([var_gex, var_cdr3ab])
    X = np.concatenate((adata.X.toarray(), np.array(cdr3a_af), np.array(cdr3b_af)), axis=1)
    adata_gex_cdr3ab = sc.AnnData(X, obs=adata.obs, var=var_gex_cdr3ab)

    return adata_gex_cdr3ab

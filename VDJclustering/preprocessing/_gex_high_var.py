# -*- coding: utf-8 -*-
# from dataclasses import dataclass
import scanpy as sc


def highly_variable_genes_all(adata: sc.AnnData):
    """Detect variable genes can be detected across the full dataset, but then we run the risk of getting many batch-specific genes that will drive a lot of the variation."""
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    var_genes_all = adata.var.highly_variable
    print(f"Highly variable genes (across the full dataset): {sum(var_genes_all)}")
    return adata


def highly_variable_genes_sample(adata: sc.AnnData):
    """Detect variable genes in each sample indepently"""
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="sample")
    print(
        f"Highly variable genes intersection (i.e., HVGs that are common in all samples): {sum(adata.var.highly_variable_intersection)}"
    )
    print("Number of samples where gene is variable:")
    print(adata.var.highly_variable_nbatches.value_counts())
    var_genes_batch = adata.var.highly_variable_nbatches > 0
    print(f"Any batch var genes: {sum(var_genes_batch)}")
    print(f"Variable genes in all batches: {sum(adata.var.highly_variable_nbatches == 12)}")
    return adata

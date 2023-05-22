# -*- coding: utf-8 -*-
import pandas as pd


def show_sample_category(anndata, sample_config):
    samples_df = pd.read_csv(sample_config, sep="\t")
    for cat in samples_df.columns[4:]:
        print(f"Catetory: {cat}")
        print(anndata.obs[cat].value_counts())
        print("=" * 40)

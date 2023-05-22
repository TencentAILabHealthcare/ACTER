# -*- coding: utf-8 -*-
import argparse
import os
import sys

import pandas as pd
import scanpy as sc


def read_gex_data(gex_data, gex_data_type, gex_only=True):
    print(f"reading {gex_data_type} type data {gex_data}")
    if gex_data_type == "10x_matx":
        adata = sc.read_10x_mtx(gex_data, gex_only=gex_only)
    elif gex_data_type == "10x_h5":
        adata = sc.read_10x_h5(gex_data, gex_only=gex_only)
    elif gex_data_type == "h5ad":
        adata = sc.read_h5ad(gex_data)
    else:
        print("Unrecognized gex_data_type: {gex_data_type}. ['h5ad', '10x_mtx', '10x_h5']")
        sys.exit()
    return adata


parser = argparse.ArgumentParser(
    description="Merge multiple datasets and generate a single clones_file and gex_data file."
)
parser.add_argument(
    "--samples",
    help="tsvfile with 4 columns: sample name, clones file, gex data, gex data type",
    default=None,
    required=True,
)
parser.add_argument("--output_merged_clones_file", default=None, required=True)
parser.add_argument("--output_merged_gex_data", default = None, required=True)
args = parser.parse_args()

all_data = []

df = pd.read_csv(args.samples, sep="\t")
for line in df.itertuples():
    if not os.path.exists(line.clones_file):
        print(f"Not found: {line.clones_file}")
        sys.exit()
    bcmap_file = line.clones_file + ".barcode_mapping.tsv"
    if not os.path.exists(bcmap_file):
        print(f"Not found: {bcmap_file}")
        sys.exit()

    clones_df = pd.read_csv(line.clones_file, sep="\t")
    bcmap_df = pd.read_csv(bcmap_file, sep="\t")
    adata = read_gex_data(line.gex_data, line.gex_data_type)

    if df.shape[0] > 1:
        suffix = "-" + str(line.Index)
        clones_df.clone_id += suffix
        bcmap_df.clone_id += suffix
        new_barcodes = []
        for bcs in bcmap_df.barcodes:
            barcodes = bcs.split(",")
            new_barcodes.append(",".join([x + suffix for x in barcodes]))
        bcmap_df.barcodes = new_barcodes

    adata.obs["batch_gex_data"] = line.gex_data
    adata.obs["batch_clones_file"] = line.clones_file

    all_data.append([clones_df, bcmap_df, adata])

new_clones_df = pd.concat([x[0] for x in all_data])
new_bcmap_df = pd.concat([x[1] for x in all_data])
if len(all_data) == 1:
    new_adata = all_data[0][2]
else:
    for x in all_data:
        x[2].var_names_make_unique()
    new_adata = all_data[0][2].concatenate(*[x[2] for x in all_data[1:]])

print(f"Writing {new_clones_df.shape[0]} clonotypes to merged clones_file {args.output_merged_clones_file}")
new_clones_df.to_csv(args.output_merged_clones_file, sep="\t", index=False)
new_bcmap_df.to_csv(args.output_merged_clones_file + ".barcode_mapping.tsv", sep="\t", index=False)

print(f"Writing merged anndata object of shape {new_adata.shape} to file {args.output_merged_gex_data}")
new_adata.write_h5ad(args.output_merged_gex_data)

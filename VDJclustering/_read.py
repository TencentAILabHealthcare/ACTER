# -*- coding: utf-8 -*-
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
import scanpy as sc

MIN_CDR3_LEN = 6
MAX_CDR3_LEN = 30


def _read_gex_data(gex_data: str, gex_data_type: str, gex_only=True) -> sc.AnnData:
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


def _get_ab_from_10x_chain(chain):
    """Returns None if the chain is not in TRA or TRB"""
    if chain in ["TRA", "TRB"]:
        return chain[2]
    else:
        return None


def _read_tcr_data(contig_annotations_csvfile):
    """Parse tcr data, only taking 'productive' tcrs
    Returns: clonotype2tcrs, clonotype2barcodes
    """
    df = pd.read_csv(contig_annotations_csvfile)
    print(f"read {df.shape[0]} lines from {contig_annotations_csvfile}")
    df["productive"] = df["productive"].astype(str)
    clonotype2barcodes = {}
    clonotype2tcrs = {}

    for l in df.itertuples():
        bc = l.barcode
        clonotype = l.raw_clonotype_id
        assert l.productive in ["None", "False", "True", "TRUE", "FALSE"]
        if clonotype == "None":
            print(f"skip line of clonotype==None: {l}")
            continue
        if clonotype not in clonotype2barcodes:
            clonotype2barcodes[clonotype] = []
        if bc in clonotype2barcodes[clonotype]:
            pass
        else:
            clonotype2barcodes[clonotype].append(bc)

        if l.productive.lower() != "true":
            print(f"skip line of productive!=TRUE: {l}")
            continue
        if l.cdr3.lower() == "none" or l.cdr3_nt.lower() == "none":
            print(f"skip line of bad_cdr3: {l}")
            continue

        chain = l.chain
        ab = _get_ab_from_10x_chain(chain)
        if ab is None:
            print(f"skip line of ab=None: {l}")
            continue
        if clonotype not in clonotype2tcrs:
            clonotype2tcrs[clonotype] = {"A": Counter(), "B": Counter()}

        vg = l.v_gene
        jg = l.j_gene

        tcr_chain = (vg, jg, l.cdr3, l.cdr3_nt.lower())
        clonotype2tcrs[clonotype][ab][tcr_chain] += int(l.umis)

    return clonotype2tcrs, clonotype2barcodes


def _make_clones_file(clonotype2tcrs, clonotype2barcodes, sample, outfile, verbose=False):
    """Make a clones file with information parsed from the 10X csv files
    organism is one of ['human','mouse']
    outfile is the name of the clones file to be created
    """

    # tmpfile = outfile+'.tmp' # a temporary intermediate file
    bc_mapfile = outfile + ".barcode_mapping.tsv"
    outmap = open(bc_mapfile, "w")
    outmap.write("clone_id\tbarcodes\n")

    outfields = "clone_id sample clone_size va_gene ja_gene va2_gene ja2_gene vb_gene jb_gene cdr3a cdr3a_nucseq cdr3a2 cdr3a2_nucseq cdr3b cdr3b_nucseq".split()
    extra_fields = "alpha_umi alpha2_umi beta_umi num_alphas num_betas".split()
    outfields += extra_fields
    out = open(outfile, "w")
    out.write("\t".join(outfields) + "\n")

    for clonotype in sorted(clonotype2tcrs.keys()):
        tcrs = clonotype2tcrs[clonotype]
        if len(tcrs["A"]) >= 1 and len(tcrs["B"]) >= 1:
            atcrs = tcrs["A"].most_common()
            btcrs = tcrs["B"].most_common()
            if len(atcrs) > 1:
                # multiple alphas, pick top umi, also store the 2nd alphas
                atcr2, atcr2_umi = atcrs[1]
            else:
                atcr2, atcr2_umi = ("", "", "", ""), 0
            if len(btcrs) > 1:
                # also possible multiple betas, picking top umi, but not store the 2nd betas
                pass
            atcr, atcr_umi = atcrs[0]
            btcr, btcr_umi = btcrs[0]
            line = {}
            line["clone_id"] = clonotype
            line["sample"] = sample
            line["clone_size"] = len(clonotype2barcodes[clonotype])
            line["va_gene"] = atcr[0]
            line["ja_gene"] = atcr[1]
            line["cdr3a"] = atcr[2]
            line["cdr3a_nucseq"] = atcr[3]
            line["alpha_umi"] = str(atcr_umi)
            line["va2_gene"] = atcr2[0]
            line["ja2_gene"] = atcr2[1]
            line["cdr3a2"] = atcr2[2]
            line["cdr3a2_nucseq"] = atcr2[3]
            line["alpha2_umi"] = str(atcr2_umi)
            line["num_alphas"] = str(len(atcrs))
            line["vb_gene"] = btcr[0]
            line["jb_gene"] = btcr[1]
            line["cdr3b"] = btcr[2]
            line["cdr3b_nucseq"] = btcr[3]
            line["beta_umi"] = str(btcr_umi)
            line["num_betas"] = str(len(btcrs))
            if len(line["cdr3a"]) < MIN_CDR3_LEN or len(line["cdr3b"]) < MIN_CDR3_LEN:
                if verbose:
                    print(
                        "Warning: skipping clonotype with short cdr3s: {} {} {}".format(
                            clonotype, line["cdr3a"], line["cdr3b"]
                        )
                    )
                continue
            if len(line["cdr3a"]) > MAX_CDR3_LEN or len(line["cdr3b"]) > MAX_CDR3_LEN:
                if verbose:
                    print(
                        "Warning: skipping clonotype with long cdr3s: {} {} {}".format(
                            clonotype, line["cdr3a"], line["cdr3b"]
                        )
                    )
                continue
            out.write("\t".join(str(line[x]) for x in outfields) + "\n")
            outmap.write("{}\t{}\n".format(clonotype, ",".join(clonotype2barcodes[clonotype])))
    out.close()
    outmap.close()
    clones = pd.read_csv(outfile, sep="\t")
    print(f"output {clones.shape[0]} clonotypes in {outfile}")


def _merge_gex_vdj_desc(sample_config, clones_folder):
    all_data = []
    samples = pd.read_csv(sample_config, sep="\t")

    for line in samples.itertuples():
        clones_file = os.path.join(clones_folder, line.sample + ".clones.tsv")
        if not os.path.exists(clones_file):
            print(f"Not found: {clones_file}")
            sys.exit()
        bcmap_file = clones_file + ".barcode_mapping.tsv"
        if not os.path.exists(bcmap_file):
            print(f"Not found: {bcmap_file}")
            sys.exit()

        clones_df = pd.read_csv(clones_file, sep="\t")
        bcmap_df = pd.read_csv(bcmap_file, sep="\t")
        adata = _read_gex_data(line.gex_data, line.gex_data_type)

        suffix = "-" + str(line.Index)
        clones_df.clone_id += suffix
        bcmap_df.clone_id += suffix
        new_barcodes = []
        for bcs in bcmap_df.barcodes:
            barcodes = bcs.split(",")
            new_barcodes.append(",".join([x + suffix for x in barcodes]))
        bcmap_df.barcodes = new_barcodes

        adata.obs["batch_gex_data"] = line.gex_data
        adata.obs["batch_clones_file"] = clones_file

        all_data.append([clones_df, bcmap_df, adata])

    new_clones_df = pd.concat([x[0] for x in all_data])
    new_bcmap_df = pd.concat([x[1] for x in all_data])

    if len(all_data) == 1:
        new_adata = all_data[0][2]
    else:
        for x in all_data:
            x[2].var_names_make_unique()
        new_adata = all_data[0][2].concatenate(*[x[2] for x in all_data[1:]])

    output_merged_clones_file = os.path.join(clones_folder, "merged_clones.tsv")
    print(f"Writing {new_clones_df.shape[0]} clonotypes to merged clones_file {output_merged_clones_file}")
    new_clones_df.to_csv(output_merged_clones_file, sep="\t", index=False)
    new_bcmap_df.to_csv(output_merged_clones_file + ".barcode_mapping.tsv", sep="\t", index=False)

    if len(all_data) > 1:
        new_adata.obs["sample"] = new_adata.obs["batch"].cat.rename_categories(
            {str(i): samples.loc[i, "sample"] for i in range(samples.shape[0])}
        )
    samples.index = samples["sample"]
    for label in samples.columns.tolist()[4:]:
        new_adata.obs[label] = [samples.loc[i, label] for i in new_adata.obs["sample"]]

    return new_adata


def _integrate_vdj_gex(out_path, adata_gex_vdj_desc):
    clones_file = os.path.join(out_path, "merged_clones.tsv")
    clones_df = pd.read_csv(clones_file, sep="\t")
    bcmap_df = pd.read_csv(clones_file + ".barcode_mapping.tsv", sep="\t")

    id2tcr = {}
    tcr2id = {}
    for line in clones_df.itertuples(index=False):
        tcra = (line.va_gene, line.ja_gene, line.cdr3a, line.cdr3a_nucseq)
        tcrb = (line.vb_gene, line.jb_gene, line.cdr3b, line.cdr3b_nucseq)
        tcr = (tcra, tcrb)
        id2tcr[line.clone_id] = tcr
        tcr2id[tcr] = line.clone_id

    barcode2tcr = {}
    for line in bcmap_df.itertuples(index=False):
        barcodes = line.barcodes.split(",")
        clone_id = line.clone_id
        tcr = id2tcr[clone_id]
        for bc in barcodes:
            barcode2tcr[bc] = tcr

    mask = [x in barcode2tcr for x in adata_gex_vdj_desc.obs.index]
    print(
        f"Reducing to the {np.sum(mask)} barcodes (out of {adata_gex_vdj_desc.shape[0]}) with paired TCR sequence data"
    )
    adata = adata_gex_vdj_desc[mask, :].copy()

    # add tcr/vdj info to adata.obs
    tcrs = [barcode2tcr[x] for x in adata.obs.index]
    tcr_keys = "va ja cdr3a cdr3a_nucseq vb jb cdr3b cdr3b_nucseq".split()
    tcr_indices = ((x, y) for x in range(2) for y in range(4))
    for tag, (i, j) in zip(tcr_keys, tcr_indices):
        adata.obs[tag] = [x[i][j] for x in tcrs]
    adata.obs["cdr3a_nucseq"] = adata.obs.cdr3a_nucseq.str.lower()
    adata.obs["cdr3b_nucseq"] = adata.obs.cdr3b_nucseq.str.lower()

    output_merged_gex_vdj = os.path.join(out_path, "merged_gex_vdj.h5ad")
    print(f"Writing merged anndata object of shape {adata.shape} to file {output_merged_gex_vdj} in h5ad format")
    adata.write_h5ad(output_merged_gex_vdj)

    return adata


def read_gex_vdj(sample_config, out_path):
    if os.path.exists(out_path) and os.listdir(out_path):
        print(f"{out_path} is NOT empty. Please clear it manually and try again.")
        sys.exit()

    print("Start making clones file")
    samples = pd.read_csv(sample_config, sep="\t")
    for line in samples.itertuples():
        clonotype2tcrs, clonotype2barcodes = _read_tcr_data(line.vdj_data)
        clones_file = os.path.join(out_path, line.sample + ".clones.tsv")
        _make_clones_file(clonotype2tcrs, clonotype2barcodes, line.sample, clones_file)

    print(f"Start merging {samples.shape[0]} samples")
    adata = _merge_gex_vdj_desc(sample_config, out_path)
    adata = _integrate_vdj_gex(out_path, adata)
    print("Success!!!")
    return adata

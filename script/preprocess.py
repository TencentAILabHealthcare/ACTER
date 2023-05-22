# -*- coding: utf-8 -*-
import argparse
import sys
import os
import pandas as pd
from collections import Counter

MIN_CDR3_LEN = 6  # otherwise tcrdist has some problems; actually we could also set this to 5 and be OK


def get_ab_from_10x_chain(chain, organism):
    """Returns None if the chain is not valid for this 'organism'"""
    if organism in ["human", "mouse"]:
        if chain in ["TRA", "TRB"]:
            return chain[2]
        else:
            return None


def read_tcr_data(organism, contig_annotations_csvfile):
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
        ab = get_ab_from_10x_chain(chain, organism)
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


def make_clones_file(clonotype2tcrs, clonotype2barcodes, outfile, verbose=False):
    """Make a clones file with information parsed from the 10X csv files
    organism is one of ['human','mouse']
    outfile is the name of the clones file to be created
    """

    # tmpfile = outfile+'.tmp' # a temporary intermediate file
    bc_mapfile = outfile + ".barcode_mapping.tsv"
    outmap = open(bc_mapfile, "w")
    outmap.write("clone_id\tbarcodes\n")

    outfields = "clone_id subject clone_size va_gene ja_gene va2_gene ja2_gene vb_gene jb_gene cdr3a cdr3a_nucseq cdr3a2 cdr3a2_nucseq cdr3b cdr3b_nucseq".split()
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
            outl = {}
            outl["clone_id"] = clonotype
            outl["subject"] = "UNK_S"
            outl["clone_size"] = len(clonotype2barcodes[clonotype])
            outl["va_gene"] = atcr[0]
            outl["ja_gene"] = atcr[1]
            outl["cdr3a"] = atcr[2]
            outl["cdr3a_nucseq"] = atcr[3]
            outl["alpha_umi"] = str(atcr_umi)
            outl["va2_gene"] = atcr2[0]
            outl["ja2_gene"] = atcr2[1]
            outl["cdr3a2"] = atcr2[2]
            outl["cdr3a2_nucseq"] = atcr2[3]
            outl["alpha2_umi"] = str(atcr2_umi)
            outl["num_alphas"] = str(len(atcrs))
            outl["vb_gene"] = btcr[0]
            outl["jb_gene"] = btcr[1]
            outl["cdr3b"] = btcr[2]
            outl["cdr3b_nucseq"] = btcr[3]
            outl["beta_umi"] = str(btcr_umi)
            outl["num_betas"] = str(len(btcrs))
            if len(outl["cdr3a"]) < MIN_CDR3_LEN or len(outl["cdr3b"]) < MIN_CDR3_LEN:
                if verbose:
                    print(
                        "Warning: skipping clonotype with short cdr3s: {} {} {}".format(
                            clonotype, outl["cdr3a"], outl["cdr3b"]
                        )
                    )
                continue
            out.write("\t".join(str(outl[x]) for x in outfields) + "\n")
            outmap.write("{}\t{}\n".format(clonotype, ",".join(clonotype2barcodes[clonotype])))
    out.close()
    outmap.close()
    clones = pd.read_csv(outfile, sep='\t')
    print(f"output {clones.shape[0]} clonotypes in {outfile}")



parser = argparse.ArgumentParser(description='Preprocess the contig annotation file to generate the clonotype file')
parser.add_argument("--organism", choices=["human", "mouse"], required=True)
parser.add_argument("--filtered_contig_annotations_csvfile", help="Generated from cellranger vdj", required=True)

args = parser.parse_args()
organism = args.organism
contig_annotations_csvfile = args.filtered_contig_annotations_csvfile
outfile = contig_annotations_csvfile + ".clones.tsv"

if not os.path.exists(contig_annotations_csvfile):
    print("filtered_contig_annotations_csvfile not exist")
    sys.exit()
clonotype2tcrs, clonotype2barcodes = read_tcr_data(organism, contig_annotations_csvfile)
make_clones_file(clonotype2tcrs, clonotype2barcodes, outfile)

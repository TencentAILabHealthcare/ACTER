# -*- coding: utf-8 -*-
import os
from itertools import combinations, permutations

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn import preprocessing


def _nw(x, y, match=1, mismatch=1, gap=1):
    nx = len(x)
    ny = len(y)
    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, -nx, nx + 1)
    F[0, :] = np.linspace(0, -ny, ny + 1)
    # Pointers to trace through an optimal aligment.
    P = np.zeros((nx + 1, ny + 1))
    P[:, 0] = 3
    P[0, :] = 4
    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if x[i] == y[j]:
                t[0] = F[i, j] + match
            else:
                t[0] = F[i, j] - mismatch
            t[1] = F[i, j + 1] - gap
            t[2] = F[i + 1, j] - gap
            tmax = np.max(t)
            F[i + 1, j + 1] = tmax
            if t[0] == tmax:
                P[i + 1, j + 1] += 2
            if t[1] == tmax:
                P[i + 1, j + 1] += 3
            if t[2] == tmax:
                P[i + 1, j + 1] += 4
    # Trace through an optimal alignment.
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i, j] in [2, 5, 6, 9]:
            rx.append(x[i - 1])
            ry.append(y[j - 1])
            i -= 1
            j -= 1
        elif P[i, j] in [3, 5, 7, 9]:
            rx.append(x[i - 1])
            ry.append("-")
            i -= 1
        elif P[i, j] in [4, 6, 7, 9]:
            rx.append("-")
            ry.append(y[j - 1])
            j -= 1
    # Reverse the strings.
    rx = "".join(rx)[::-1]
    ry = "".join(ry)[::-1]
    return rx, ry


def _af_encode(peptide):
    """Calculate Atchley factors"""
    wd = os.path.dirname(__file__)
    af_df = pd.read_csv(os.path.join(wd, "data/Atchley_factors.csv"), index_col="aa")
    res = []
    for aa in peptide:
        if aa == "-":
            res.extend([0] * 5)
            continue
        res.extend(list(af_df.loc[aa]))
    return res


def centers_corrs(adata, key):
    centers = {}
    corrs_dict = {}
    data = adata.obs.copy()

    for cluster in sorted(data[key].unique()):
        print(f"Processing {cluster}")
        corrs = []
        sub = data.loc[data[key] == cluster, :]
        for comb in combinations(sub.index, 2):
            # print(f"\t {comb}")
            i_cdr3a = sub.loc[comb[0], "cdr3a"]
            i_cdr3b = sub.loc[comb[0], "cdr3b"]
            j_cdr3a = sub.loc[comb[1], "cdr3a"]
            j_cdr3b = sub.loc[comb[1], "cdr3b"]

            if i_cdr3a == j_cdr3a and i_cdr3b == j_cdr3b:
                af1 = _af_encode(i_cdr3a + i_cdr3b)
                af2 = _af_encode(j_cdr3a + j_cdr3b)
            elif i_cdr3a == j_cdr3a:
                cdr3b = _nw(i_cdr3b, j_cdr3b)
                af1 = _af_encode(i_cdr3a + cdr3b[0])
                af2 = _af_encode(j_cdr3a + cdr3b[1])
            elif i_cdr3b == j_cdr3b:
                cdr3a = _nw(i_cdr3a, j_cdr3a)
                af1 = _af_encode(cdr3a[0] + i_cdr3b)
                af2 = _af_encode(cdr3a[1] + j_cdr3b)
            else:
                cdr3a = _nw(i_cdr3a, j_cdr3a)
                cdr3b = _nw(i_cdr3b, j_cdr3b)
                af1 = _af_encode(cdr3a[0] + cdr3b[0])
                af2 = _af_encode(cdr3a[1] + cdr3b[1])

            # corr[comb] = np.corrcoef(af1, af2)[0, 1]
            corr = pd.DataFrame([af1, af2]).T.corr().iloc[0, 1]
            corrs.append([comb[0], comb[1], corr])

        corrs_df = pd.DataFrame(corrs, columns=["bc1", "bc2", "corr"])
        corrs_dict[cluster] = corrs_df

        corrs_sum = {}
        for bc in sub.index:
            corrs_filt = corrs_df[(bc == corrs_df.bc1) | (bc == corrs_df.bc2)]
            corrs_sum[bc] = sum(corrs_filt["corr"])
        center = max(corrs_sum, key=corrs_sum.get)

        centers[cluster] = center

    return centers, corrs_dict


def centers_corrs_noalign(adata, key):
    data = adata.obs.copy()
    data["corr"] = 0
    data["center"] = None
    adata = adata[:, adata.var["feature_types"] == "Atchley Factor"]

    centers = {}
    corrs_dict = {}

    for cluster in sorted(data[key].unique()):
        print(f"Processing {cluster}")
        sub = adata[adata.obs[key] == cluster, :]
        af = pd.DataFrame(sub.X.toarray().T, columns=sub.obs.index)
        max_cdr3a = max(data["cdr3a"].str.len())
        max_cdr3b = max(data["cdr3b"].str.len())
        max_nonzero_cdr3a = max([sum(af.iloc[0 : max_cdr3a * 5, col] != 0) for col in range(af.shape[1])])
        max_nonzero_cdr3b = max([sum(af.iloc[0 : max_cdr3b * 5, col] != 0) for col in range(af.shape[1])])
        max_nonzero_af = pd.concat([af.iloc[0:max_nonzero_cdr3a, :], af.iloc[max_cdr3a:max_nonzero_cdr3b, :]], axis=0)
        mt = max_nonzero_af.corr(method="pearson")
        idx = mt.sum().idxmax()
        max_record = mt[idx]
        data.loc[max_record.index, "corr"] = max_record
        data.loc[max_record.index, "center"] = idx

        centers[cluster] = idx

    for cluster in sorted(data[key].unique()):
        data_cluster = data[data[key] == cluster]
        corrs_df = pd.DataFrame({"bc1": data_cluster["center"], "bc2": data_cluster.index, "corr": data_cluster["corr"]})
        corrs_dict[cluster] = corrs_df
    
    return centers, corrs_dict


def filter_cell(adata, key, thr, centers, corrs):
    data = adata.obs.copy()
    # centers, corrs = centers_corrs(data, key)

    for cluster in sorted(data[key].unique()):
        data.loc[data[key] == cluster, key + "_center"] = centers[cluster]
        data.loc[data[key] == cluster, key + "_corr"] = None
        for l in data.loc[data[key] == cluster, :].itertuples():
            df = corrs[cluster]
            if l.Index == centers[cluster]:
                data.loc[l.Index, key + "_corr"] = 1.0
            else:
                corr = float(
                    df.loc[
                        (
                            ((df["bc1"] == centers[cluster]) & (df["bc2"] == l.Index))
                            | ((df["bc2"] == centers[cluster]) & (df["bc1"] == l.Index))
                        ),
                        "corr",
                    ]
                )
                data.loc[l.Index, key + "_corr"] = corr

    data[key + "_corr_scaled"] = preprocessing.scale(data[key + "_corr"])
    right = norm.isf(thr)
    data_keep = data[data[key + "_corr_scaled"] >= right]

    cluster_size = data_keep[key].value_counts()
    large_cluster = list(cluster_size[cluster_size > 1].index)
    data_keep = data_keep.loc[data_keep[key].isin(large_cluster)]

    filt = list(data_keep.index)
    if filt:
        adata_filt = adata[filt, :]
    else:
        print("There is no cell left!")
        adata_filt = None
    return adata_filt


def filter_cluster(adata, key, thr, centers, corrs):
    data = adata.obs.copy()
    # centers, corrs = centers_corrs(data, key)

    for cluster in sorted(data[key].unique()):
        data.loc[data[key] == cluster, key + "_center"] = centers[cluster]
        data.loc[data[key] == cluster, key + "_corr"] = None
        for l in data.loc[data[key] == cluster, :].itertuples():
            df = corrs[cluster]
            if l.Index == centers[cluster]:
                data.loc[l.Index, key + "_corr"] = 1.0
            else:
                corr = float(
                    df.loc[
                        (
                            ((df["bc1"] == centers[cluster]) & (df["bc2"] == l.Index))
                            | ((df["bc2"] == centers[cluster]) & (df["bc1"] == l.Index))
                        ),
                        "corr",
                    ]
                )
                data.loc[l.Index, key + "_corr"] = corr

    data[key + "_corr_scaled"] = preprocessing.scale(data[key + "_corr"])
    grouped = data.groupby(key)
    cluster_mean = grouped[key + "_corr_scaled"].agg(np.mean)
    for cluster in cluster_mean.index:
        data.loc[data[key] == cluster, "cluster_mean"] = cluster_mean[cluster]

    cluster_size = data[key].value_counts()
    large_cluster = list(cluster_size[cluster_size > 1].index)
    data = data.loc[data[key].isin(large_cluster)]

    right = norm.isf(thr)
    data_keep = data[data["cluster_mean"] >= right]
    filt = list(data_keep.index)
    if filt:
        adata_filt = adata[filt, :]
    else:
        print("There is no cell left!")
        adata_filt = None
    return adata_filt



# The following two functions are obsolete.
def filter1(adata, key, thr, centers, corrs):
    data = adata.obs.copy()
    # centers, corrs = centers_corrs(data, key)

    for cluster in sorted(data[key].unique()):
        data.loc[data[key] == cluster, key + "_center"] = centers[cluster]
        data.loc[data[key] == cluster, key + "_corr"] = None
        for l in data.loc[data[key] == cluster, :].itertuples():
            df = corrs[cluster]
            if l.Index == centers[cluster]:
                data.loc[l.Index, key + "_corr"] = 1.0
            else:
                corr = float(
                    df.loc[
                        (
                            ((df["bc1"] == centers[cluster]) & (df["bc2"] == l.Index))
                            | ((df["bc2"] == centers[cluster]) & (df["bc1"] == l.Index))
                        ),
                        "corr",
                    ]
                )
                data.loc[l.Index, key + "_corr"] = corr

    data[key + "_corr_scaled"] = preprocessing.scale(data[key + "_corr"])
    right = norm.isf(thr)
    data_rm = data.loc[data[key + "_corr_scaled"] < right, :]
    filt = list(set(adata.obs.index) - set(data_rm.index))
    if filt:
        adata_filt = adata[filt, :]
    else:
        print("There is no cell left!")
        adata_filt = None
    return adata_filt

def filter2(adata, key, thr, centers, corrs):
    data = adata.obs.copy()
    # centers, corrs = centers_corrs(data, key)

    for cluster in sorted(data[key].unique()):
        data.loc[data[key] == cluster, key + "_center"] = centers[cluster]
        data.loc[data[key] == cluster, key + "_corr"] = None
        for l in data.loc[data[key] == cluster, :].itertuples():
            df = corrs[cluster]
            if l.Index == centers[cluster]:
                data.loc[l.Index, key + "_corr"] = 1.0
            else:
                corr = float(
                    df.loc[
                        (
                            ((df["bc1"] == centers[cluster]) & (df["bc2"] == l.Index))
                            | ((df["bc2"] == centers[cluster]) & (df["bc1"] == l.Index))
                        ),
                        "corr",
                    ]
                )
                data.loc[l.Index, key + "_corr"] = corr

    data[key + "_corr_scaled"] = preprocessing.scale(data[key + "_corr"])
    left = norm.ppf(thr)
    data_rm = data.loc[data[key + "_corr_scaled"] < left, :]
    filt = list(set(adata.obs.index) - set(data_rm.index))
    if filt:
        adata_filt = adata[filt, :]
    else:
        print("There is no cell left!")
        adata_filt = None
    return adata_filt

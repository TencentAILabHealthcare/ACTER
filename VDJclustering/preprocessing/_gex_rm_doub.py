# -*- coding: utf-8 -*-

import scrublet as scr


def gex_remove_doublets(adata):
    scrub = scr.Scrublet(adata.X)
    adata.obs["doublet_scores"], adata.obs["predicted_doublets"] = scrub.scrub_doublets()
    adata.obs["doublet_info"] = adata.obs["predicted_doublets"].astype(str)
    adata = adata[adata.obs["doublet_info"] == "False", :]
    print(f"Remaining {adata.n_obs} cells")
    return adata

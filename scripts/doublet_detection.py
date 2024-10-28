#!/usr/bin/env python3

import argparse
import scanpy as sc
import scrublet as scr
from pathlib import Path
from scipy.stats import median_abs_deviation

def mad_outlier(adata, metric, n_mads, upper_only=False):
    metric = adata.obs[metric]
    median = metric.median()
    mad_cutoff = n_mads * median_abs_deviation(metric)

    if upper_only:
        return (metric > median + mad_cutoff)
    
    return ((metric < median - mad_cutoff) | (metric > median + mad_cutoff))


parser = argparse.ArgumentParser()
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--n_mads', type=int, required=True)
parser.add_argument('--n_mads_mt', type=int, required=True)

args = parser.parse_args()
input, output, n_mads, n_mads_mt = args.input, args.output, args.n_mads, args.n_mads_mt

adata = sc.read_h5ad(input)

# Identify outliers using MAD
is_outlier = (
    mad_outlier(adata, "log1p_total_counts", n_mads) |
    mad_outlier(adata, "log1p_n_genes_by_counts", n_mads) |
    mad_outlier(adata, "pct_counts_in_top_20_genes", n_mads) |
    mad_outlier(adata, "pct_counts_mt", n_mads_mt, upper_only=True)
)

# Remove outliers
adata = adata[~is_outlier]  # Keep only those where is_outlier == False
adata.uns["cells_removed"] = sum(is_outlier)

# Scrublet for doublet detection
## TODO look into scrublet settings
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
adata.obs["doublet_scores"] = doublet_scores

# Remove predicted doublets
adata = adata[~predicted_doublets, :]
adata.uns["doublets_removed"] = sum(predicted_doublets)

# TODO make scrublet UMAPS?
# scrub.set_embedding("UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
# scrub.plot_embedding('UMAP', order_points=True)

adata.write_h5ad(output)
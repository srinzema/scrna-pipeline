#!/usr/bin/env python3

import argparse
import scanpy as sc
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--min_genes', type=int, required=True)

args = parser.parse_args()
input, output, min_genes = args.input, args.output, args.min_genes

adata = sc.read_10x_h5(input)

sc.pp.filter_cells(adata, min_genes=min_genes)

# Mark mitochondrial, ribosomal, and hemoglobin genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

# Use these vars to calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

# Remove some useless columns
remove = ["total_counts_mt", "log1p_total_counts_mt", "total_counts_ribo", "log1p_total_counts_ribo", "total_counts_hb", "log1p_total_counts_hb"]
adata.obs = adata.obs[[x for x in adata.obs.columns if x not in remove]]

adata.write_h5ad(output)

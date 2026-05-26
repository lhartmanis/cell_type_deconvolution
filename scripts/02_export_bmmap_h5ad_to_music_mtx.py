#!/usr/bin/env python3
"""
Export a BoneMarrowMap AnnData reference to the MTX bundle expected by run_music.R.

Output files:
  <out_prefix>_counts.mtx      sparse MatrixMarket file, genes x cells
  <out_prefix>_genes.txt       one gene per row, matching MTX rows
  <out_prefix>_barcodes.txt    one cell barcode per row, matching MTX columns
  <out_prefix>_meta.tsv        cell metadata, rownames matching barcodes

run_music.R expects the count matrix as genes x cells. AnnData stores matrices as
cells x genes, so this script exports X.T.
"""

from pathlib import Path
import argparse

import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmwrite


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export AnnData reference to MuSiC-compatible MTX bundle."
    )
    parser.add_argument("--h5ad", required=True, help="Input .h5ad reference file.")
    parser.add_argument(
        "--out-prefix",
        required=True,
        help=(
            "Output prefix. Files are written as <prefix>_counts.mtx, "
            "<prefix>_genes.txt, <prefix>_barcodes.txt, and <prefix>_meta.tsv."
        ),
    )
    parser.add_argument(
        "--celltype-col",
        default="CellType_collapse",
        help="Cell type annotation column required by run_music.R. Default: CellType_collapse.",
    )
    parser.add_argument(
        "--donor-col",
        default="Donor",
        help="Donor/sample column required by MuSiC. Default: Donor.",
    )
    parser.add_argument(
        "--extra-metadata-cols",
        nargs="*",
        default=["CellType_Broad", "CellType"],
        help="Additional metadata columns to export if present.",
    )
    parser.add_argument(
        "--use-layer",
        default="counts",
        help="Layer to use for raw counts if present. Default: counts. Falls back to adata.X.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    h5ad_path = Path(args.h5ad)
    out_prefix = Path(args.out_prefix)

    print(f"Reading AnnData reference: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"AnnData shape: {adata.n_obs} cells x {adata.n_vars} genes")

    if not adata.obs_names.is_unique:
        raise ValueError("Cell/barcode names are not unique.")
    if not adata.var_names.is_unique:
        raise ValueError("Gene names are not unique.")

    required_cols = [args.celltype_col, args.donor_col]
    missing = [col for col in required_cols if col not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Missing required metadata columns: {missing}. "
            f"Available columns include: {adata.obs.columns[:20].tolist()}"
        )

    if args.use_layer in adata.layers:
        print(f"Using adata.layers['{args.use_layer}'] as count matrix.")
        X = adata.layers[args.use_layer]
    else:
        print(
            f"Layer '{args.use_layer}' not found. Falling back to adata.X. "
            "Make sure adata.X contains raw counts."
        )
        X = adata.X

    if not sp.issparse(X):
        print("Converting dense matrix to sparse CSR matrix.")
        X = sp.csr_matrix(X)

    if X.shape[0] != adata.n_obs:
        raise ValueError(f"Matrix rows do not match cells: {X.shape[0]} vs {adata.n_obs}")
    if X.shape[1] != adata.n_vars:
        raise ValueError(f"Matrix columns do not match genes: {X.shape[1]} vs {adata.n_vars}")

    metadata_cols = [args.celltype_col, args.donor_col]
    for col in args.extra_metadata_cols:
        if col in adata.obs.columns and col not in metadata_cols:
            metadata_cols.append(col)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    counts_path = str(out_prefix) + "_counts.mtx"
    genes_path = str(out_prefix) + "_genes.txt"
    barcodes_path = str(out_prefix) + "_barcodes.txt"
    meta_path = str(out_prefix) + "_meta.tsv"

    print(f"Writing counts MTX: {counts_path}")
    print(f"Exported MTX orientation: genes x cells = {adata.n_vars} x {adata.n_obs}")
    mmwrite(counts_path, X.T)

    print(f"Writing genes: {genes_path}")
    pd.Series(adata.var_names).to_csv(genes_path, index=False, header=False)

    print(f"Writing barcodes: {barcodes_path}")
    pd.Series(adata.obs_names).to_csv(barcodes_path, index=False, header=False)

    print(f"Writing metadata: {meta_path}")
    adata.obs[metadata_cols].to_csv(meta_path, sep="\t")

    print("Done. Wrote MuSiC reference bundle with prefix:", out_prefix)
    print("Metadata columns:", metadata_cols)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Create an AnnData reference object from BoneMarrowMap count and metadata files.

Input:
  - Tab-delimited single-cell count matrix with genes as rows and cells as columns
  - Tab-delimited cell metadata with cells as rows

Output:
  - Compressed .h5ad file with cells x genes matrix and selected metadata in .obs

Notes:
  AnnData expects X to be cells x genes. The input count matrix is assumed to be
  genes x cells, so it is transposed after loading.
"""

from pathlib import Path
import argparse

import pandas as pd
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create BoneMarrowMap AnnData reference from counts and metadata."
    )
    parser.add_argument(
        "--counts",
        required=True,
        help="Path to tab-delimited count matrix. Expected format: genes x cells.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Path to tab-delimited cell metadata. Rows should correspond to cells.",
    )
    parser.add_argument(
        "--metadata-index-col",
        default="Cell",
        help="Column in metadata containing cell IDs/barcodes. Default: Cell.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output .h5ad path.",
    )
    parser.add_argument(
        "--metadata-cols",
        nargs="+",
        default=["CellType_Broad", "CellType", "CellType_collapse", "Donor"],
        help=(
            "Metadata columns to keep if present. "
            "Default: CellType_Broad CellType CellType_collapse Donor."
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()

    counts_path = Path(args.counts)
    metadata_path = Path(args.metadata)
    output_path = Path(args.output)

    print(f"Reading counts from: {counts_path}")
    print("Expected input count matrix orientation: genes x cells")

    adata = sc.read_text(
        counts_path,
        delimiter="\t",
        first_column_names=True,
    ).T

    print(f"AnnData shape after transpose: {adata.n_obs} cells x {adata.n_vars} genes")

    print(f"Reading metadata from: {metadata_path}")
    meta = pd.read_csv(metadata_path, sep="\t")

    if args.metadata_index_col in meta.columns:
        meta = meta.set_index(args.metadata_index_col)
    else:
        # If the metadata was saved with cell IDs as the row index, use the first column.
        meta = pd.read_csv(metadata_path, sep="\t", index_col=0)

    if not adata.obs_names.is_unique:
        raise ValueError("Cell names in count matrix are not unique.")

    if not adata.var_names.is_unique:
        raise ValueError("Gene names in count matrix are not unique.")

    missing_cells = adata.obs_names.difference(meta.index)
    if len(missing_cells) > 0:
        raise ValueError(
            f"{len(missing_cells)} cells from count matrix are missing in metadata. "
            f"Example: {missing_cells[:5].tolist()}"
        )

    meta = meta.loc[adata.obs_names].copy()
    keep_cols = [col for col in args.metadata_cols if col in meta.columns]

    if len(keep_cols) == 0:
        raise ValueError(
            "None of the requested metadata columns were found. "
            f"Requested: {args.metadata_cols}. "
            f"Available examples: {meta.columns[:20].tolist()}"
        )

    print(f"Keeping metadata columns: {keep_cols}")
    adata.obs = meta[keep_cols].copy()

    # Store raw counts explicitly for downstream export.
    adata.layers["counts"] = adata.X.copy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")

    print(f"Wrote: {output_path}")
    print(f"Final AnnData shape: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()

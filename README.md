# MuSiC cell type deconvolution for AML bulk RNA-seq

This repository contains a reproducible workflow for estimating cell type proportions in AML bulk RNA-seq samples using [MuSiC](https://github.com/xuranw/MuSiC) and a BoneMarrowMap-derived single-cell reference.

The workflow is designed for local execution so that sensitive patient-level bulk RNA-seq data do not need to be uploaded to external services.

## Overview

The pipeline has two parts:

1. Prepare a BoneMarrowMap single-cell reference as a MuSiC-compatible MatrixMarket bundle.
2. Run MuSiC on each bulk RNA-seq batch using `run_music.R`.

```text
BoneMarrowMap count matrix + metadata
        ↓
scripts/01_create_bmmap_h5ad.py
        ↓
BM_map_ref.h5ad
        ↓
scripts/02_export_bmmap_h5ad_to_music_mtx.py
        ↓
BM_map_ref_counts.mtx
BM_map_ref_genes.txt
BM_map_ref_barcodes.txt
BM_map_ref_meta.tsv
        ↓
scripts/run_music.R
        ↓
MuSiC cell type proportion estimates
```

## Repository structure

```text
music-cell-type-deconvolution/
├── README.md
├── environment.yml
├── LICENSE
├── .gitignore
└── scripts/
    ├── 01_create_bmmap_h5ad.py
    ├── 02_export_bmmap_h5ad_to_music_mtx.py
    ├── check_music_mtx_reference.R
    ├── install_packages.R
    └── run_music.R
```

## Installation

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate music_deconvolution
```

Then install/check R dependencies:

```bash
Rscript scripts/install_packages.R
```

The R pipeline uses:

- `optparse`
- `Matrix`
- `SingleCellExperiment`
- `S4Vectors`
- `Biobase`
- `TOAST`
- `MuSiC`
- base R `parallel`

The Python reference-preparation scripts use:

- `scanpy`
- `pandas`
- `scipy`

No `zellkonverter` dependency is required because the final R workflow reads a MatrixMarket reference bundle directly rather than reading `.h5ad` files in R.

## Input data

### Bulk RNA-seq counts

`run_music.R` expects a raw, unnormalized bulk count matrix in tab-delimited format, optionally gzipped:

```text
gene    sample_1    sample_2    sample_3
A1BG    10          0           4
A2M     50          22          31
...
```

Requirements:

- rows = genes
- columns = bulk samples
- first column = gene names
- values = raw counts

### Single-cell reference bundle

`run_music.R` expects a reference prefix passed with `-p` / `--sc_ref_prefix`.

For a prefix such as:

```text
/path/to/BM_map_ref/BM_map_ref
```

these files must exist:

```text
/path/to/BM_map_ref/BM_map_ref_counts.mtx
/path/to/BM_map_ref/BM_map_ref_genes.txt
/path/to/BM_map_ref/BM_map_ref_barcodes.txt
/path/to/BM_map_ref/BM_map_ref_meta.tsv
```

The required format is:

```text
BM_map_ref_counts.mtx      sparse MatrixMarket count matrix, genes x cells
BM_map_ref_genes.txt       one gene per line, same order as MTX rows
BM_map_ref_barcodes.txt    one cell barcode per line, same order as MTX columns
BM_map_ref_meta.tsv        cell metadata, row names matching barcodes
```

The metadata file must contain the columns used for:

- cell type labels, for example `CellType_collapse`
- donor/sample labels, for example `Donor`

## Preparing the BoneMarrowMap reference

### 1. Create an AnnData reference

The first Python script converts a tab-delimited BoneMarrowMap count matrix and cell metadata table into an `.h5ad` file.

The expected count matrix orientation is genes x cells. The script transposes it to AnnData's standard cells x genes orientation.

```bash
python scripts/01_create_bmmap_h5ad.py \
  --counts /path/to/BM_map_sc_counts.txt \
  --metadata /path/to/BM_map_meta.txt \
  --metadata-index-col Cell \
  --output /path/to/BM_map_ref/BM_map_ref.h5ad
```

By default, the script keeps these metadata columns if present:

```text
CellType_Broad
CellType
CellType_collapse
Donor
```

`CellType_collapse` can also be created before this step in a separate metadata-preparation script.

### 2. Export the AnnData object to a MuSiC-compatible MTX bundle

```bash
python scripts/02_export_bmmap_h5ad_to_music_mtx.py \
  --h5ad /path/to/BM_map_ref/BM_map_ref.h5ad \
  --out-prefix /path/to/BM_map_ref/BM_map_ref \
  --celltype-col CellType_collapse \
  --donor-col Donor
```

This writes:

```text
BM_map_ref_counts.mtx
BM_map_ref_genes.txt
BM_map_ref_barcodes.txt
BM_map_ref_meta.tsv
```

Important orientation note:

- AnnData stores `X` as cells x genes.
- `run_music.R` expects the MTX matrix as genes x cells.
- Therefore `02_export_bmmap_h5ad_to_music_mtx.py` exports `X.T`.

### 3. Check the reference bundle

```bash
Rscript scripts/check_music_mtx_reference.R /path/to/BM_map_ref/BM_map_ref
```

This verifies that:

- the MTX row count matches the number of genes
- the MTX column count matches the number of cell barcodes
- all barcodes are present in the metadata table

## Running MuSiC

The current workflow calls `run_music.R` directly. `runner.R` is not required.

Example for one bulk RNA-seq batch:

```bash
Rscript scripts/run_music.R \
  -b /path/to/sl_inex_rename.txt.gz \
  -p /path/to/BM_map_ref/BM_map_ref \
  -c CellType_collapse \
  -d Donor \
  -o outputs/SL_music_celltypecollapse_hvgmarkers_boot100 \
  --gene_mode hvg_markers \
  --n_genes 4000 \
  --bootstrap 100 \
  --bootstrap_frac 0.8 \
  --n_cores 12 \
  --checkpoint_every 5 \
  --resume \
  --seed 15
```

The same command can be repeated for separate RNA-seq batches, changing only the bulk count matrix and output prefix.

See:

```text
examples/run_sl_clinseq_bootstrap.sh
```

## Gene selection modes

`run_music.R` intersects genes between the bulk and single-cell reference, then optionally selects a subset of genes before deconvolution.

Available modes:

```text
all           use all overlapping genes
hvg           select highly variable genes across reference cell type profiles
markers       select marker-like genes based on cluster-specific expression
hvg_markers   union of hvg and marker-like genes
```

Example:

```bash
--gene_mode hvg_markers --n_genes 4000
```

For `hvg_markers`, the script takes the union of the top HVG-like and marker-like genes. The final number of selected genes can therefore exceed `--n_genes`.

## Bootstrap gene subsampling

The optional bootstrap mode assesses robustness of the cell type proportion estimates by repeatedly subsampling genes and rerunning MuSiC.

Example:

```bash
--bootstrap 100 \
--bootstrap_frac 0.8 \
--n_cores 12 \
--checkpoint_every 5 \
--resume \
--seed 15
```

This means:

- run 100 gene-subsampling iterations
- use 80% of selected genes per iteration
- use 12 parallel workers
- checkpoint progress every 5 completed iterations
- resume from checkpoint files if present
- use deterministic per-iteration seeding

To avoid CPU over-subscription, the example shell script also sets BLAS threading variables to 1:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

## Outputs

If `--bootstrap 0`, the script writes:

```text
<output>.csv
<output>.rds
```

If bootstrapping is enabled, the script writes:

```text
<output>.csv          mean estimated proportions
<output>_ci025.csv    2.5% bootstrap quantile
<output>_ci975.csv    97.5% bootstrap quantile
<output>_sd.csv       bootstrap standard deviation
<output>.rds          full result object, including all bootstrap runs
```

Checkpoint files are also written next to the output RDS while bootstrapping:

```text
<output>_reslist.rds
<output>_done_idx.rds
<output>_progress.log
```

## Methods summary

Cell type proportions in AML bulk RNA-seq samples were estimated using MuSiC, a reference-based deconvolution method. Analyses were performed separately for each RNA-seq batch using raw, unnormalized count matrices. The BoneMarrowMap single-cell reference was prepared locally from a tab-delimited genes-by-cells count matrix and matching cell metadata. The reference matrix was loaded with Scanpy, transposed to AnnData cell-by-gene format, annotated with BoneMarrowMap cell type labels and donor information, and exported to a MuSiC-compatible MatrixMarket bundle. The R pipeline reconstructed a `SingleCellExperiment` object from this bundle, intersected genes between the bulk and single-cell datasets, selected highly variable and marker-like genes, and ran MuSiC using the collapsed BoneMarrowMap cell type labels and donor IDs. Robustness was assessed using bootstrap gene subsampling.

## References

- MuSiC: Wang et al., *Nature Communications*, 2019.
- BoneMarrowMap: single-cell transcriptional atlas of human hematopoiesis.
- Scanpy: Wolf et al., *Genome Biology*, 2018.

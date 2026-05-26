suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/check_music_mtx_reference.R <reference_prefix>")
}

prefix <- args[[1]]

genes_path <- paste0(prefix, "_genes.txt")
cells_path <- paste0(prefix, "_barcodes.txt")
meta_path  <- paste0(prefix, "_meta.tsv")
mtx_path   <- paste0(prefix, "_counts.mtx")

genes <- readLines(genes_path)
cells <- readLines(cells_path)
meta <- read.csv(meta_path, sep = "\t", row.names = 1, check.names = FALSE)
X <- readMM(mtx_path)

stopifnot(nrow(X) == length(genes))
stopifnot(ncol(X) == length(cells))
stopifnot(all(cells %in% rownames(meta)))

cat("Reference bundle OK\n")
cat("Genes:", length(genes), "\n")
cat("Cells:", length(cells), "\n")
cat("Matrix:", nrow(X), "x", ncol(X), "\n")
cat("Metadata columns:\n")
print(colnames(meta))

if ("CellType_collapse" %in% colnames(meta)) {
  cat("Cell type counts:\n")
  print(sort(table(meta$CellType_collapse), decreasing = TRUE))
}

# Install R dependencies for the MuSiC cell type deconvolution pipeline.
# Tested with R 4.2 / Bioconductor 3.15.

required_cran <- c(
  "optparse",
  "remotes"
)

required_bioc <- c(
  "Matrix",
  "SingleCellExperiment",
  "S4Vectors",
  "Biobase",
  "TOAST"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(version = "3.15", ask = FALSE)

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

if (!requireNamespace("MuSiC", quietly = TRUE)) {
  remotes::install_github("xuranw/MuSiC")
}

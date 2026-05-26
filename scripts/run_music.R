suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(MuSiC)
  library(parallel)
})

# -----------------------------
# helpers
# -----------------------------
msg <- function(...) { message(...); flush.console() }

read_sc_reference_mtx <- function(prefix) {
  genes_path <- paste0(prefix, "_genes.txt")
  cells_path <- paste0(prefix, "_barcodes.txt")
  meta_path  <- paste0(prefix, "_meta.tsv")
  mtx_path   <- paste0(prefix, "_counts.mtx")

  msg("Reading sc reference files:")
  msg("  ", mtx_path)
  msg("  ", genes_path)
  msg("  ", cells_path)
  msg("  ", meta_path)

  genes <- readLines(genes_path)
  cells <- readLines(cells_path)

  meta <- read.csv(meta_path, sep="\t", row.names=1, check.names=FALSE)
  stopifnot(all(cells %in% rownames(meta)))
  meta <- meta[cells, , drop=FALSE]

  X <- readMM(mtx_path)
  X <- as(X, "dgCMatrix")  # genes x cells

  stopifnot(nrow(X) == length(genes))
  stopifnot(ncol(X) == length(cells))

  rownames(X) <- genes
  colnames(X) <- cells

  sce <- SingleCellExperiment(
    assays = list(counts = X),
    colData = DataFrame(meta)
  )
  sce
}

# Gene selection by cluster profiles (fast, sparse->dense only at genes x clusters)
select_genes <- function(sce, clusters, mode="all", n_genes=3000) {
  if (mode == "all") return(rownames(sce))

  ct <- as.character(colData(sce)[[clusters]])
  stopifnot(length(ct) == ncol(sce))

  counts <- assay(sce, "counts")
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")

  f <- factor(ct)
  msg(sprintf("  Gene selection: aggregating to clusters (n_clusters=%d)", nlevels(f)))

  # cells x clusters indicator
  i <- seq_along(f)
  j <- as.integer(f)
  M <- sparseMatrix(i=i, j=j, x=1, dims=c(length(f), nlevels(f)))

  # genes x clusters sums
  cl_sum <- counts %*% M

  # CPM normalize per cluster, log1p
  lib <- Matrix::colSums(cl_sum)
  lib[lib == 0] <- 1
  cl_cpm <- cl_sum %*% Diagonal(x = 1e4 / lib)
  cl_log <- cl_cpm
  cl_log@x <- log1p(cl_log@x)

  cl_dense <- as.matrix(cl_log)  # genes x clusters (clusters small)
  v <- apply(cl_dense, 1, var)

  hvg <- character(0)
  mk  <- character(0)

  if (mode %in% c("hvg", "hvg_markers")) {
    hvg <- names(sort(v, decreasing=TRUE))[seq_len(min(n_genes, length(v)))]
  }
  if (mode %in% c("markers", "hvg_markers")) {
    effect <- apply(cl_dense, 1, function(x) max(x) - median(x))
    mk <- names(sort(effect, decreasing=TRUE))[seq_len(min(n_genes, length(effect)))]
  }

  if (mode == "hvg") return(hvg)
  if (mode == "markers") return(mk)
  if (mode == "hvg_markers") return(union(hvg, mk))

  stop("Unknown gene_mode: ", mode)
}

run_once <- function(bulk_mtx, sce_obj, clusters, samples) {
  MuSiC::music_prop(
    bulk.mtx = bulk_mtx,
    sc.sce   = sce_obj,
    clusters = clusters,
    samples  = samples
  )$Est.prop.weighted
}

# --- threading hygiene (avoid 30 cores * BLAS threads = meltdown) ---
set_blas_threads_1 <- function() {
  Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
             VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
}

# --- checkpoint paths ---
ckpt_paths <- function(out_rds) {
  base <- sub("\\.rds$", "", out_rds, ignore.case=TRUE)
  list(
    reslist = paste0(base, "_reslist.rds"),
    done    = paste0(base, "_done_idx.rds"),
    log     = paste0(base, "_progress.log")
  )
}

# --- safe save (atomic write) ---
saveRDS_atomic <- function(obj, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, path)
}

# --- run one subsample with its own RNG stream ---
run_one_iter <- function(i, G, m_sub, bulk_use, sce_use, clusters, samples, seed0) {
  set.seed(seed0 + i)  # deterministic per-iter
  g_i <- sample(G, size=m_sub, replace=FALSE)
  MuSiC::music_prop(
    bulk.mtx = bulk_use[g_i, , drop=FALSE],
    sc.sce   = sce_use[g_i, ],
    clusters = clusters,
    samples  = samples
  )$Est.prop.weighted
}

# -----------------------------
# Options
# -----------------------------
option_list <- list(
  make_option(c("-b", "--bulk_counts"), type="character",
              help="Path to bulk counts (tsv.gz). genes x samples, rownames=gene."),
  make_option(c("-p", "--sc_ref_prefix"), type="character",
              help="Prefix to sc reference MTX bundle (expects *_counts.mtx, *_genes.txt, *_barcodes.txt, *_meta.tsv)."),
  make_option(c("-c", "--clusters"), type="character",
              help="colData column for clusters (cell types)."),
  make_option(c("-d", "--samples"), type="character",
              help="colData column for samples/donors."),
  make_option(c("-o", "--output"), type="character",
              help="Output base path (writes .csv and .rds)."),

  make_option(c("--gene_mode"), type="character", default="all",
              help="Gene selection: all | hvg | markers | hvg_markers [default %default]"),
  make_option(c("--n_genes"), type="integer", default=3000,
              help="Top N genes for hvg/markers selection [default %default]"),

  make_option(c("--bootstrap"), type="integer", default=0,
              help="Number of gene subsampling runs (0=off) [default %default]"),
  make_option(c("--bootstrap_frac"), type="double", default=0.8,
              help="Fraction of selected genes per run [default %default]"),
  make_option(c("--seed"), type="integer", default=1,
              help="Random seed [default %default]"),
  make_option(c("--n_cores"), type="integer", default=1,
            help="Number of parallel workers for bootstrap [default %default]"),
  make_option(c("--checkpoint_every"), type="integer", default=5,
            help="Write checkpoint every N completed bootstraps [default %default]"),
  make_option(c("--resume"), action="store_true", default=FALSE,
            help="Resume from checkpoint if present [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

stopifnot(!is.null(opt$bulk_counts),
          !is.null(opt$sc_ref_prefix),
          !is.null(opt$clusters),
          !is.null(opt$samples),
          !is.null(opt$output))

set.seed(opt$seed)

# -----------------------------
# Load bulk
# -----------------------------
msg("Reading bulk counts...")
bulk <- read.csv(gzfile(opt$bulk_counts), sep="\t", row.names=1, check.names=FALSE)
bulk <- as.matrix(bulk)
msg(sprintf("  bulk dim: %d genes x %d samples", nrow(bulk), ncol(bulk)))
msg(paste("  bulk genes head:", paste(head(rownames(bulk)), collapse=" ")))
msg(paste("  bulk samples head:", paste(head(colnames(bulk)), collapse=" ")))

# -----------------------------
# Load sc reference
# -----------------------------
msg("Reading sc reference (MTX)...")
sce <- read_sc_reference_mtx(opt$sc_ref_prefix)

stopifnot(opt$clusters %in% colnames(colData(sce)))
stopifnot(opt$samples  %in% colnames(colData(sce)))

msg(sprintf("  sc dim: %d genes x %d cells", nrow(sce), ncol(sce)))
msg(paste("  sc genes head:", paste(head(rownames(sce)), collapse=" ")))
msg(paste0("  cluster col: ", opt$clusters, " | sample col: ", opt$samples))
msg(paste0("  counts class: ", class(assay(sce, "counts"))[1]))

# -----------------------------
# Intersect genes
# -----------------------------
msg("Intersecting genes between bulk and sc reference...")
common <- intersect(rownames(bulk), rownames(sce))
msg(paste("  overlapping genes:", length(common)))

if (length(common) < 500) {
  stop("Too few overlapping genes (", length(common), "). Check gene naming (symbols vs Ensembl, suffixes, etc.).")
}

bulk <- bulk[common, , drop=FALSE]
sce  <- sce[common, ]
gc()

# -----------------------------
# Gene selection
# -----------------------------

msg("  cells per cluster:")
print(sort(table(colData(sce)[[opt$clusters]]), decreasing=TRUE))

if (any(is.na(colData(sce)[[opt$clusters]])))
  stop("NA values found in cluster column: ", opt$clusters)

if (any(is.na(colData(sce)[[opt$samples]])))
  stop("NA values found in sample column: ", opt$samples)

msg(paste0("Selecting genes (mode=", opt$gene_mode, ", n_genes=", opt$n_genes, ") ..."))
genes_use <- select_genes(sce, clusters=opt$clusters, mode=opt$gene_mode, n_genes=opt$n_genes)
genes_use <- intersect(genes_use, rownames(bulk))

msg(paste("  selected genes:", length(genes_use)))
if (length(genes_use) < 200) stop("Selected too few genes after intersecting with bulk.")

bulk_use <- bulk[genes_use, , drop=FALSE]
sce_use  <- sce[genes_use, ]
gc()

# -----------------------------
# Output paths
# -----------------------------
out_csv <- opt$output
if (!grepl("\\.csv$", out_csv, ignore.case=TRUE)) out_csv <- paste0(out_csv, ".csv")
out_rds <- sub("\\.csv$", ".rds", out_csv, ignore.case=TRUE)

# -----------------------------
# Run MuSiC
# -----------------------------
msg("Running MuSiC ...")
t0 <- Sys.time()

if (opt$bootstrap <= 0) {
  est <- run_once(bulk_use, sce_use, opt$clusters, opt$samples)

  write.csv(est, file=out_csv, quote=FALSE)
  saveRDS(list(
    Est.prop.weighted = est,
    gene_mode = opt$gene_mode,
    n_genes = opt$n_genes,
    genes_used = rownames(bulk_use),
    seed = opt$seed,
    bootstrap = 0
  ), file=out_rds)

} else {
  set_blas_threads_1()

  B <- opt$bootstrap
  frac <- opt$bootstrap_frac
  G <- rownames(bulk_use)

  m_sub <- max(200, floor(length(G) * frac))
  m_sub <- min(m_sub, length(G))

  n_cores <- max(1, opt$n_cores)
  ckpt <- ckpt_paths(out_rds)

  msg(sprintf("  Gene subsampling: B=%d, frac=%.2f, genes/run=%d", B, frac, m_sub))
  msg(sprintf("  Parallel workers: %d", n_cores))
  msg(sprintf("  Checkpoint every: %d", opt$checkpoint_every))
  msg(sprintf("  Resume: %s", opt$resume))

  # ---- initialize / resume ----
  res_list <- vector("list", B)
  done_idx <- integer(0)

  if (opt$resume && file.exists(ckpt$reslist) && file.exists(ckpt$done)) {
    msg("  Resuming from checkpoint...")
    res_list <- readRDS(ckpt$reslist)
    done_idx <- readRDS(ckpt$done)
    done_idx <- done_idx[!is.na(done_idx)]
    done_idx <- unique(done_idx)
    msg(sprintf("  Found %d/%d completed iterations.", length(done_idx), B))
  }

  remaining <- setdiff(seq_len(B), done_idx)
  if (length(remaining) == 0) {
    msg("  Nothing to do (all iterations already completed).")
  } else {
    msg(sprintf("  Running %d remaining iterations...", length(remaining)))
  }

  t_loop <- Sys.time()
  completed_since_ckpt <- 0L

  # ---- chunked scheduling to allow periodic checkpoints ----
  # Each chunk submits up to (n_cores * checkpoint_every) jobs.
  chunk_size <- max(n_cores, n_cores * opt$checkpoint_every)

  while (length(remaining) > 0) {
    this <- head(remaining, chunk_size)
    remaining <- setdiff(remaining, this)

    # run a chunk in parallel
    chunk_res <- mclapply(
      this,
      function(i) {
        run_one_iter(i, G, m_sub, bulk_use, sce_use, opt$clusters, opt$samples, opt$seed)
      },
      mc.cores = n_cores
    )

    # store results
    for (k in seq_along(this)) {
      i <- this[k]
      res_list[[i]] <- chunk_res[[k]]
      done_idx <- c(done_idx, i)
      completed_since_ckpt <- completed_since_ckpt + 1L
    }

    # progress / ETA
    done_n <- length(unique(done_idx))
    elapsed <- as.numeric(difftime(Sys.time(), t_loop, units="secs"))
    rate <- elapsed / max(done_n, 1)
    eta <- rate * (B - done_n)

    msg(sprintf("  completed %d/%d | elapsed %.1fs | ETA %.1fs",
                done_n, B, elapsed, eta))

    # log too (useful when screen/tmux scrollback is limited)
    cat(sprintf("[%s] completed %d/%d | elapsed %.1fs | ETA %.1fs\n",
                format(Sys.time(), "%F %T"), done_n, B, elapsed, eta),
        file=ckpt$log, append=TRUE)

    # checkpoint every N *completed* iters (across chunks)
    if (completed_since_ckpt >= opt$checkpoint_every || length(remaining) == 0) {
      msg("  Writing checkpoint...")
      saveRDS_atomic(res_list, ckpt$reslist)
      saveRDS_atomic(unique(done_idx), ckpt$done)
      completed_since_ckpt <- 0L
    }
  }

  # ---- finalize: keep only successful runs (sanity) ----
  keep <- which(vapply(res_list, function(x) !is.null(x), logical(1)))
  if (length(keep) < 2) stop("Too few completed bootstrap runs to summarize.")

  # align columns across runs
  all_ct <- Reduce(union, lapply(res_list[keep], colnames))
  res_list_aligned <- lapply(res_list[keep], function(mat) {
    mat2 <- matrix(0, nrow=nrow(mat), ncol=length(all_ct),
                   dimnames=list(rownames(mat), all_ct))
    mat2[, colnames(mat)] <- mat
    mat2
  })

  # stack into 3D array: samples x celltypes x B_done
  B_done <- length(res_list_aligned)
  arr <- array(NA_real_,
               dim=c(nrow(res_list_aligned[[1]]), length(all_ct), B_done),
               dimnames=list(rownames(res_list_aligned[[1]]), all_ct, paste0("b", keep)))
  for (i in seq_len(B_done)) arr[,,i] <- res_list_aligned[[i]]

  mean_est <- apply(arr, c(1,2), mean, na.rm=TRUE)
  lo_est   <- apply(arr, c(1,2), quantile, probs=0.025, na.rm=TRUE)
  hi_est   <- apply(arr, c(1,2), quantile, probs=0.975, na.rm=TRUE)
  sd_est   <- apply(arr, c(1,2), sd, na.rm=TRUE)

  write.csv(mean_est, file=out_csv, quote=FALSE)
  write.csv(lo_est, file=sub("\\.csv$", "_ci025.csv", out_csv, ignore.case=TRUE), quote=FALSE)
  write.csv(hi_est, file=sub("\\.csv$", "_ci975.csv", out_csv, ignore.case=TRUE), quote=FALSE)
  write.csv(sd_est, file=sub("\\.csv$", "_sd.csv", out_csv, ignore.case=TRUE), quote=FALSE)

  saveRDS(list(
    mean = mean_est, ci025 = lo_est, ci975 = hi_est, sd = sd_est,
    all_runs = arr,
    gene_mode = opt$gene_mode,
    n_genes = opt$n_genes,
    genes_used = rownames(bulk_use),
    seed = opt$seed,
    bootstrap = B,
    bootstrap_done = keep,
    bootstrap_frac = frac,
    n_cores = n_cores,
    checkpoint_files = ckpt
  ), file=out_rds)
}

msg(paste0("Done. Runtime: ", round(difftime(Sys.time(), t0, units="mins"), 2), " min"))

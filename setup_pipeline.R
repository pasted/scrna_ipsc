#!/usr/bin/env Rscript
# =============================================================================
# Setup script for iPSC-microglia scRNA-seq pipeline
# Installs core + optional CRAN/Bioconductor packages
# Added: glmGamPoi (fast SCT) and ashr (LFC shrinker fallback)
# =============================================================================

message("\n==> Using R ", getRversion(), " on ", R.version$platform)

# --- helpers ------------------------------------------------------------------
cran_repo <- getOption("repos")
if (is.null(cran_repo) || is.na(cran_repo["CRAN"]) || cran_repo["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

np <- getOption("Ncpus", 1L)
if (is.null(np) || is.na(np) || !is.numeric(np)) np <- 1L
message("Will use Ncpus = ", np)

install_if_missing_cran <- function(pkgs) {
  pkgs <- setdiff(pkgs, rownames(installed.packages()))
  if (length(pkgs)) {
    message("\nInstalling CRAN packages: ", paste(pkgs, collapse=", "))
    install.packages(pkgs, Ncpus = np, quiet = TRUE)
  } else {
    message("\nAll requested CRAN packages already installed.")
  }
}

ensure_biocmanager <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", Ncpus = np, quiet = TRUE)
  }
  suppressPackageStartupMessages(library(BiocManager))
}

install_if_missing_bioc <- function(pkgs) {
  ensure_biocmanager()
  pkgs <- setdiff(pkgs, rownames(installed.packages()))
  if (length(pkgs)) {
    message("\nInstalling Bioconductor packages: ", paste(pkgs, collapse=", "))
    BiocManager::install(pkgs, Ncpus = np, ask = FALSE, update = TRUE)
  } else {
    message("\nAll requested Bioconductor packages already installed.")
  }
}

install_doubletfinder <- function() {
  # DoubletFinder is sometimes on CRAN; if not, fall back to GitHub
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    message("\nInstalling DoubletFinder (GitHub fallback)...")
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", quiet = TRUE)
    tryCatch(
      remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = "never", Ncpus = np),
      error = function(e) {
        message("  - GitHub install failed: ", conditionMessage(e))
      }
    )
  }
}

# --- packages -----------------------------------------------------------------

cran_core <- c(
  # Seurat + friends
  "Seurat", "sctransform", "harmony",
  # tidy / util
  "Matrix", "dplyr", "purrr", "ggplot2", "ggrepel", "readr", "tidyr",
  "stringr", "tibble", "scales", "magrittr",
  # conflicts resolver
  "conflicted",         # <— added
  # stats / modeling
  "glmnet", "pROC", "DirichletReg",
  # plotting helpers possibly used by enrichplot
  "ggraph", "igraph",
  # infra
  "future",
  # msigdbr lives on CRAN
  "msigdbr",
  # LFC shrinker fallback (CRAN)
  "ashr"
)


cran_optional <- c(
  "patchwork", "cowplot", "fastDummies"
)

bioc_core <- c(
  # DE + RNA-seq
  "DESeq2", "edgeR", "limma",
  # shrinker (preferred)
  "apeglm",
  # fast GLM for SCTransform v2
  "glmGamPoi",
  # GSEA / enrichment
  "fgsea", "clusterProfiler", "enrichplot", "DOSE",
  # annotation
  "org.Hs.eg.db", "GO.db",
  # variance partitioning
  "variancePartition", "BiocParallel",
  # single-cell containers / label transfer support
  "SingleCellExperiment", "zellkonverter",
  # ambient RNA (optional in pipeline flags)
  "SoupX"
)

# --- install ------------------------------------------------------------------
install_if_missing_cran(cran_core)
install_if_missing_cran(cran_optional)
install_if_missing_bioc(bioc_core)
install_doubletfinder()

# --- quick load test & report -------------------------------------------------
pkgs_test <- c(
  cran_core,
  "DoubletFinder",
  bioc_core
)

loaded <- vapply(pkgs_test, function(p) {
  suppressWarnings(suppressPackageStartupMessages(require(p, character.only = TRUE)))
}, logical(1))

message("\nLoad summary:")
print(sort(loaded))

if (!all(loaded)) {
  not_loaded <- names(loaded)[!loaded]
  message("\nThe following packages failed to load. You can re-run the script or install them manually:\n  - ",
          paste(not_loaded, collapse = "\n  - "))
} else {
  message("\nAll requested packages loaded successfully.")
}

# --- Bioconductor sanity check (optional) -------------------------------------
if (requireNamespace("BiocManager", quietly = TRUE)) {
  v <- tryCatch(BiocManager::valid(), error = function(e) e)
  if (inherits(v, "error")) {
    message("\nBiocManager::valid() check skipped: ", conditionMessage(v))
  } else {
    message("\nBiocManager::valid() summary:")
    print(v)
  }
}

message("\nSetup complete.")
# =============================================================================

#!/usr/bin/env Rscript

# =========================
# Setup for scRNA pipeline
# - Installs CRAN, Bioconductor, and GitHub deps
# - Adds m3addon (v3 DoubletFinder helpers)
# =========================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = max(1L, parallel::detectCores() - 1L))
options(timeout = max(300, getOption("timeout", 60)))

msg <- function(...) cat(sprintf("[setup %s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse=" ")))

install_if_missing <- function(pkgs) {
  inst <- installed.packages()[, "Package"]
  need <- setdiff(pkgs, inst)
  if (length(need)) {
    msg("Installing CRAN packages: ", paste(need, collapse = ", "))
    install.packages(need, dependencies = TRUE)
  }
}

install_bioc_if_missing <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  inst <- installed.packages()[, "Package"]
  need <- setdiff(pkgs, inst)
  if (length(need)) {
    msg("Installing Bioconductor packages: ", paste(need, collapse = ", "))
    BiocManager::install(need, update = FALSE, ask = FALSE)
  }
}

# ---------------- CRAN ----------------
cran_pkgs <- c(
  "Seurat", "harmony", "Matrix", "dplyr", "purrr", "ggplot2", "ggrepel",
  "readr", "tidyr", "stringr", "glmnet", "pROC", "DirichletReg",
  "future", "conflicted", "magrittr", "remotes"
)
install_if_missing(cran_pkgs)

# ---------------- Bioconductor ----------------
bioc_pkgs <- c(
  "DESeq2", "fgsea", "org.Hs.eg.db", "clusterProfiler", "enrichplot",
  "variancePartition", "edgeR", "limma", "SingleCellExperiment",
  "zellkonverter", "GO.db", "glmGamPoi", "Biobase", "SummarizedExperiment", "batchelor"
)
install_bioc_if_missing(bioc_pkgs)

# ---------------- Optional shrinkers ----------------
# apeglm (Bioc) + ashr (CRAN)
install_bioc_if_missing("apeglm")
install_if_missing("ashr")

# ---------------- DoubletFinder + m3addon ----------------
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

#Monocle3 dependencies
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr'))

#BPCells package
remotes::install_github("bnprks/BPCells/r")

#Batchelor - required by m3addon, in bioc_pkgs above

#Monocle3 
devtools::install_github('cole-trapnell-lab/monocle3')

# Install/refresh DoubletFinder (use GitHub to avoid old mirrors)
try({
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder",
                          upgrade = "never", dependencies = TRUE)
}, silent = TRUE)

# Install m3addon (exposes paramSweep_v3/doubletFinder_v3 in some forks)
if (!requireNamespace("m3addon", quietly = TRUE)) {
  remotes::install_github("scfurl/m3addon", upgrade = "never", dependencies = TRUE)
}

# ---------------- Sanity checks ----------------
ok <- list(
  DoubletFinder = requireNamespace("DoubletFinder", quietly = TRUE),
  m3addon       = requireNamespace("m3addon", quietly = TRUE)
)
msg("DF installed? ", ok$DoubletFinder, " | m3addon installed? ", ok$m3addon)

# Probe for helper availability
probe <- function(pkg, fn) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(FALSE)
  exists(fn, envir = asNamespace(pkg), inherits = FALSE)
}
report <- c(
  DF_paramSweep_v3   = probe("DoubletFinder", "paramSweep_v3"),
  DF_doubletFinder_v3= probe("DoubletFinder", "doubletFinder_v3"),
  M3_paramSweep_v3   = probe("m3addon",       "paramSweep_v3"),
  M3_doubletFinder_v3= probe("m3addon",       "doubletFinder_v3"),
  DF_paramSweep      = probe("DoubletFinder", "paramSweep"),
  DF_doubletFinder   = probe("DoubletFinder", "doubletFinder")
)
print(report)

msg("Setup complete.")

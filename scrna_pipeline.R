#!/usr/bin/env Rscript
# =============================================================================
# iPSC-microglia scRNA-seq full pipeline
# Seurat v5: adaptive integration (try SCT, else LogNormalize + RPCA)
# Dimensionality: PCA, with optional Harmony (chip) via USE_HARMONY flag
# Pseudobulk DESeq2 -> TWO VARIANTS: labelMG & MGclusters (+ MA & Volcano each)
# Hallmark GSEA -> TWO VARIANTS
# GO enrichment -> TWO VARIANTS (H4/H8 top100 + DE up/down)
# Add-ons: Label transfer (default), DoubletFinder (opt), SoupX (opt),
#          Variance partitioning, Subset composition (DirichletReg),
#          HEK cross-chip diagnostics, Minimal marker panel (glmnet)
#          Restored: per-cluster DE CSVs + top genes figure
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(tidyr)
  library(stringr)
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(glmnet)
  library(pROC)
  library(DirichletReg)
  library(variancePartition)     # pulls limma, edgeR
  library(future)
  library(conflicted)
  suppressWarnings(requireNamespace("glmGamPoi", quietly = TRUE))
  suppressPackageStartupMessages(library(magrittr))
})

# ------------------------------ Conflicts prefs --------------------------------
conflicted::conflict_prefer("filter",    "dplyr")
conflicted::conflict_prefer("select",    "dplyr")
conflicted::conflict_prefer("lag",       "dplyr")
conflicted::conflict_prefer("count",     "dplyr")
conflicted::conflict_prefer("as.factor", "base")
conflicted::conflict_prefer("factor",    "base")
conflicted::conflict_prefer("setdiff",   "base")
conflicted::conflict_prefer("intersect", "base")
conflicted::conflict_prefer("union",     "base")

# ------------------------------ Future / memory --------------------------------
plan(sequential)
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB

# ------------------------------ Logging helper ---------------------------------
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     paste0(..., collapse=" ")))

# ------------------------------ Working directory ------------------------------
log_msg("Initial working directory:", normalizePath(getwd(), mustWork = FALSE))
args   <- commandArgs(trailingOnly = TRUE)
wd_arg <- sub("^--wd=", "", args[grepl("^--wd=", args)])
if (length(wd_arg) == 1 && dir.exists(wd_arg)) { setwd(wd_arg); log_msg("Set working directory from --wd:", normalizePath(getwd())) }

required_dirs <- c("H4_A","H4_C","H4_D","H8_A","H8_C","H8_D","H4_HEK","H8_HEK")
find_root <- function(start = getwd(), markers = required_dirs, max_up = 6){
  cand <- normalizePath(start, mustWork = FALSE)
  for (i in 0:max_up) { if (all(dir.exists(file.path(cand, markers)))) return(cand); cand <- dirname(cand) }
  NULL
}
if (!all(dir.exists(required_dirs))) {
  root <- find_root()
  if (!is.null(root)) { setwd(root); log_msg("Auto-detected project root; setwd:", normalizePath(getwd())) }
  else log_msg("Auto-detect failed. Use --wd=/path/to/project or call setwd() before running.")
}
log_msg("Final working directory:", normalizePath(getwd(), mustWork = FALSE))

# ------------------------- Project paths & output ------------------------------
PATHS <- list(
  H4_A = "./H4_A/H4_A_CogentDS_Analysis_Processed.rds",
  H4_C = "./H4_C/H4_C_CogentDS_Analysis_Processed.rds",
  H4_D = "./H4_D/H4_D_CogentDS_Analysis_Processed.rds",
  H8_A = "./H8_A/H8_A_CogentDS_analysis_Processed.rds",
  H8_C = "./H8_C/H8_C_CogentDS_analysis_Processed.rds",
  H8_D = "./H8_D/H8_D_CogentDS_analysis_Processed.rds"
)
PATHS_HEK <- list(
  H4_HEK = "./H4_HEK/H4_HEK_CogentDS_Analysis_Processed.rds",
  H8_HEK = "./H8_HEK/H8_HEK_CogentDS_Analysis_Processed.rds"
)

OUT_ROOT <- "results"; FIG_ROOT <- "figures"
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_ROOT, showWarnings = FALSE, recursive = TRUE)
make_dir <- function(path){ if(!dir.exists(path)) dir.create(path, recursive=TRUE); path }

# --------------------------- Helpers / Themes ---------------------------------
resolve_rds <- function(p) {
  if (file.exists(p) && !dir.exists(p) && grepl("\\.rds$", p, ignore.case = TRUE)) return(p)
  if (dir.exists(p)) {
    cand <- list.files(p, pattern="\\.rds$", full.names=TRUE, ignore.case=TRUE)
    if (!length(cand)) stop("No .rds file found in: ", p)
    message("Note: multiple .rds in ", p, "; using: ", basename(cand[1]))
    return(cand[1])
  }
  stop("Path not found: ", p)
}

theme_basic <- function(base=12){
  theme_minimal(base_size=base) +
    theme(panel.grid.major = element_line(linewidth=.25, colour="grey85"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face="bold"),
          strip.text = element_text(face="bold"))
}
save_plot <- function(p, file, w=9, h=7){ ggsave(file, p, width=w, height=h, dpi=300, bg="white"); log_msg("Saved:", file) }

to_entrez <- function(genes, from="SYMBOL"){
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if(!length(genes)) return(character())
  suppressMessages(suppressWarnings(
    bitr(genes, fromType=from, toType="ENTREZID", OrgDb=org.Hs.eg.db)
  )) |>
    dplyr::distinct(ENTREZID) |>
    dplyr::pull(ENTREZID)
}

run_go_and_plots <- function(gene_symbols, set_label, universe_symbols=NULL,
                             show_n=35, ont="BP"){
  subdir <- make_dir(file.path(FIG_ROOT, "GO", set_label))
  eg  <- to_entrez(gene_symbols, "SYMBOL")
  uni <- if(!is.null(universe_symbols)) to_entrez(universe_symbols, "SYMBOL") else NULL
  if(length(eg) < 10){ log_msg("GO: skipping", set_label, "mapped:", length(eg)); return(invisible(NULL)) }
  
  ego <- enrichGO(gene=eg, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ont,
                  universe=uni, pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
  ego_df <- as.data.frame(ego)
  if (!nrow(ego_df)) { log_msg("GO: none for", set_label); return(invisible(NULL)) }
  ego_sim <- tryCatch(enrichplot::pairwise_termsim(ego), error=function(e) NULL)
  
  ncat   <- min(show_n, nrow(ego_df))
  h_dot  <- max(7, min(18, 0.38 * ncat + 3))
  
  p_dot <- enrichplot::dotplot(ego, showCategory = ncat, label_format = 45) +
    ggtitle(sprintf("%s • GO:%s", set_label, ont)) +
    theme_basic(13) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
    guides(y = guide_axis(n.dodge = 2)) +
    theme(axis.text.y = element_text(size = 9, lineheight = 0.95),
          plot.margin = grid::unit(c(6, 6, 6, 16), "pt"))
  
  has_treemap <- "treemapplot" %in% getNamespaceExports("enrichplot")
  has_cp_bar  <- "barplot"     %in% getNamespaceExports("clusterProfiler")
  if (has_treemap) {
    p_tree <- enrichplot::treemapplot(ego, showCategory = ncat) +
      ggtitle(sprintf("%s • GO:%s Treemap", set_label, ont)) + theme_basic(12)
  } else if (has_cp_bar) {
    p_tree <- clusterProfiler::barplot(ego, showCategory = ncat) +
      ggtitle(sprintf("%s • GO:%s (barplot)", set_label, ont)) + theme_basic(12)
  } else {
    top_df <- ego_df |>
      dplyr::slice_head(n = ncat) |>
      dplyr::mutate(neglogFDR = -log10(p.adjust))
    p_tree <- ggplot(top_df, aes(x = reorder(Description, neglogFDR), y = neglogFDR, fill = p.adjust)) +
      geom_col() + coord_flip() +
      labs(x = NULL, y = "-log10(FDR)") +
      ggtitle(sprintf("%s • GO:%s (custom barplot)", set_label, ont)) +
      theme_basic(12)
  }
  
  p_emap <- NULL
  if (!is.null(ego_sim) && nrow(ego_df) > 1 && "emapplot" %in% getNamespaceExports("enrichplot")) {
    p_emap <- tryCatch(
      enrichplot::emapplot(ego_sim, showCategory = ncat) +
        ggtitle(sprintf("%s • GO:%s Similarity", set_label, ont)) + theme_basic(12),
      error = function(e) NULL
    )
  }
  
  save_plot(p_dot,  file.path(subdir, "dotplot_BP.png"),        w = 11, h = h_dot)
  save_plot(p_tree, file.path(subdir, "treemap_or_bar_BP.png"), w = 12, h = 7.5)
  if (!is.null(p_emap)) save_plot(p_emap, file.path(subdir, "emap_BP.png"), w = 12, h = 8)
  readr::write_csv(ego_df, file.path(subdir, "go_table_BP.csv"))
  invisible(ego)
}

# ------------------------------ Flags -----------------------------------------
have_DF   <- requireNamespace("DoubletFinder", quietly = TRUE)
have_Soup <- requireNamespace("SoupX", quietly = TRUE)
have_zell <- requireNamespace("zellkonverter", quietly = TRUE) &&
  requireNamespace("SingleCellExperiment", quietly = TRUE)
RUN_DOUBLETS <- TRUE
RUN_SOUPX    <- FALSE
RUN_BULK_OPT <- FALSE

# >>>>>>>>>> Harmony control <<<<<<<<<<
USE_HARMONY <- FALSE  # IMPORTANT: chip is confounded with harvest in this design

# --------------------------- Resolve inputs / audit ----------------------------
PATHS     <- lapply(PATHS, resolve_rds)
PATHS_HEK <- lapply(PATHS_HEK, resolve_rds)
audit <- tibble::tibble(
  sample = c(names(PATHS), names(PATHS_HEK)),
  rds    = c(unlist(PATHS), unlist(PATHS_HEK)),
  exists = file.exists(c(unlist(PATHS), unlist(PATHS_HEK)))
)
print(audit)
if (!all(audit$exists)) stop("Missing RDS file(s): ", paste(audit$sample[!audit$exists], collapse=", "))

# --------------------------- Load & annotate ----------------------------------
load_one <- function(path, nm){
  x <- readRDS(resolve_rds(path))
  if (!inherits(x, "Seurat")) stop("RDS for ", nm, " is not a Seurat object.")
  x$sample  <- nm
  x$harvest <- ifelse(grepl("^H4", nm), "H4", "H8")
  x$culture <- ifelse(grepl("_(A|C|D)$", nm), sub("^.*_", "", nm),
                      ifelse(grepl("HEK", nm), "HEK", NA))
  x$chip    <- ifelse(x$harvest=="H4", "chip_A", "chip_B")
  x
}
log_msg("Loading microglia samples..."); seur_list <- purrr::imap(PATHS, load_one)
log_msg("Loading HEK controls...");      hek_list  <- purrr::imap(PATHS_HEK, load_one)

# ---------------------- (Optional) Doublets + SoupX per sample -----------------
if (RUN_DOUBLETS && have_DF){
  log_msg("Running DoubletFinder per sample (quick settings)...")
  suppressPackageStartupMessages(library(DoubletFinder))
  seur_list <- lapply(seur_list, function(obj){
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj, verbose=FALSE) |> FindVariableFeatures(verbose=FALSE) |> ScaleData(verbose=FALSE)
    obj <- RunPCA(obj, npcs=30, verbose=FALSE) |> RunUMAP(dims=1:20, verbose=FALSE)
    sweep.res   <- paramSweep_v3(obj, PCs=1:20, sct=FALSE, num.cores=1)
    sweep.stats <- summarizeSweep(sweep.res, GT=FALSE)
    pk_tbl      <- find.pK(sweep.stats); pk <- as.numeric(as.character(pk_tbl$pK[which.max(pk_tbl$BCmetric)]))
    if (is.na(pk)) pk <- 0.01
    nExp <- round(ncol(obj)*0.05)
    obj <- doubletFinder_v3(obj, PCs=1:20, pN=0.25, pK=pk, nExp=nExp, reuse.pANN=FALSE, sct=FALSE)
    df_col <- grep("DF.classifications", colnames(obj@meta.data), value=TRUE)[1]
    obj$doublet <- obj@meta.data[[df_col]]
    obj
  })
} else if (RUN_DOUBLETS) log_msg("DoubletFinder not installed; skipping.")

if (RUN_SOUPX && have_Soup){
  log_msg("Running SoupX (placeholder: requires raw droplets). Skipping.")
}

# --------------------------- Merge --------------------------------------------
log_msg("Merging...")
obj_mg  <- Reduce(function(a,b) merge(a,b), seur_list)
obj_hek <- Reduce(function(a,b) merge(a,b), hek_list)
obj_all <- merge(obj_mg, obj_hek)
if (!"percent.mt" %in% colnames(obj_all@meta.data)) obj_all[["percent.mt"]] <- PercentageFeatureSet(obj_all, pattern="^MT-")
if ("doublet" %in% colnames(obj_mg@meta.data)){ keep <- colnames(obj_mg)[obj_mg$doublet!="Doublet"]; obj_mg <- obj_mg[, keep]; log_msg("Removed doublets; remaining:", ncol(obj_mg)) }

# --------------------------- HEK QC mini-panel --------------------------------
make_dir(file.path(FIG_ROOT, "HEK_QC"))
hek_df <- obj_hek@meta.data |> tibble::as_tibble()
save_plot(ggplot(hek_df, aes(sample, after_stat(count))) + geom_bar() + labs(y="Cells") + theme_basic(),
          file.path(FIG_ROOT, "HEK_QC", "cell_counts.png"), 7, 5)
save_plot(ggplot(hek_df, aes(sample, percent.mt)) + geom_boxplot(outlier.alpha=.2) + labs(y="Mito %") + theme_basic(),
          file.path(FIG_ROOT, "HEK_QC", "mito_boxplot.png"), 7, 5)

# --------------------------- SCTransform / LogNorm Integration -----------------
log_msg("Preparing per-sample objects for integration...")
obj_list_raw <- SplitObject(obj_mg, split.by="sample")
safe_ncol <- function(x){ nc <- tryCatch(ncol(x), error=function(e) NA_real_); if (is.na(nc)) 0L else as.integer(nc) }
cell_counts <- vapply(obj_list_raw, safe_ncol, integer(1))
obj_list_raw <- obj_list_raw[cell_counts >= 50]
if (length(obj_list_raw) < 2) stop("Fewer than two samples with >=50 cells after filtering.")
assay_names <- function(x) tryCatch(names(Assays(x)), error=function(e) character(0))

log_msg("Attempting SCT integration...")
obj_list_sct <- lapply(obj_list_raw, function(x){
  tryCatch({ DefaultAssay(x) <- "RNA"; SCTransform(x, vst.flavor="v2", return.only.var.genes=FALSE, verbose=FALSE) },
           error=function(e){ attr(x,"sct_error") <- conditionMessage(e); x })
})
sct_ok <- vapply(obj_list_sct, function(x) "SCT" %in% assay_names(x), logical(1))

if (all(sct_ok)) {
  log_msg("SCT created successfully for all samples — proceeding with SCT (RPCA) integration.")
  obj_list <- lapply(obj_list_sct, function(x){
    DefaultAssay(x) <- "SCT"
    if (length(VariableFeatures(x)) == 0) x <- FindVariableFeatures(x, selection.method="vst", nfeatures=3000, assay="SCT")
    RunPCA(x, npcs=50, verbose=FALSE)
  })
  features <- tryCatch(
    SelectIntegrationFeatures(object.list=obj_list, nfeatures=3000, assay=rep("SCT", length(obj_list))),
    error=function(e){
      message("[SelectIntegrationFeatures] retry: ", e$message)
      common <- Reduce(intersect, lapply(obj_list, \(x) rownames(x[["SCT"]])))
      if (length(common) < 2000) stop("Too few common genes across samples (", length(common), ").")
      ranked <- lapply(obj_list, \(x) VariableFeatures(x, assay="SCT")[VariableFeatures(x, assay="SCT") %in% common])
      union_ranked <- unique(unlist(ranked)); if (length(union_ranked) < 3000) head(union_ranked, max(2000, length(union_ranked))) else union_ranked[1:3000]
    }
  )
  log_msg("FindIntegrationAnchors (SCT, RPCA)...")
  anchors <- FindIntegrationAnchors(object.list=obj_list, normalization.method="SCT",
                                    anchor.features=features, reduction="rpca", dims=1:50, k.filter=50, verbose=FALSE)
  log_msg("IntegrateData (SCT, RPCA)...")
  obj_int <- IntegrateData(anchorset=anchors, dims=1:50, normalization.method="SCT", verbose=FALSE)
  int_assay  <- "integrated"; norm_method <- "SCT"
} else {
  bad <- names(obj_list_sct)[!sct_ok]
  errs <- vapply(obj_list_sct[bad], function(x) if (!is.null(attr(x,"sct_error"))) attr(x,"sct_error") else "no SCT assay created", character(1))
  log_msg("SCT not for all samples — fallback LogNormalize (RPCA). Offenders: ", paste(paste0(bad," (",errs,")"), collapse="; "))
  obj_list <- lapply(obj_list_raw, function(x){
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, normalization.method="LogNormalize", scale.factor=1e4, verbose=FALSE)
    x <- FindVariableFeatures(x, selection.method="vst", nfeatures=3000, verbose=FALSE)
    x <- ScaleData(x, verbose=FALSE)
    RunPCA(x, npcs=50, verbose=FALSE)
  })
  features <- SelectIntegrationFeatures(object.list=obj_list, nfeatures=3000)
  log_msg("FindIntegrationAnchors (LogNormalize, RPCA)...")
  anchors <- FindIntegrationAnchors(object.list=obj_list, anchor.features=features,
                                    normalization.method="LogNormalize", reduction="rpca",
                                    dims=1:50, k.filter=50, verbose=FALSE)
  log_msg("IntegrateData (LogNormalize, RPCA)...")
  obj_int <- IntegrateData(anchorset=anchors, dims=1:50, verbose=FALSE)
  int_assay  <- "integrated"; norm_method <- "LogNormalize"
}

# --------------------------- PCA + (optional) Harmony --------------------------
log_msg("PCA ...")
DefaultAssay(obj_int) <- int_assay

feat_names <- rownames(GetAssayData(obj_int, assay=int_assay, slot="data"))
stopifnot(!is.null(feat_names), length(feat_names) > 0)
vf <- VariableFeatures(obj_int)
if (is.null(vf) || length(vf) == 0) { if (exists("features") && length(features) > 0) vf <- intersect(features, feat_names) }
if (is.null(vf) || length(vf) < 500) { obj_int <- FindVariableFeatures(obj_int, selection.method="vst", nfeatures=3000, assay=int_assay); vf <- VariableFeatures(obj_int) }
VariableFeatures(obj_int) <- vf

obj_int <- ScaleData(obj_int, features=vf, verbose=FALSE)
npc <- min(50, length(vf))
obj_int <- RunPCA(obj_int, features=vf, npcs=npc, verbose=FALSE)
use_dims <- 1:min(50, ncol(obj_int@reductions$pca))

reduction_used <- "pca"
if (isTRUE(USE_HARMONY)) {
  log_msg("Running Harmony on 'chip' (note: confounded with harvest in this design)")
  rh_args <- list(object=obj_int, group.by.vars="chip", dims.use=use_dims, verbose=FALSE)
  h_formals <- names(formals(harmony::RunHarmony))
  if ("assay.use" %in% h_formals)    rh_args$assay.use <- int_assay
  if ("reduction" %in% h_formals)    rh_args$reduction <- "pca"
  if ("project.dim" %in% h_formals)  rh_args$project.dim <- FALSE
  if ("epsilon.cluster" %in% h_formals) rh_args$epsilon.cluster <- 1e-5
  obj_int <- do.call(harmony::RunHarmony, rh_args)
  reduction_used <- "harmony"
} else {
  log_msg("Harmony disabled (reduction_used = 'pca').")
}

obj_int <- RunUMAP(obj_int, reduction = reduction_used,
                   dims = 1:min(30, length(use_dims)), verbose=FALSE)

# ---------- Clustering (on chosen reduction) ----------
log_msg("Clustering on ", reduction_used, " …")
obj_int <- Seurat::FindNeighbors(obj_int, reduction = reduction_used,
                                 dims = 1:min(30, length(use_dims)), verbose = FALSE)
obj_int <- Seurat::FindClusters(obj_int, resolution = 0.5, algorithm = 1, verbose = FALSE)
Seurat::Idents(obj_int) <- "seurat_clusters"

saveRDS(obj_int, file.path(OUT_ROOT, "seurat_integrated.rds"))

# UMAPs
umap_dir <- make_dir(file.path(FIG_ROOT, "UMAP"))
lab_suffix <- if (reduction_used == "harmony") "(integrated; Harmony)" else "(integrated; PCA)"
p_umap1 <- DimPlot(obj_int, reduction="umap", group.by="harvest") + theme_basic() + ggtitle(paste("UMAP • Harvest", lab_suffix))
p_umap2 <- DimPlot(obj_int, reduction="umap", group.by="culture") + theme_basic() + ggtitle(paste("UMAP • Culture", lab_suffix))
p_umap3 <- DimPlot(obj_int, reduction="umap", group.by="chip")    + theme_basic() + ggtitle(paste("UMAP • Chip", lab_suffix))
p_umap_clusters <- Seurat::DimPlot(obj_int, reduction = "umap", group.by = "seurat_clusters",
                                   label = TRUE, repel = TRUE) + theme_basic() +
  ggplot2::ggtitle(paste("UMAP • Seurat clusters", lab_suffix))
save_plot(p_umap1, file.path(umap_dir, "umap_harvest.png"), 7.5, 6.5)
save_plot(p_umap2, file.path(umap_dir, "umap_culture.png"), 7.5, 6.5)
save_plot(p_umap3, file.path(umap_dir, "umap_chip.png"),    7.5, 6.5)
save_plot(p_umap_clusters, file.path(umap_dir, "umap_clusters.png"), 8.5, 7)

# --------------------------- Label transfer (DEFAULT) --------------------------
log_msg("Cell identity: label transfer or fallback scoring...")
label_dir <- make_dir(file.path(OUT_ROOT, "labels"))
DefaultAssay(obj_int) <- int_assay

get_ref <- function(){
  if (file.exists("refs/microglia_ref_seurat.rds")) { ref <- readRDS("refs/microglia_ref_seurat.rds"); attr(ref,"type") <- "seurat"; return(ref) }
  if (have_zell && file.exists("refs/olah_2020_microglia.h5ad")) {
    suppressPackageStartupMessages({ library(zellkonverter); library(SingleCellExperiment) })
    ref_sce <- zellkonverter::readH5AD("refs/olah_2020_microglia.h5ad")
    ref <- as.Seurat(ref_sce); attr(ref,"type") <- "h5ad"; return(ref)
  }
  NULL
}
ref <- get_ref()

if (!is.null(ref) && "celltype" %in% colnames(ref@meta.data)) {
  anchors_ref <- FindTransferAnchors(reference=ref, query=obj_int, normalization.method=norm_method)
  preds <- TransferData(anchorset=anchors_ref, refdata=ref$celltype)
  obj_int$label <- preds$predicted.id
  obj_int$label_score <- preds$prediction.score.max
  saveRDS(obj_int, file.path(OUT_ROOT, "seurat_with_labels.rds"))
  write_csv(obj_int@meta.data |> tibble::rownames_to_column("cell"),
            file.path(label_dir, "labels_with_scores.csv"))
} else {
  mg_markers <- c("CX3CR1","P2RY12","TMEM119","CSF1R","TREM2","AIF1","ITGAM")
  obj_int <- AddModuleScore(obj_int, features=list(mg_markers), name="MG")
  obj_int$MG_Score <- obj_int@meta.data$MG1
  obj_int$label <- ifelse(obj_int$MG_Score > quantile(obj_int$MG_Score, 0.4), "Microglia_like", "Non_MG_like")
  obj_int$label_score <- scales::rescale(obj_int$MG_Score)
  saveRDS(obj_int, file.path(OUT_ROOT, "seurat_with_labels_coarse.rds"))
  write_csv(obj_int@meta.data |> tibble::rownames_to_column("cell"),
            file.path(label_dir, "labels_coarse_scores.csv"))
}

# Composition stacks
label_sum <- obj_int@meta.data |>
  tibble::as_tibble() |>
  group_by(sample, harvest, label) |>
  summarise(n=n(), .groups="drop") |>
  group_by(sample, harvest) |>
  mutate(frac = n/sum(n))
make_dir(file.path(FIG_ROOT, "Celltypes"))
save_plot(
  label_sum |> ggplot(aes(sample, frac, fill=label)) + geom_col(width=0.7) +
    facet_wrap(~harvest, scales="free_x") + scale_y_continuous(labels=scales::percent) +
    labs(title="Cell type composition by sample & harvest", y="% cells", x="Sample") + theme_basic(),
  file.path(FIG_ROOT, "Celltypes", "stacks_by_harvest.png"), 10, 6
)
write_csv(label_sum, file.path(OUT_ROOT, "labels", "celltype_proportions_by_sample.csv"))

# ---------- Olah labels by cluster: dotplot + majority for UMAP ----------
stopifnot("seurat_clusters" %in% colnames(obj_int@meta.data))
stopifnot("label" %in% colnames(obj_int@meta.data))

meta_df <- obj_int@meta.data |>
  tibble::as_tibble(rownames = "cell") |>
  dplyr::mutate(seurat_clusters = as.character(.data$seurat_clusters),
                label = as.character(.data$label))

lab_counts <- meta_df |>
  dplyr::count(seurat_clusters, label, name = "n") |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::mutate(frac = n / sum(n)) |>
  dplyr::ungroup()

maj_labels <- lab_counts |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::slice_max(frac, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::transmute(seurat_clusters, maj_label = as.character(label))

p_dot_lab <- ggplot2::ggplot(lab_counts,
                             ggplot2::aes(x = seurat_clusters, y = label, size = frac, fill = frac)) +
  ggplot2::geom_point(shape = 21, alpha = 0.95, color = "grey20", stroke = 0.2) +
  ggplot2::scale_size(range = c(1.5, 8), breaks = c(0.1, 0.25, 0.5, 0.75, 1.0)) +
  ggplot2::scale_fill_gradient(low = "white", high = "firebrick") +
  theme_basic() +
  ggplot2::labs(title = "Olah label distribution by cluster",
                x = "Cluster", y = "Transferred label", size = "Frac", fill = "Frac")
celltypes_dir <- make_dir(file.path(FIG_ROOT, "Celltypes"))
save_plot(p_dot_lab, file.path(celltypes_dir, "olah_labels_by_cluster_dotplot.png"), 10, 8)

# ---------- UMAP with majority Olah label per cluster ----------
umap_mat  <- Seurat::Embeddings(obj_int, "umap")
stopifnot(!is.null(dim(umap_mat)), ncol(umap_mat) >= 2)
umap_cols <- colnames(umap_mat)[1:2]

emb <- as.data.frame(umap_mat) |>
  tibble::rownames_to_column("cell") |>
  dplyr::rename(U1 = !!umap_cols[1], U2 = !!umap_cols[2]) |>
  dplyr::left_join(meta_df[, c("cell","seurat_clusters")], by = "cell") |>
  dplyr::left_join(maj_labels, by = "seurat_clusters")

centers <- emb |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::summarise(
    U1 = stats::median(U1, na.rm = TRUE),
    U2 = stats::median(U2, na.rm = TRUE),
    maj_label = dplyr::first(maj_label),
    .groups = "drop"
  )

emb$seurat_clusters     <- base::factor(emb$seurat_clusters)
centers$seurat_clusters <- base::factor(centers$seurat_clusters, levels = levels(emb$seurat_clusters))

p_umap_maj <- ggplot2::ggplot(emb, ggplot2::aes(U1, U2, color = seurat_clusters)) +
  ggplot2::geom_point(size = 0.25, alpha = 0.65) +
  theme_basic() +
  ggplot2::guides(color = ggplot2::guide_legend(title = "Cluster",
                                                override.aes = list(size = 3, alpha = 1))) +
  ggrepel::geom_label_repel(
    data = centers,
    ggplot2::aes(label = maj_label, fill = seurat_clusters),
    color = "white", size = 3, label.size = 0.2, segment.size = 0.2,
    alpha = 0.7, show.legend = FALSE
  ) +
  ggplot2::labs(title = paste0("UMAP • Cluster majority (Olah) labels ", lab_suffix))
umap_dir <- make_dir(file.path(FIG_ROOT, "UMAP"))
save_plot(p_umap_maj, file.path(umap_dir, "olah_cluster_majority_labels_on_umap.png"), 9, 7)

# ---------- UMAP: Seurat clusters with H8/H4 composition labels ----------
emb2 <- as.data.frame(umap_mat) |>
  tibble::rownames_to_column("cell") |>
  dplyr::rename(U1 = !!umap_cols[1], U2 = !!umap_cols[2]) |>
  dplyr::left_join(
    obj_int@meta.data |> tibble::as_tibble(rownames = "cell") |>
      dplyr::select(cell, seurat_clusters, harvest),
    by = "cell"
  )

comp <- emb2 %>%
  dplyr::count(seurat_clusters, harvest, name = "n") %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup() %>%
  tidyr::complete(seurat_clusters,
                  harvest = factor(c("H4","H8"), levels = c("H4","H8")),
                  fill = list(frac = 0, n = 0)) %>%
  tidyr::pivot_wider(names_from = harvest, values_from = frac, values_fill = 0) %>%
  dplyr::mutate(H4 = dplyr::coalesce(H4, 0),
                H8 = dplyr::coalesce(H8, 0))

centers2 <- emb2 %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(U1 = stats::median(U1, na.rm = TRUE),
                   U2 = stats::median(U2, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(comp,       by = "seurat_clusters") %>%
  dplyr::left_join(maj_labels, by = "seurat_clusters") %>%
  dplyr::mutate(
    lab = ifelse(
      is.na(maj_label),
      sprintf("C%s\nH8 %d%% / H4 %d%%", seurat_clusters, round(100*H8), round(100*H4)),
      sprintf("C%s — %s\nH8 %d%% / H4 %d%%", seurat_clusters, maj_label, round(100*H8), round(100*H4))
    )
  )

emb2$seurat_clusters     <- base::factor(emb2$seurat_clusters)
centers2$seurat_clusters <- base::factor(centers2$seurat_clusters, levels = levels(emb2$seurat_clusters))

xr <- diff(range(emb2$U1, na.rm = TRUE)); yr <- diff(range(emb2$U2, na.rm = TRUE))
pad_x <- 0.15 * xr; pad_y <- 0.15 * yr

p_umap_comp <- ggplot2::ggplot(emb2, ggplot2::aes(U1, U2, color = seurat_clusters)) +
  ggplot2::geom_point(size = 0.25, alpha = 0.65) +
  theme_basic() +
  ggplot2::guides(color = ggplot2::guide_legend(title = "Cluster",
                                                override.aes = list(size = 3, alpha = 1))) +
  ggrepel::geom_label_repel(
    data = centers2,
    ggplot2::aes(label = lab, fill = seurat_clusters),
    color = "black", size = 3, label.size = 0.2,
    alpha = 0.7,
    force = 5, force_pull = 1,
    box.padding = grid::unit(0.6, "lines"),
    point.padding = grid::unit(0.2, "lines"),
    max.iter = 15000, max.time = 6,
    min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.7,
    max.overlaps = Inf, seed = 42, show.legend = FALSE
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = pad_x)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = pad_y)) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(plot.margin = grid::unit(c(8, 80, 8, 8), "pt")) +
  ggplot2::labs(title = paste0("Seurat clusters with H8/H4 composition ", lab_suffix),
                x = colnames(umap_mat)[1], y = colnames(umap_mat)[2])
save_plot(p_umap_comp, file.path(FIG_ROOT, "UMAP", "umap_clusters_with_harvest_composition.png"), 9, 7.5)

# --------------------------- Per-cluster DE (H8 vs H4) ------------------------
de_assay <- if ("SCT" %in% names(Seurat::Assays(obj_int))) "SCT" else "RNA"
Seurat::DefaultAssay(obj_int) <- de_assay
if (de_assay == "RNA") { try({ obj_int <- Seurat::NormalizeData(obj_int, assay = "RNA", verbose = FALSE) }, silent = TRUE) }
Seurat::Idents(obj_int) <- "seurat_clusters"

de_cluster_dir <- make_dir(file.path(OUT_ROOT, "DE_cluster"))

meta <- obj_int@meta.data |> tibble::as_tibble(rownames = "cell") |>
  dplyr::select(cell, seurat_clusters, harvest)
cell_tab <- meta |>
  dplyr::count(seurat_clusters, harvest, name = "n") |>
  tidyr::pivot_wider(names_from = harvest, values_from = n, values_fill = 0)
readr::write_csv(cell_tab, file.path(de_cluster_dir, "cells_per_cluster_by_harvest.csv"))

min_cells <- 10
valid_clusters <- cell_tab |>
  dplyr::filter(H4 >= min_cells, H8 >= min_cells) |>
  dplyr::pull(seurat_clusters)

if (length(valid_clusters) == 0) {
  log_msg("Per-cluster DE skipped: no clusters with >= ", min_cells, " cells in BOTH H4 and H8.")
} else {
  results_list <- list()
  top_list <- list()
  topn <- 5
  
  for (cl in valid_clusters) {
    cl_cells <- meta |> dplyr::filter(seurat_clusters == cl) |> dplyr::pull(cell)
    obj_cl <- obj_int[, cl_cells]
    Seurat::Idents(obj_cl) <- obj_cl$harvest
    
    de_res <- tryCatch(
      Seurat::FindMarkers(
        obj_cl,
        ident.1 = "H8", ident.2 = "H4",
        assay = de_assay,
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0.05,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(de_res) && nrow(de_res) > 0) {
      de_tbl <- de_res |>
        as.data.frame() |>
        tibble::rownames_to_column("gene") |>
        dplyr::mutate(cluster = as.character(cl))
      readr::write_csv(de_tbl, file.path(de_cluster_dir, paste0("cluster_", cl, "_H8_vs_H4.csv")))
      results_list[[cl]] <- de_tbl
      
      top_cl <- dplyr::bind_rows(
        de_tbl |> dplyr::arrange(dplyr::desc(avg_log2FC)) |> dplyr::slice_head(n = topn) |> dplyr::mutate(direction = "Up in H8"),
        de_tbl |> dplyr::arrange(avg_log2FC)             |> dplyr::slice_head(n = topn) |> dplyr::mutate(direction = "Down in H8")
      )
      top_list[[cl]] <- top_cl
    } else {
      log_msg("Cluster ", cl, ": no DE (insufficient cells or no expression differences).")
    }
  }
  
  full_tbl <- dplyr::bind_rows(results_list)
  top_tbl  <- dplyr::bind_rows(top_list)
  
  if (!is.null(full_tbl) && nrow(full_tbl) > 0) {
    readr::write_csv(full_tbl, file.path(de_cluster_dir, "all_markers_by_cluster_H8_vs_H4.csv"))
  }
  if (!is.null(top_tbl) && nrow(top_tbl) > 0) {
    readr::write_csv(top_tbl, file.path(de_cluster_dir, "top_updown_by_cluster_H8_vs_H4.csv"))
    gene_union <- unique(top_tbl$gene)
    Seurat::DefaultAssay(obj_int) <- de_assay
    de_fig <- make_dir(file.path(FIG_ROOT, "DE"))
    dp <- Seurat::DotPlot(obj_int, features = gene_union, group.by = "seurat_clusters") +
      ggplot2::scale_colour_gradient(low = "grey80", high = "firebrick") +
      ggplot2::scale_size(range = c(1, 6)) +
      theme_basic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = paste("Top genes by cluster (H8 vs H4)", lab_suffix), x = "Cluster", y = "Gene")
    save_plot(dp, file.path(de_fig, "top_genes_by_cluster_H8vsH4.png"), 12, 8)
  }
}

# --------------------------- Pseudobulk DESeq2 (TWO VARIANTS) ------------------
log_msg("Pseudobulk counts per (culture, harvest) ...")
DefaultAssay(obj_int) <- "RNA"
md2 <- obj_int@meta.data |> tibble::as_tibble(rownames="cell")

run_pseudobulk_de <- function(cells_subset, tag,
                              out_dir = file.path(OUT_ROOT, "deseq2"),
                              fig_dir = file.path(FIG_ROOT, "DE"),
                              alpha = 0.05, lfc_thr = 1) {
  
  if (length(cells_subset) < 100) {
    log_msg("[", tag, "] Too few cells (", length(cells_subset), "); skipping.")
    return(NULL)
  }
  
  pb <- AggregateExpression(
    obj_int[, cells_subset],
    group.by = c("culture","harvest"),
    assays   = "RNA",
    slot     = "counts"
  )$RNA
  
  coldata <- tibble(sample = colnames(pb)) |>
    tidyr::separate(sample, into = c("culture","harvest"), sep = "_",
                    remove = FALSE, extra = "merge", fill = "right") |>
    dplyr::mutate(across(c(culture, harvest), ~ base::as.factor(.x)))
  rownames(coldata) <- coldata$sample
  
  if (length(unique(coldata$harvest)) < 2) {
    log_msg("[", tag, "] Only one harvest present in pseudobulk; skipping.")
    return(NULL)
  }
  
  dds <- DESeqDataSetFromMatrix(round(pb), colData = coldata, design = ~ culture + harvest)
  dds <- DESeq(dds)
  
  coef_name <- "harvest_H8_vs_H4"
  if (!coef_name %in% resultsNames(dds)) {
    stop("[", tag, "] Coefficient '", coef_name, "' not found. Available:\n  ",
         paste(resultsNames(dds), collapse = "\n  "))
  }
  
  lfc_type <- if (requireNamespace("apeglm", quietly = TRUE)) {
    "apeglm"
  } else if (requireNamespace("ashr", quietly = TRUE)) {
    "ashr"
  } else {
    "normal"
  }
  log_msg("[", tag, "] Using LFC shrinker:", lfc_type)
  
  res <- tryCatch(
    lfcShrink(dds, coef = coef_name, type = lfc_type),
    error = function(e) {
      log_msg("[", tag, "] lfcShrink failed (", conditionMessage(e), "); using unshrunken Wald results().")
      results(dds, name = coef_name)
    }
  )
  
  res_tbl <- as.data.frame(res) |>
    tibble::rownames_to_column("gene") |>
    dplyr::arrange(padj)
  
  # Save table
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_csv <- file.path(out_dir, paste0("pseudobulk_H8_vs_H4_", tag, ".csv"))
  readr::write_csv(res_tbl, out_csv)
  
  # Plots
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  # MA
  p_ma <- ggplot(res_tbl, aes(baseMean, log2FoldChange, color = padj < alpha)) +
    geom_point(alpha = .6, size = .8) +
    scale_x_log10() +
    scale_color_manual(values = c("grey70","firebrick")) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", linewidth = .4) +
    labs(title = paste0("DESeq2 MA • H8 vs H4 [", tag, "]"),
         x = "baseMean (log10)", y = "log2FC") + theme_basic()
  save_plot(p_ma, file.path(fig_dir, paste0("MA_", tag, ".png")), 7.5, 6)
  
  # Volcano (colored + top labels)
  volc_df <- res_tbl |>
    dplyr::mutate(
      sig = dplyr::case_when(
        !is.na(padj) & padj < alpha & log2FoldChange >=  lfc_thr ~ "Up",
        !is.na(padj) & padj < alpha & log2FoldChange <= -lfc_thr ~ "Down",
        TRUE ~ "NS"
      ),
      nlp = -log10(pmax(padj, 1e-300))
    )
  
  top_lab <- volc_df |>
    dplyr::filter(sig != "NS") |>
    dplyr::arrange(dplyr::desc(nlp)) |>
    dplyr::slice_head(n = 15)
  
  p_volc <- ggplot2::ggplot(volc_df, ggplot2::aes(log2FoldChange, nlp, color = sig)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", linewidth = 0.4) +
    ggrepel::geom_text_repel(data = top_lab, ggplot2::aes(label = gene),
                             size = 2.3, max.overlaps = Inf,
                             box.padding = 0.3, min.segment.length = 0) +
    ggplot2::scale_color_manual(values = c(Down = "royalblue", NS = "grey75", Up = "firebrick")) +
    theme_basic() +
    ggplot2::labs(title = paste0("Volcano • H8 vs H4 [", tag, "]"),
                  x = "log2FC (H8/H4)", y = "-log10(FDR)", color = "")
  save_plot(p_volc, file.path(fig_dir, paste0("Volcano_labeled_", tag, ".png")), 8.5, 7)
  
  list(res_tbl = res_tbl, coldata = coldata, pb = pb)
}

# Variant 1: ORIGINAL (label-based microglia)
mg_labels <- c("Microglia","microglia","Microglia_like")
cells_labelMG <- md2 |> dplyr::filter(label %in% mg_labels) |> dplyr::pull(cell)
if (length(cells_labelMG) < 100) {
  log_msg("[labelMG] Few MG-labelled cells; falling back to culture A/C/D regardless of label.")
  cells_labelMG <- md2 |> dplyr::filter(!is.na(culture) & culture %in% c("A","C","D")) |> dplyr::pull(cell)
}
pb1 <- run_pseudobulk_de(cells_labelMG, tag = "labelMG")

# Variant 2: MGclusters (clusters where ≥50% cells are microglia-labelled)
cl_frac <- md2 |>
  dplyr::mutate(is_mg = label %in% mg_labels) |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::summarise(frac_mg = mean(is_mg), n = dplyr::n(), .groups = "drop")
mg_clusters <- cl_frac |> dplyr::filter(frac_mg >= 0.5) |> dplyr::pull(seurat_clusters)
if (length(mg_clusters) == 0) {
  log_msg("[MGclusters] No clusters with ≥50% microglia_like; skipping this variant.")
  pb2 <- NULL
} else {
  cells_MGclusters <- md2 |> dplyr::filter(seurat_clusters %in% mg_clusters) |> dplyr::pull(cell)
  pb2 <- run_pseudobulk_de(cells_MGclusters, tag = "MGclusters")
  de_meta_dir <- make_dir(file.path(OUT_ROOT, "deseq2"))
  readr::write_csv(cl_frac, file.path(de_meta_dir, "cluster_microglia_fraction.csv"))
  writeLines(paste(sort(unique(mg_clusters)), collapse = ", "),
             con = file.path(de_meta_dir, "clusters_used_in_MGclusters.txt"))
}

# Keep ORIGINAL as default for any downstream step that expects a single table
if (!is.null(pb1)) {
  res_tbl_default <- pb1$res_tbl
} else if (!is.null(pb2)) {
  log_msg("Original labelMG failed; falling back to MGclusters for downstream default.")
  res_tbl_default <- pb2$res_tbl
} else {
  stop("No valid pseudobulk result available for downstream analysis.")
}

# --------------------------- Hallmark GSEA (TWO VARIANTS) ----------------------
log_msg("Hallmark GSEA (fgsea) for two variants ...")
msig <- msigdbr(species="Homo sapiens", category="H"); stopifnot(nrow(msig) > 0)
hallmark <- split(msig$gene_symbol, msig$gs_name); hallmark <- lapply(hallmark, unique)

run_gsea_variant <- function(res_tbl, tag,
                             out_dir = file.path(OUT_ROOT, "gsea"),
                             fig_dir = file.path(FIG_ROOT, "GSEA")) {
  if (is.null(res_tbl)) { log_msg("[GSEA/", tag, "] No res_tbl; skipping."); return(NULL) }
  res_rank <- res_tbl %>%
    mutate(padj = ifelse(is.na(padj), 1, padj),
           log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
           rank = sign(log2FoldChange) * -log10(pmax(padj, 1e-300))) %>%
    filter(is.finite(rank)) %>%
    group_by(gene) %>% summarise(rank = mean(rank), .groups="drop")
  genes_in_sets <- unique(unlist(hallmark, use.names = FALSE))
  res_rank <- res_rank %>% filter(gene %in% genes_in_sets)
  ranks <- res_rank$rank; names(ranks) <- res_rank$gene; ranks <- sort(ranks, decreasing=TRUE)
  
  if (length(ranks) < 50) {
    log_msg("[GSEA/", tag, "] Too few ranked genes (", length(ranks), "); skipping.")
    return(NULL)
  }
  fg <- fgsea(pathways=hallmark, stats=ranks, minSize=15, maxSize=500, nperm=5000)
  fg_tbl <- arrange(as_tibble(fg), padj)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(fg_tbl, file.path(out_dir, paste0("hallmark_fgsea_", tag, ".csv")))
  p_gsea <- fg_tbl %>% slice_min(padj, n=25) %>%
    ggplot(aes(reorder(pathway, NES), NES, size=size, fill=padj)) +
    geom_point(shape=21, alpha=.9, stroke=.2) + coord_flip() +
    scale_fill_continuous(trans="log10", name="FDR") +
    labs(title=paste("Hallmark GSEA • H8 vs H4 [", tag, "]", lab_suffix), x="Pathway", y="NES") + theme_basic()
  save_plot(p_gsea, file.path(fig_dir, paste0("hallmark_dotplot_", tag, ".png")), 8.5, 7)
  invisible(fg_tbl)
}

invisible(run_gsea_variant(if (!is.null(pb1)) pb1$res_tbl else NULL, "labelMG"))
invisible(run_gsea_variant(if (!is.null(pb2)) pb2$res_tbl else NULL, "MGclusters"))

# --------------------------- GO enrichment (TWO VARIANTS) ----------------------
log_msg("GO enrichment: two variants (top100 H4/H8 + DE up/down)")

run_go_variant <- function(cells_subset, res_tbl, tag,
                           alpha = 0.05, lfc_thr = 1) {
  if (is.null(res_tbl) || length(cells_subset) < 50) {
    log_msg("[GO/", tag, "] Not enough cells or missing res_tbl; skipping.")
    return(invisible(NULL))
  }
  DefaultAssay(obj_int) <- "RNA"
  avg_by_harvest <- tryCatch(
    AverageExpression(obj_int[, cells_subset], group.by="harvest", assays="RNA", slot="data")$RNA,
    error = function(e) NULL
  )
  if (is.null(avg_by_harvest) || !all(c("H4","H8") %in% colnames(avg_by_harvest))) {
    log_msg("[GO/", tag, "] AverageExpression failed or missing H4/H8; skipping top100.")
    genes_H4 <- genes_H8 <- character()
  } else {
    genes_H4 <- head(rownames(avg_by_harvest[order(avg_by_harvest[,"H4"], decreasing=TRUE), , drop=FALSE]), 100)
    genes_H8 <- head(rownames(avg_by_harvest[order(avg_by_harvest[,"H8"], decreasing=TRUE), , drop=FALSE]), 100)
  }
  
  up_genes   <- res_tbl |> dplyr::filter(!is.na(padj), padj < alpha, log2FoldChange >=  lfc_thr) |> dplyr::pull(gene) |> unique()
  down_genes <- res_tbl |> dplyr::filter(!is.na(padj), padj < alpha, log2FoldChange <= -lfc_thr) |> dplyr::pull(gene) |> unique()
  universe_symbols <- rownames(obj_int)
  
  if (length(genes_H4))   invisible(run_go_and_plots(genes_H4,   paste0("variant_", tag, "/H4_top100"),    universe_symbols))
  if (length(genes_H8))   invisible(run_go_and_plots(genes_H8,   paste0("variant_", tag, "/H8_top100"),    universe_symbols))
  if (length(up_genes))   invisible(run_go_and_plots(up_genes,   paste0("variant_", tag, "/DE_up_in_H8"),  universe_symbols))
  if (length(down_genes)) invisible(run_go_and_plots(down_genes, paste0("variant_", tag, "/DE_down_in_H8"),universe_symbols))
}

# Run GO for both variants (only if the variant existed)
if (!is.null(pb1)) run_go_variant(cells_labelMG,   pb1$res_tbl, "labelMG")
if (!is.null(pb2)) run_go_variant(cells_MGclusters, pb2$res_tbl, "MGclusters")

# --------------------------- Variance Partitioning -----------------------------
log_msg("Variance partitioning: culture vs harvest (pseudobulk voom)...")
# Use the default pseudobulk matrix if available; rebuild from labelMG cells if not
if (!is.null(pb1) && !is.null(pb1$pb)) {
  pb_mat <- pb1$pb
  coldata <- tibble(sample=colnames(pb_mat)) |>
    tidyr::separate(sample, into=c("culture","harvest"), sep="_", remove=FALSE, extra="merge", fill="right") |>
    dplyr::mutate(across(c(culture,harvest), ~ base::as.factor(.x)))
  rownames(coldata) <- coldata$sample
} else {
  # fallback
  DefaultAssay(obj_int) <- "RNA"
  md <- obj_int@meta.data |> tibble::as_tibble(rownames="cell")
  cells_fallback <- md |> dplyr::filter(label %in% mg_labels) |> dplyr::pull(cell)
  pb_mat <- AggregateExpression(obj_int[, cells_fallback], group.by=c("culture","harvest"), assays="RNA", slot="counts")$RNA
  coldata <- tibble(sample=colnames(pb_mat)) |>
    tidyr::separate(sample, into=c("culture","harvest"), sep="_", remove=FALSE, extra="merge", fill="right") |>
    dplyr::mutate(across(c(culture,harvest), ~ base::as.factor(.x)))
  rownames(coldata) <- coldata$sample
}

library(edgeR); library(limma)
dge <- DGEList(counts=round(pb_mat))
keep <- filterByExpr(dge, group=coldata$harvest)
dge <- dge[keep,, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
v   <- voom(dge, design = model.matrix(~ culture + harvest, data=coldata))

vp   <- fitExtractVarPartModel(v$E, ~ culture + harvest, coldata)
vp_tbl <- as.data.frame(vp)
var_out <- make_dir(file.path(OUT_ROOT, "variance"))
write_csv(tibble::rownames_to_column(vp_tbl, "gene"), file.path(var_out, "variance_partition.csv"))

vp_long <- vp_tbl |> tibble::rownames_to_column("gene") |> tidyr::pivot_longer(-gene, names_to="factor", values_to="variance")
p_vp <- vp_long |>
  group_by(factor) |>
  summarise(median_var = median(variance, na.rm=TRUE)) |>
  ggplot(aes(reorder(factor, median_var), median_var)) +
  geom_col(width=.7) + coord_flip() + scale_y_continuous(labels=scales::percent) +
  labs(title="Median variance explained per factor", x="", y="Median % variance") + theme_basic()
save_plot(p_vp, file.path(FIG_ROOT, "DE", "variance_partition.png"), 7.5, 6)

# --------------------------- Subset composition (Dirichlet) --------------------
log_msg("Subset composition test (Dirichlet regression)...")
subset_prop <- obj_int@meta.data |>
  tibble::as_tibble() |>
  dplyr::filter(!is.na(culture), culture %in% c("A","C","D")) |>
  dplyr::count(sample, culture, harvest, label, name = "n") |>
  dplyr::group_by(sample, culture, harvest) |>
  dplyr::mutate(total = sum(n), prop = n/total) |>
  dplyr::ungroup()

comp_wide <- subset_prop |>
  dplyr::select(sample, culture, harvest, label, prop) |>
  tidyr::pivot_wider(names_from = label, values_from = prop, values_fill = 0)
comp_cols <- base::setdiff(colnames(comp_wide), c("sample","culture","harvest"))

comp_mat <- as.matrix(comp_wide[, comp_cols])
comp_mat[is.na(comp_mat)] <- 0
comp_mat[comp_mat == 0] <- 1e-6
rs <- rowSums(comp_mat)
comp_mat <- sweep(comp_mat, 1, rs, "/")

comp_wide$Y <- DirichletReg::DR_data(comp_mat)
comp_wide$harvest <- base::factor(comp_wide$harvest, levels = c("H4","H8"))
fit_dir <- DirichletReg::DirichReg(Y ~ harvest, data = comp_wide)

comp_out <- make_dir(file.path(OUT_ROOT, "composition"))
sink(file.path(comp_out, "dirichlet_summary.txt")); print(summary(fit_dir)); sink()

top_labels <- subset_prop |>
  dplyr::group_by(label) |>
  dplyr::summarise(mean_prop = mean(prop)) |>
  dplyr::arrange(dplyr::desc(mean_prop)) |>
  dplyr::slice_head(n = 8) |>
  dplyr::pull(label)

comp_plot <- subset_prop |>
  dplyr::filter(label %in% top_labels) |>
  dplyr::group_by(harvest, label) |>
  dplyr::summarise(prop = mean(prop), .groups = "drop") |>
  ggplot2::ggplot(ggplot2::aes(harvest, prop, group = label, color = label)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Average subset proportions by harvest", y = "% cells", x = "") +
  theme_basic()
comp_fig <- make_dir(file.path(FIG_ROOT, "Composition"))
save_plot(comp_plot, file.path(comp_fig, "subset_change_lines.png"), 9, 6)

# --------------------------- HEK cross-chip diagnostics ------------------------
log_msg("HEK cross-chip diagnostics...")
DefaultAssay(obj_hek) <- "RNA"
hek_pb <- AggregateExpression(
  obj_hek, group.by = c("harvest","chip"), assays = "RNA", slot = "counts"
)$RNA
hek_pb_mat <- as.matrix(hek_pb)

if (!is.null(dim(hek_pb_mat)) && ncol(hek_pb_mat) >= 2) {
  if (is.null(colnames(hek_pb_mat))) { colnames(hek_pb_mat) <- paste0("grp", seq_len(ncol(hek_pb_mat))) }
  else { colnames(hek_pb_mat) <- make.unique(colnames(hek_pb_mat), sep = "_") }
  cor_mat <- stats::cor(log1p(hek_pb_mat), method = "pearson", use = "pairwise.complete.obs")
  cor_tbl <- cor_mat |>
    as.data.frame() |>
    tibble::rownames_to_column("sample1") |>
    tidyr::pivot_longer(-sample1, names_to = "sample2", values_to = "r")
  hek_out <- make_dir(file.path(OUT_ROOT, "HEK_QC"))
  readr::write_csv(cor_tbl, file.path(hek_out, "hek_crosschip_correlation.csv"))
  p_cor <- ggplot2::ggplot(cor_tbl, ggplot2::aes(sample1, sample2, fill = r)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 3) +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
    theme_basic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "HEK cross-chip log1p(count) correlation", x = "", y = "")
  save_plot(p_cor, file.path(FIG_ROOT, "HEK_QC", "hek_crosschip_cor.png"), 7.5, 6.5)
} else {
  log_msg("Not enough HEK pseudo-bulk columns for correlation plot.")
}

# --------------------------- Minimal marker panel (glmnet) ---------------------
log_msg("Minimal marker panel via glmnet (per-cell, stratified by culture)...")
DefaultAssay(obj_int) <- int_assay
hvf <- VariableFeatures(obj_int)
if (length(hvf) < 1000) { obj_int <- FindVariableFeatures(obj_int, nfeatures=3000, assay=int_assay); hvf <- VariableFeatures(obj_int) }
X <- t(GetAssayData(obj_int, slot="data")[hvf, ])
cul <- ifelse(is.na(obj_int$culture), "UNK", as.character(obj_int$culture))
y <- factor(obj_int$harvest, levels=c("H4","H8"))  # H8 positive class
foldid <- as.integer(base::factor(cul))
cvfit <- cv.glmnet(X, y, family="binomial", standardize=TRUE, type.measure="auc", foldid=foldid)
fit  <- glmnet(X, y, family="binomial", lambda=cvfit$lambda.min, standardize=TRUE)
coefs <- as.matrix(coef(fit))
sel_genes <- base::setdiff(rownames(coefs)[coefs[, 1] != 0], "(Intercept)")
panel_tbl <- tibble::tibble(gene=sel_genes) |> dplyr::mutate(coef=as.numeric(coefs[match(gene, rownames(coefs)),1])) |> dplyr::arrange(dplyr::desc(abs(coef)))
mark_out <- make_dir(file.path(OUT_ROOT, "markers")); write_csv(panel_tbl, file.path(mark_out, "minimal_panel_glmnet.csv"))

prob <- predict(fit, newx=X, type="response")[,1]
roc_obj <- roc(response=y, predictor=as.numeric(prob), levels=c("H4","H8"), direction="<")
auc_val <- as.numeric(auc(roc_obj))
p_roc <- ggplot(data.frame(tpr=roc_obj$sensitivities, fpr=1-roc_obj$specificities), aes(fpr, tpr)) +
  geom_path(linewidth=1) + geom_abline(linetype="dashed") + coord_equal() + theme_basic() +
  labs(title=paste0("ROC (glmnet panel) — AUC=", sprintf("%.3f", auc_val)), x="False positive rate", y="True positive rate")
markers_fig <- make_dir(file.path(FIG_ROOT, "markers"))
save_plot(p_roc, file.path(markers_fig, "roc_glmnet.png"), 6.5, 6.5)

# --------------------------- OPTIONAL: Bulk RNA-seq ----------------------------
if (RUN_BULK_OPT){ log_msg("Bulk RNA-seq integration enabled, but no data configured. Skipping.") }

log_msg("Pipeline complete. Outputs in 'results/' and 'figures/'.")
# =============================================================================

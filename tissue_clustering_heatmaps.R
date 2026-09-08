# ==============================================================================
# Tissue clustering heatmaps for GTEx and PaxDb expression/abundance data
#
# R translation of Sections 7-8 of filtering_source_exploration.ipynb:
#   7. Within-source hierarchical clustering of tissues (Spearman correlation)
#   8. Top 50 most variable genes/proteins, hierarchically clustered
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(ggplot2)
  library(patchwork)
})
 
REPO_ROOT     <- normalizePath(".")
EXPR_FILE     <- file.path(REPO_ROOT, "data", "expression_by_tissue.gct")
PAXDB_FILE    <- file.path(REPO_ROOT, "data", "paxdb_all_tissues.tsv")
UBERON_GROUPS_FILE <- file.path(REPO_ROOT, "data", "tissue_uberon_groups.tsv")
OUT_ROOT      <- file.path(REPO_ROOT, "sanity_check_outputs")
GTEX_OUT      <- file.path(OUT_ROOT, "gtex")
PAXDB_OUT     <- file.path(OUT_ROOT, "paxdb")
dir.create(OUT_ROOT, showWarnings = FALSE)
dir.create(GTEX_OUT, showWarnings = FALSE)
dir.create(PAXDB_OUT, showWarnings = FALSE)

# ── Shared color scales (declared once, reused by every heatmap below) ───────

VALUE_PALETTE <- colorRampPalette(brewer.pal(9, "YlGnBu"))          # raw/log10 value heatmaps
CORR_PALETTE  <- colorRampPalette(c("#3B4CC0", "white", "#B40426")) # correlation heatmaps
SOURCE_COLORS <- c(PaxDb = "#FF7F0E", GTEx = "#1F77B4")             # GTEx vs. PaxDb annotation

# ── Data loading ──────────────────────────────────────────────────────────────

# GTEx GCT: header on line 3, genes x tissues, median TPM. Duplicate gene
# symbols are collapsed with the row-wise median (matches the Python notebook).
read_gtex_median_tpm <- function(path) {
  raw <- read_tsv(path, skip = 2, show_col_types = FALSE) %>%
    rename(ensembl = Name, symbol = Description)
  tissue_cols <- setdiff(colnames(raw), c("ensembl", "symbol"))

  raw %>%
    group_by(symbol) %>%
    summarise(across(all_of(tissue_cols), \(x) median(x, na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("symbol") %>%
    as.matrix()
}

# PaxDb wide TSV: id + one column per tissue, ppm. Missing entries mean
# "not detected" (~0 abundance), not "unknown", so fill rather than drop.
read_paxdb <- function(path) {
  mat <- read_tsv(path, show_col_types = FALSE) %>%
    column_to_rownames("id") %>%
    as.matrix()
  mat[is.na(mat)] <- 0
  mat
}

gtex_mat  <- read_gtex_median_tpm(EXPR_FILE)
paxdb_mat <- read_paxdb(PAXDB_FILE)

gtex_log  <- log2(gtex_mat + 1)
paxdb_log <- log2(paxdb_mat + 1)

# ── Tissue grouping ───────────────────────────────────────────────────────────
# Tissue -> group assignments come from group_tissues_uberon.py, which climbs
# each tissue's UBERON is_a/part_of ancestry to the nearest organ-level term
# (e.g. every GTEx "Brain - ..." subregion and PaxDb's BRAIN/CEREBRAL_CORTEX/
# FRONTAL_CORTEX all collapse to a single "brain" group), then climbs once
# more to the nearest organ-system-level term (e.g. "brain" and "spinal cord"
# both collapse to "nervous system"). Both tiers are shown as separate
# annotation strips below. Run that script to regenerate
# data/tissue_uberon_groups.tsv after the underlying tissue lists change.

uberon_groups <- read_tsv(UBERON_GROUPS_FILE, show_col_types = FALSE)

build_lookup <- function(source_name, level_col) {
  rows <- filter(uberon_groups, source == source_name)
  setNames(rows[[level_col]], rows$tissue)
}

TISSUE_GROUP_LOOKUPS <- list(
  GTEx  = list(organ = build_lookup("GTEx", "organ_group"), system = build_lookup("GTEx", "organ_system")),
  PaxDB = list(organ = build_lookup("PaxDB", "organ_group"), system = build_lookup("PaxDB", "organ_system"))
)

# Tissues absent from the UBERON grouping (e.g. PaxDb body fluids like
# AMNIOTIC_FLUID that aren't organs) fall back to their own name as a
# singleton group, so they still get a distinct label/color instead of NA.
tissue_group_name <- function(tissue_names, source = c("GTEx", "PaxDB"), level = c("organ", "system")) {
  source <- match.arg(source)
  level <- match.arg(level)
  lookup <- TISSUE_GROUP_LOOKUPS[[source]][[level]]
  group <- unname(lookup[tissue_names])
  ifelse(is.na(group), tissue_names, group)
}

# Two fixed group -> color mappings (one per level), each built once from the
# union of GTEx and PaxDb groups at that level, so the same group (e.g.
# "brain", "nervous system") always gets the same color across every heatmap
# instead of a per-call, per-source assignment. The two levels use distinct
# palette families so their annotation strips/legends stay visually separable.
all_organ_groups  <- sort(unique(c(
  tissue_group_name(colnames(gtex_mat), "GTEx", "organ"),
  tissue_group_name(colnames(paxdb_mat), "PaxDB", "organ")
)))
all_system_groups <- sort(unique(c(
  tissue_group_name(colnames(gtex_mat), "GTEx", "system"),
  tissue_group_name(colnames(paxdb_mat), "PaxDB", "system")
)))
organ_group_colors <- setNames(
  colorRampPalette(c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3")))(length(all_organ_groups)),
  all_organ_groups
)
SYSTEM_PALETTE <- c(
 "#E58606","#5D69B1","#52BCA3","#99C945","#CC61B0","#24796C","#DAA51B", "#2F8AC4" ,"#764E9F","#ED645A","#CC3A8E","#A5AA99"
)
# Some PaxDb sample types (biofluids, whole organism) have no UBERON organ
# mapping and fall back to their own name as a singleton "system" group (see
# tissue_group_name() above). There are always exactly 12 real organ systems,
# matching SYSTEM_PALETTE 1:1, so those get the fixed palette directly rather
# than through colorRampPalette() interpolation, which - once fallback groups
# are mixed in - can generate two interpolated colors close enough to be
# visually indistinguishable (e.g. cardiovascular vs. renal system). Fallback
# groups get their own grey ramp so they never compete with organ-system hues.
canonical_systems  <- sort(unique(uberon_groups$organ_system))
fallback_groups    <- setdiff(all_system_groups, canonical_systems)
system_group_colors <- c(
  setNames(SYSTEM_PALETTE[seq_along(canonical_systems)], canonical_systems),
  setNames(colorRampPalette(c("grey35", "grey75"))(length(fallback_groups)), fallback_groups)
)

# Single-column annotation data frame keyed by tissue name ("Organ system" =
# system-level), for pheatmap's annotation_row/annotation_col. Heatmaps only
# show the organ-system tier; the organ-level tier is still used by the PCA
# plots in Section 9 via tissue_group_name() directly.
tissue_annotation_df <- function(tissue_names, source) {
  data.frame(
    "Organ system" = tissue_group_name(tissue_names, source, "system"),
    row.names = tissue_names,
    check.names = FALSE
  )
}

# Subset of the global system-level color map restricted to the groups
# actually present among tissue_names. pheatmap draws a legend entry for
# every color it is handed, so passing the full palette to a GTEx-only (or
# PaxDb-only) heatmap would list groups that never appear in it.
tissue_annotation_colors <- function(tissue_names, source) {
  list(
    "Organ system" = system_group_colors[sort(unique(tissue_group_name(tissue_names, source, "system")))]
  )
}

# Draw a pheatmap result to a PNG with "Tissue" / row-count axis titles added
# (pheatmap has no built-in xlab/ylab), matching the axis labels used in the
# Python exploration notebook / playground.ipynb heatmaps.
save_heatmap_with_axis_labels <- function(pheatmap_result, out_path, width, height,
                                           xlab = NULL, ylab = NULL) {
  png(out_path, width = width, height = height, units = "in", res = 200)
  grid.newpage()
  grid.draw(pheatmap_result$gtable)
  if (!is.null(xlab)) {
    grid.text(xlab, x = 0.5, y = 0.01, gp = gpar(fontsize = 12))
  }
  if (!is.null(ylab)) {
    grid.text(ylab, x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 12))
  }
  dev.off()
}


# ── Section 8: top 50 most variable genes/proteins, hierarchically clustered ──

top_variable_heatmap <- function(mat, title, out_file, ann_col, ann_colors, value_label,
                                  out_dir = OUT_ROOT, exclude_mito = TRUE, n_top = 50, log_norm = TRUE) {
  data <- mat
  data[is.na(data)] <- 0
  if (exclude_mito) {
    data <- data[!startsWith(rownames(data), "MT-"), , drop = FALSE]
  }

  row_var <- apply(data, 1, var)
  top <- data[order(row_var, decreasing = TRUE)[seq_len(min(n_top, nrow(data)))], , drop = FALSE]

  if (log_norm) {
    max_val <- max(top)
    breaks <- unique(c(0, 10^seq(0, ceiling(log10(max_val)))))
    breaks[length(breaks)] <- max_val
    colors <- VALUE_PALETTE(length(breaks) - 1)
    legend_breaks <- breaks
    legend_labels <- formatC(breaks, format = "fg", big.mark = ",")
  } else {
    breaks <- seq(min(top), max(top), length.out = 101)
    colors <- VALUE_PALETTE(100)
    legend_breaks <- pretty(breaks)
    legend_labels <- legend_breaks
  }

  result <- pheatmap(
    top,
    color = colors,
    breaks = breaks,
    legend_breaks = legend_breaks,
    legend_labels = legend_labels,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "average",
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    annotation_legend = TRUE,
    annotation_names_col = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 6,
    fontsize_col = 6,
    angle_col = 45,
    treeheight_row = 80,
    treeheight_col = 120,
    main = title,
    silent = TRUE
  )

  row_label <- if (exclude_mito) {
    sprintf("Top %d non-mitochondrial genes/proteins by variance", n_top)
  } else {
    sprintf("Top %d genes/proteins by variance", n_top)
  }

  save_heatmap_with_axis_labels(
    result, file.path(out_dir, out_file), width = 24, height = 14,
    xlab = "Tissue", ylab = row_label
  )
}

top_variable_heatmap(
  gtex_mat,
  "GTEx - top 50 most variable genes across tissues (raw median TPM)",
  "gtex_top50_variable_clustermap_R.png",
  ann_col = tissue_annotation_df(colnames(gtex_mat), "GTEx"),
  ann_colors = tissue_annotation_colors(colnames(gtex_mat), "GTEx"),
  value_label = "TPM",
  out_dir = GTEX_OUT
)

top_variable_heatmap(
  paxdb_mat,
  "PaxDb - top 50 most variable proteins across tissues (raw ppm)",
  "paxdb_top50_variable_clustermap_R.png",
  ann_col = tissue_annotation_df(colnames(paxdb_mat), "PaxDB"),
  ann_colors = tissue_annotation_colors(colnames(paxdb_mat), "PaxDB"),
  value_label = "ppm",
  out_dir = PAXDB_OUT
)

# Same two heatmaps, but on log2(x + 1)-transformed values: top-50 variance
# genes/proteins are reselected on the log2 matrix (not the raw top-50 list),
# and plotted on a linear color scale since the values are already log2.
top_variable_heatmap(
  gtex_log,
  "GTEx - top 50 most variable genes across tissues (log2 TPM + 1)",
  "gtex_top50_variable_clustermap_log_R.png",
  ann_col = tissue_annotation_df(colnames(gtex_log), "GTEx"),
  ann_colors = tissue_annotation_colors(colnames(gtex_log), "GTEx"),
  value_label = "log2(TPM + 1)",
  log_norm = FALSE,
  out_dir = GTEX_OUT
)

top_variable_heatmap(
  paxdb_log,
  "PaxDb - top 50 most variable proteins across tissues (log2 ppm + 1)",
  "paxdb_top50_variable_clustermap_log_R.png",
  ann_col = tissue_annotation_df(colnames(paxdb_log), "PaxDB"),
  ann_colors = tissue_annotation_colors(colnames(paxdb_log), "PaxDB"),
  value_label = "log2(ppm + 1)",
  log_norm = FALSE,
  out_dir = PAXDB_OUT
)

# Combined GTEx + PaxDb: join on shared genes (all tissues, prefixed by
# source), log2-transform, then reselect the top-50 variable rows within the
# joined log2 matrix (not a union of the two individual top-50 lists).
gtex_prefixed  <- gtex_mat
paxdb_prefixed <- paxdb_mat
colnames(gtex_prefixed)  <- paste0("GTEx_", colnames(gtex_prefixed))
colnames(paxdb_prefixed) <- paste0("PaxDb_", colnames(paxdb_prefixed))

shared_genes  <- intersect(rownames(paxdb_prefixed), rownames(gtex_prefixed))
combined_mat  <- cbind(paxdb_prefixed[shared_genes, ], gtex_prefixed[shared_genes, ])
combined_log  <- log2(combined_mat + 1)

source_of <- ifelse(startsWith(colnames(combined_log), "PaxDb_"), "PaxDb", "GTEx")
source_ann_col <- data.frame(
  Source = source_of, row.names = colnames(combined_log), check.names = FALSE
)

top_variable_heatmap(
  combined_log,
  "GTEx + PaxDb combined - top 50 most variable genes/proteins (log2 value + 1)",
  "combined_top50_variable_clustermap_R.png",
  ann_col = source_ann_col,
  ann_colors = list(Source = SOURCE_COLORS),
  value_label = "log2(value + 1)",
  log_norm = FALSE
)

# ── Section 9: PCA of tissue profiles, colored by tissue group / organ system ──
# R equivalent of Section 6 of filtering_source_exploration.ipynb (PCA of each
# source's log-transformed tissue profiles), but colored with the same
# two-tier UBERON tissue taxonomy used for the heatmap annotation strips above
# (tissue_group_name()/organ_group_colors/system_group_colors) instead of that
# notebook's ad-hoc leading-token grouping, so tissue groupings are consistent
# across every plot in this script. For each source, both colorings of the
# same PCA are placed side by side via patchwork.

# One row per tissue (observations), one column per gene/protein (features);
# matches the log2 matrices already built above (gtex_log, paxdb_log).
pca_scores_for_source <- function(log_mat, source_name) {
  mat <- log_mat
  mat[is.na(mat)] <- 0
  data <- t(mat)
  keep <- apply(data, 2, function(col) sd(col) > 0)
  data <- data[, keep, drop = FALSE]

  pca <- prcomp(data, center = TRUE, scale. = FALSE)
  var_explained <- (pca$sdev^2 / sum(pca$sdev^2))[1:2]

  scores <- as.data.frame(pca$x[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  scores$tissue <- rownames(scores)
  scores[["Tissue group"]] <- tissue_group_name(scores$tissue, source_name, "organ")
  scores[["Organ system"]] <- tissue_group_name(scores$tissue, source_name, "system")
  attr(scores, "var_explained") <- var_explained
  scores
}

plot_pca_panel <- function(scores, color_col, color_map, subtitle) {
  var_explained <- attr(scores, "var_explained")
  ggplot(scores, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.4) +
    geom_point(size = 2.8, alpha = 0.9) +
    scale_color_manual(values = color_map, name = color_col) +
    labs(
      x = sprintf("PC1 (%.1f%% variance)", var_explained[1] * 100),
      y = sprintf("PC2 (%.1f%% variance)", var_explained[2] * 100),
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.35, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8, face = "bold")
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.6)))
}

plot_source_pca_pair <- function(log_mat, source_name, title, out_file, out_dir) {
  scores <- pca_scores_for_source(log_mat, source_name)

  present_organ  <- sort(unique(scores[["Tissue group"]]))
  present_system <- sort(unique(scores[["Organ system"]]))

  p_organ  <- plot_pca_panel(scores, "Tissue group", organ_group_colors[present_organ], "Colored by tissue group")
  p_system <- plot_pca_panel(scores, "Organ system", system_group_colors[present_system], "Colored by organ system")

  combined <- (p_organ | p_system) +
    plot_annotation(title = title, theme = theme(plot.title = element_text(size = 16, face = "bold")))

  ggsave(file.path(out_dir, out_file), combined, width = 16, height = 6.5, dpi = 200)
  combined
}

plot_source_pca_pair(
  gtex_log, "GTEx",
  "PCA of GTEx median TPM by tissue",
  "gtex_pca_tissue_groups_R.png",
  out_dir = GTEX_OUT
)

plot_source_pca_pair(
  paxdb_log, "PaxDB",
  "PCA of PaxDb ppm by tissue",
  "paxdb_pca_tissue_groups_R.png",
  out_dir = PAXDB_OUT
)

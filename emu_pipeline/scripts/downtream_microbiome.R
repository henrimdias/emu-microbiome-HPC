# ============================================================
# Microbiome downstream analyses (rarefaction, composition, alpha, beta)
# - Robust I/O (metadata, raw counts, relative abundance)
# - Rarefaction curves (raw counts)
# - Taxonomic composition at chosen rank (re-aggregate + renormalize)
# - Alpha diversity: S_obs (richness), Chao1 (raw), Shannon, Simpson, Pielou
#   + global tests (ANOVA/Kruskal) + post-hoc (Tukey/Dunn) + CSV exports
# - Beta diversity: Bray, Jaccard, Aitchison (CLR) distances
#   + PCoA (safe fallback), NMDS, PERMANOVA, betadisper
# - Consistent figure export helper (PNG + SVG) with explicit sizes per plot
# - Preview plots are printed before saving
# ============================================================

# -------------------------
# Packages
# -------------------------
pkgs <- c(
  "data.table","dplyr","tidyr","ggplot2","forcats","vegan","stringr","ggrepel",
  "ggpubr","rstatix","broom","tools","svglite","ape"
)
inst <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(inst)) install.packages(inst, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# -------------------------
# CONFIG (edit for your data)
# -------------------------
metadata_path   <- "/Users/westriveragresearch/Library/CloudStorage/OneDrive-SouthDakotaStateUniversity-SDSU/Postdoctoral_Project/Project_06_Protocol_emu/testing_R/metadata_native_grass.txt"                 # metadata table (samples in rows)
raw_counts_path <- "/Users/westriveragresearch/Library/CloudStorage/OneDrive-SouthDakotaStateUniversity-SDSU/Postdoctoral_Project/Project_06_Protocol_emu/testing_R/ng_raw_counts-emu-combined-genus-kraken.tsv"      # features x samples raw counts
rel_abund_path  <- "/Users/westriveragresearch/Library/CloudStorage/OneDrive-SouthDakotaStateUniversity-SDSU/Postdoctoral_Project/Project_06_Protocol_emu/testing_R/ng_emu-combined-genus-kraken.tsv"       # optional (otherwise computed)
sample_id_col   <- "SampleID"                              # sample ID column in metadata
group_col       <- "Group"                                 # grouping column used for colors/stats
tax_level       <- "genus"                                 # composition level (case-insensitive)
topN_taxa       <- 16
raref_steps     <- 100

# Output directory
EXPORT_DIR <- "/Users/westriveragresearch/Library/CloudStorage/OneDrive-SouthDakotaStateUniversity-SDSU/Postdoctoral_Project/Project_06_Protocol_emu/testing_R/testando"
dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

# (Optional) Alpha preprocessing (filters/rarefaction on RAW counts)
filter_for_alpha <- FALSE
min_lib_alpha    <- 1000     # minimum reads per sample for alpha
min_prevalence   <- 2        # minimum presence across samples to keep a feature
min_total_count  <- 10       # minimum total counts per feature
rarefy_for_alpha <- FALSE
rare_depth_alpha <- NA       # if NA, uses min library size after filtering

# (Optional) Beta minimum reads per sample (on RAW)
min_lib_beta     <- 1        # 1 removes only zero-library samples

# ============================================================
# Helpers
# ============================================================

# Robust reader for csv/tsv/txt with auto separator
smart_read <- function(path) {
  data.table::fread(path, sep = "auto", header = TRUE, data.table = FALSE, check.names = FALSE)
}

# Save a ggplot to PNG and/or SVG with explicit sizes per call (no global defaults)
save_plot_multi <- function(plot, path_noext,
                            png = TRUE, svg = TRUE,
                            png_width, png_height, png_units = c("in","cm","mm"), png_dpi = 300,
                            svg_width, svg_height, svg_units = c("cm","in","mm")) {
  dir.create(dirname(path_noext), showWarnings = FALSE, recursive = TRUE)
  
  if (isTRUE(png)) {
    if (missing(png_width) || missing(png_height))
      stop("For PNG, provide 'png_width' and 'png_height'.")
    png_units <- match.arg(png_units)
    ggplot2::ggsave(paste0(path_noext, ".png"), plot,
                    width = png_width, height = png_height, units = png_units, dpi = png_dpi)
  }
  if (isTRUE(svg)) {
    if (missing(svg_width) || missing(svg_height))
      stop("For SVG, provide 'svg_width' and 'svg_height'.")
    svg_units <- match.arg(svg_units)
    ggplot2::ggsave(paste0(path_noext, ".svg"), plot,
                    device = svglite::svglite,
                    width = svg_width, height = svg_height, units = svg_units)
  }
}

# Detect taxonomic columns (case-insensitive; treat superkingdom/domain as kingdom)
get_tax_cols <- function(df) {
  ranks_all <- c("kingdom","phylum","class","order","family","genus","species","superkingdom","domain")
  ln <- tolower(colnames(df))
  keep_idx <- which(ln %in% ranks_all)
  if (!length(keep_idx)) return(character(0))
  cols <- colnames(df)[keep_idx]
  canon_order <- c("kingdom","phylum","class","order","family","genus","species")
  ln_keep <- tolower(cols)
  canon_key <- ifelse(ln_keep %in% c("superkingdom","domain"), "kingdom", ln_keep)
  cols[order(match(canon_key, canon_order), na.last = TRUE)]
}

# Resolve desired taxonomic level: try exact; otherwise climb up, then down hierarchy
resolve_tax_level <- function(df, level) {
  ranks <- c("kingdom","phylum","class","order","family","genus","species")
  cn_raw <- colnames(df); cn <- tolower(trimws(cn_raw))
  level <- tolower(trimws(level))
  if (level %in% c("superkingdom","domain")) level <- "kingdom"
  
  hit <- which(cn == level)
  if (length(hit)) return(cn_raw[hit[1]])
  
  pos <- match(level, ranks)
  if (!is.na(pos)) {
    for (i in seq(pos-1, 1, by = -1)) { hit <- which(cn == ranks[i]); if (length(hit)) return(cn_raw[hit[1]]) }
    for (i in seq(pos+1, length(ranks), by = 1)) { hit <- which(cn == ranks[i]); if (length(hit)) return(cn_raw[hit[1]]) }
  }
  tc <- get_tax_cols(df)
  if (length(tc)) return(tc[length(tc)])
  NA_character_
}

# Chao1 (simple; with f2=0 correction)
chao1_vec <- function(x) {
  x <- as.numeric(x)
  Sobs <- sum(x > 0, na.rm = TRUE)
  f1 <- sum(x == 1, na.rm = TRUE)
  f2 <- sum(x == 2, na.rm = TRUE)
  if (f2 == 0) return(Sobs + ifelse(f1 > 0, (f1*(f1-1))/(2*(f2+1)), 0))
  Sobs + (f1^2)/(2*f2)
}

# Rarefaction (integer counts; skips zero-library samples)
compute_rarefaction_df <- function(counts_mat, steps = 60) {
  if (nrow(counts_mat) == 0) return(data.frame(SampleID=character(), Reads=numeric(), Richness=numeric()))
  out <- lapply(rownames(counts_mat), function(s) {
    v <- as.numeric(counts_mat[s, ])
    v[!is.finite(v) | v < 0] <- 0
    v <- round(v)
    n <- sum(v, na.rm = TRUE)
    if (!is.finite(n) || n <= 0) return(data.frame(SampleID = s, Reads = 0, Richness = 0))
    m_seq <- unique(round(seq(1, n, length.out = steps)))
    m_seq <- m_seq[m_seq >= 1 & m_seq <= n]
    if (length(m_seq) == 0) m_seq <- n
    rich <- suppressWarnings(sapply(m_seq, function(m) vegan::rarefy(v, sample = m)))
    data.frame(SampleID = s, Reads = m_seq, Richness = as.numeric(rich))
  })
  do.call(rbind, out)
}

# Beta helpers
percent_from_eig <- function(eig) {
  p <- eig / sum(abs(eig))
  round(100 * p[1:2], 1)
}
.clr <- function(x, pseudo = 1e-6) {
  x <- sweep(x, 1, pmax(rowSums(x), 1e-12), "/")
  x[x <= 0] <- pseudo
  l <- log(x)
  sweep(l, 1, rowMeans(l), "-")
}
pairwise_permanova <- function(dist_obj, grouping, permutations = 999) {
  gr <- factor(grouping)
  lev <- levels(gr)
  res <- data.frame(group1=character(), group2=character(), F=numeric(), R2=numeric(), p=numeric())
  if (length(lev) < 2) return(res)
  for (i in 1:(length(lev)-1)) {
    for (j in (i+1):length(lev)) {
      idx <- gr %in% c(lev[i], lev[j])
      sub_meta <- data.frame(g = droplevels(gr[idx]))
      ad <- vegan::adonis2(as.dist(as.matrix(dist_obj)[idx, idx]) ~ g,
                           data = sub_meta, permutations = permutations, by = "terms")
      res <- rbind(res,
                   data.frame(group1=lev[i], group2=lev[j],
                              F=ad$F[1], R2=ad$R2[1], p=ad$`Pr(>F)`[1]))
    }
  }
  res$p_adj <- p.adjust(res$p, method = "holm")
  res[order(res$p_adj), ]
}

# ============================================================
# 1) Read & prepare data
# ============================================================
meta <- smart_read(metadata_path) %>% as.data.frame()
stopifnot(sample_id_col %in% colnames(meta))
meta[[sample_id_col]] <- as.character(meta[[sample_id_col]])
if (!group_col %in% colnames(meta)) {
  meta[[group_col]] <- "All"
  message("`", group_col, "` not found in metadata; set to 'All'.")
}
sample_ids <- unique(meta[[sample_id_col]])

counts_raw <- smart_read(raw_counts_path) %>% as.data.frame()

# Find sample columns by intersecting with metadata IDs
sample_cols <- intersect(colnames(counts_raw), sample_ids)
if (length(sample_cols) < 2) {
  stop("No sample columns in raw counts matched the metadata column `", sample_id_col, "`.")
}

# Feature ID column (if any)
feature_col <- setdiff(colnames(counts_raw), sample_cols)
feature_col <- if (length(feature_col)) feature_col[1] else "FeatureID"
if (!feature_col %in% colnames(counts_raw)) {
  counts_raw[[feature_col]] <- paste0("feat_", seq_len(nrow(counts_raw)))
}

# Build counts matrix: samples in rows
counts_mat <- counts_raw[, sample_cols, drop = FALSE] |> as.matrix() |> t()
mode(counts_mat) <- "numeric"
rownames(counts_mat) <- sample_cols
counts_mat[!is.finite(counts_mat)] <- NA
counts_mat[is.na(counts_mat)] <- 0
counts_mat[counts_mat < 0] <- 0

# Align with metadata
common_ids <- intersect(rownames(counts_mat), sample_ids)
counts_mat  <- counts_mat[common_ids, , drop = FALSE]
meta        <- meta[match(common_ids, meta[[sample_id_col]]), , drop = FALSE]
stopifnot(all(rownames(counts_mat) == meta[[sample_id_col]]))

# Relative abundance table (read or compute)
rel_abund_df <- if (file.exists(rel_abund_path)) {
  smart_read(rel_abund_path) %>% as.data.frame()
} else {
  rel <- counts_raw
  rel[, sample_cols] <- apply(rel[, sample_cols, drop = FALSE], 2, function(x) {
    s <- sum(x, na.rm = TRUE); if (s == 0) return(rep(0, length(x))); x / s
  })
  rel
}

# ============================================================
# 2) Rarefaction curves (raw counts)
# ============================================================
lib_sizes <- rowSums(counts_mat)
counts_mat_rare <- counts_mat[lib_sizes > 0, , drop = FALSE]

rare_df <- compute_rarefaction_df(counts_mat_rare, steps = raref_steps) %>%
  dplyr::left_join(meta, by = setNames(sample_id_col, "SampleID"))

p_rare <- ggplot(rare_df, aes(x = Reads, y = Richness, color = .data[[group_col]], group = SampleID)) +
  geom_line(alpha = 0.9, linewidth = 0.7) +
  labs(x = "Subsampled reads", y = "Expected richness (rarefaction)",
       color = group_col, title = "Rarefaction curves (vegan::rarefy)") +
  theme_bw()
print(p_rare)
save_plot_multi(
  p_rare, file.path(EXPORT_DIR, "rarefaction_curves"),
  png_width = 8, png_height = 5, png_units = "in",
  svg_width = 18, svg_height = 10, svg_units = "cm"
)

# ============================================================
# 3) Taxonomic composition (aggregate + renormalize per sample)
# ============================================================
sample_cols_rel <- intersect(colnames(rel_abund_df), sample_ids)
if (length(sample_cols_rel) < 1) {
  stop("No sample columns in relative abundance matched the metadata IDs.")
}
tax_cols_rel <- get_tax_cols(rel_abund_df)
tax_level_resolved <- resolve_tax_level(rel_abund_df, tax_level)
if (is.na(tax_level_resolved)) {
  if (!feature_col %in% colnames(rel_abund_df)) rel_abund_df[[feature_col]] <- counts_raw[[feature_col]]
  tax_level_resolved <- feature_col
  message("No taxonomy columns found; using `", feature_col, "`.")
}

needed_cols <- unique(c(tax_cols_rel, sample_cols_rel))
rel_long <- rel_abund_df %>%
  dplyr::select(dplyr::all_of(needed_cols)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(sample_cols_rel), names_to = "SampleID", values_to = "Abundance") %>%
  dplyr::mutate(SampleID = as.character(SampleID)) %>%
  dplyr::left_join(meta, by = setNames(sample_id_col, "SampleID"))

# Handle missing taxa labels
rel_long[[tax_level_resolved]] <- ifelse(
  is.na(rel_long[[tax_level_resolved]]) | rel_long[[tax_level_resolved]] == "",
  "Unassigned",
  rel_long[[tax_level_resolved]]
)

# Aggregate at chosen level and renormalize per sample
rel_level <- rel_long %>%
  dplyr::group_by(SampleID, .data[[tax_level_resolved]]) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop_last") %>%
  dplyr::mutate(Abundance = Abundance / pmax(sum(Abundance, na.rm = TRUE), 1e-12)) %>%
  dplyr::ungroup()

# Keep top N taxa; collapse the rest as "Other"
top_taxa <- rel_level %>%
  dplyr::group_by(.data[[tax_level_resolved]]) %>%
  dplyr::summarise(MeanAbund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(MeanAbund)) %>%
  dplyr::slice_head(n = topN_taxa) %>%
  dplyr::pull(.data[[tax_level_resolved]])

rel_level$TaxonPlot <- ifelse(rel_level[[tax_level_resolved]] %in% top_taxa,
                              rel_level[[tax_level_resolved]], "Other")
fill_lab <- gsub("_", " ", tools::toTitleCase(tolower(tax_level_resolved)))

p_comp <- rel_level %>%
  dplyr::left_join(meta, by = "SampleID") %>%
  dplyr::mutate(SampleID = factor(SampleID, levels = unique(SampleID))) %>%
  ggplot(aes(x = SampleID, y = Abundance, fill = TaxonPlot)) +
  geom_col(width = 0.9) +
  labs(x = "Sample", y = "Relative abundance", fill = fill_lab,
       title = paste("Taxonomic composition by", fill_lab)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(p_comp)
save_plot_multi(
  p_comp, file.path(EXPORT_DIR, paste0("composition_", tolower(tax_level))),
  png_width = 10, png_height = 5, png_units = "in",
  svg_width = 22, svg_height = 12, svg_units = "cm"
)

# ============================================================
# 4) Alpha diversity (Chao1 on RAW; others on relative)
# ============================================================
counts_alpha <- counts_mat
meta_alpha   <- meta

if (filter_for_alpha) {
  # Filter samples by minimum library size
  lib_sizes_alpha <- rowSums(counts_alpha, na.rm = TRUE)
  keep_samples <- lib_sizes_alpha >= min_lib_alpha
  if (any(!keep_samples)) {
    message("Alpha: removing samples with library <", min_lib_alpha, ": ",
            paste(rownames(counts_alpha)[!keep_samples], collapse = ", "))
  }
  counts_alpha <- counts_alpha[keep_samples, , drop = FALSE]
  meta_alpha   <- meta_alpha[meta_alpha[[sample_id_col]] %in% rownames(counts_alpha), , drop = FALSE]
  # Filter features by prevalence and total counts
  prev <- colSums(counts_alpha > 0, na.rm = TRUE)
  tot  <- colSums(counts_alpha,     na.rm = TRUE)
  keep_features <- (prev >= min_prevalence) & (tot >= min_total_count)
  if (any(!keep_features)) {
    message("Alpha: removing ", sum(!keep_features), " features (low prevalence/abundance).")
  }
  counts_alpha <- counts_alpha[, keep_features, drop = FALSE]
}

if (rarefy_for_alpha) {
  if (nrow(counts_alpha) == 0) stop("Alpha: no samples after filtering.")
  lib_sizes_alpha <- rowSums(counts_alpha, na.rm = TRUE)
  if (any(lib_sizes_alpha <= 0)) stop("Alpha: zero-library samples after filtering.")
  if (is.na(rare_depth_alpha)) rare_depth_alpha <- min(lib_sizes_alpha)
  set.seed(123)
  counts_alpha <- vegan::rrarefy(counts_alpha, sample = rare_depth_alpha)
  message("Alpha: rarefied all samples to ", rare_depth_alpha, " reads.")
}

# Alpha metrics
S_obs   <- rowSums(counts_alpha > 0, na.rm = TRUE)
Chao1   <- apply(counts_alpha, 1, chao1_vec)                        # RAW
rel_mat_alpha <- counts_alpha / pmax(rowSums(counts_alpha), 1)      # per-sample proportions
Shannon <- vegan::diversity(rel_mat_alpha, index = "shannon")
Simpson <- vegan::diversity(rel_mat_alpha, index = "simpson")       # Gini–Simpson (1 - D)
Pielou  <- ifelse(S_obs > 1, Shannon / log(S_obs), NA_real_)

alpha_df <- data.frame(
  SampleID = rownames(counts_alpha),
  S_obs = S_obs, Chao1 = Chao1, Shannon = Shannon, Simpson = Simpson, Pielou = Pielou,
  check.names = FALSE
) %>% dplyr::left_join(meta_alpha, by = "SampleID")
alpha_df[[group_col]] <- as.factor(alpha_df[[group_col]])

# Stats + annotated boxplots
.alpha_stats_for_plot <- function(df, metric, group_col) {
  dd <- df[, c(metric, group_col)]; names(dd) <- c("value","group")
  dd <- dd %>% dplyr::filter(is.finite(value)) %>% tidyr::drop_na()
  # Normality by group (Shapiro)
  sh <- dd %>% dplyr::group_by(group) %>%
    dplyr::summarise(n=dplyr::n(), p = ifelse(n>=3, shapiro.test(value)$p.value, NA_real_), .groups="drop")
  normal_ok <- all(is.na(sh$p) | sh$p > 0.05)
  # Homoscedasticity (Levene)
  lev <- tryCatch(rstatix::levene_test(value ~ group, data = dd, center = "median"),
                  error = function(e) NULL)
  homosked_ok <- if (!is.null(lev)) lev$p > 0.05 else FALSE
  # Global + post-hoc
  if (normal_ok && homosked_ok) {
    method <- "ANOVA"
    global <- broom::tidy(aov(value ~ group, data = dd)) %>% dplyr::slice(1) %>% dplyr::transmute(p = p.value)
    post   <- rstatix::tukey_hsd(dd, value ~ group) %>% dplyr::mutate(test = "Tukey")
  } else {
    method <- "Kruskal"
    global <- rstatix::kruskal_test(dd, value ~ group) %>% dplyr::transmute(p = p)
    post   <- rstatix::dunn_test(dd, value ~ group, p.adjust.method = "holm") %>% dplyr::mutate(test = "Dunn")
  }
  list(method=method, global_p=global$p[1], posthoc=post, shapiro=sh, levene=lev)
}

.alpha_plot_annotated <- function(df, metric, ylab, title) {
  st <- .alpha_stats_for_plot(df, metric, group_col)
  y_max <- suppressWarnings(max(df[[metric]], na.rm = TRUE)); if (!is.finite(y_max)) y_max <- 1
  g <- ggplot(df, aes(x = .data[[group_col]], y = .data[[metric]], fill = .data[[group_col]])) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.7) +
    labs(x = group_col, y = ylab, title = paste0(title, " (", st$method, ")")) +
    theme_bw() + theme(legend.position = "none") +
    ggpubr::stat_compare_means(
      method = ifelse(st$method == "ANOVA", "anova", "kruskal.test"),
      label.y = y_max * 1.05
    )
  post <- st$posthoc
  if (!is.null(post) && nrow(post) > 0) {
    if (!"p.adj" %in% names(post)) post$p.adj <- post$p
    if (!"p.adj.signif" %in% names(post))
      post$p.adj.signif <- cut(post$p.adj, breaks=c(-Inf, .001, .01, .05, .1, Inf),
                               labels=c("***","**","*",".","ns"))
    sig <- post %>% dplyr::filter(!is.na(p.adj), p.adj <= 0.05) %>%
      dplyr::select(group1, group2, p.adj, p.adj.signif)
    if (nrow(sig) > 0) {
      sig$y.position <- seq(from = y_max * 1.12, by = y_max * 0.08, length.out = nrow(sig))
      g <- g + ggpubr::stat_pvalue_manual(
        sig, label = "p.adj.signif",
        xmin = "group1", xmax = "group2",
        y.position = "y.position", tip.length = 0.01, hide.ns = TRUE
      )
    }
  }
  g
}

p_alpha1 <- .alpha_plot_annotated(alpha_df, "Chao1",  "Alpha Index Value",  "Chao1 Richness Estimator")
p_alpha2 <- .alpha_plot_annotated(alpha_df, "Shannon","Alpha Index Value",  "Shannon Diversity Index")
p_alpha3 <- .alpha_plot_annotated(alpha_df, "Simpson","Alpha Index Value",  "Simpson's Diversity Index")
p_alpha4 <- .alpha_plot_annotated(alpha_df, "Pielou", "Alpha Index Value",  "Pielou's Evenness Index")
print(p_alpha1); print(p_alpha2); print(p_alpha3); print(p_alpha4)

# Save alpha diagnostics and plots
save_alpha_diagnostics <- function(df, metric, group_col, outdir = EXPORT_DIR) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  dd <- df[, c(metric, group_col)]; names(dd) <- c("value","group")
  if (length(unique(stats::na.omit(dd$group))) < 2) {
    overall <- data.frame(metric = metric, method = NA, p_global = NA, note = "< 2 groups")
    utils::write.csv(overall, file.path(outdir, paste0("alpha_overall_", metric, ".csv")), row.names = FALSE)
    return(invisible(NULL))
  }
  st <- .alpha_stats_for_plot(df, metric, group_col)
  overall <- data.frame(metric = metric, method = st$method, p_global = st$global_p)
  utils::write.csv(overall, file.path(outdir, paste0("alpha_overall_", metric, ".csv")), row.names = FALSE)
  if (!is.null(st$posthoc) && nrow(st$posthoc) > 0)
    utils::write.csv(st$posthoc, file.path(outdir, paste0("alpha_posthoc_", metric, ".csv")), row.names = FALSE)
  if (!is.null(st$shapiro))
    utils::write.csv(st$shapiro, file.path(outdir, paste0("alpha_", metric, "_shapiro_by_group.csv")), row.names = FALSE)
  if (!is.null(st$levene))
    utils::write.csv(st$levene, file.path(outdir, paste0("alpha_", metric, "_levene.csv")), row.names = FALSE)
  invisible(st)
}
invisible(lapply(c("Chao1","Shannon","Simpson","Pielou"),
                 function(m) save_alpha_diagnostics(alpha_df, m, group_col, outdir = EXPORT_DIR)))

# Save alpha plots (explicit sizes)
save_plot_multi(p_alpha1, file.path(EXPORT_DIR, "alpha_chao1"),
                png_width = 5, png_height = 6, png_units = "in",
                svg_width = 12, svg_height = 14, svg_units = "cm"
)
save_plot_multi(p_alpha2, file.path(EXPORT_DIR, "alpha_shannon"),
                png_width = 5, png_height = 6, png_units = "in",
                svg_width = 12, svg_height = 14, svg_units = "cm"
)
save_plot_multi(p_alpha3, file.path(EXPORT_DIR, "alpha_simpson"),
                png_width = 5, png_height = 6, png_units = "in",
                svg_width = 12, svg_height = 14, svg_units = "cm"
)
save_plot_multi(p_alpha4, file.path(EXPORT_DIR, "alpha_pielou"),
                png_width = 5, png_height = 6, png_units = "in",
                svg_width = 12, svg_height = 14, svg_units = "cm"
)

# ============================================================
# 5) Beta diversity (distances, PCoA, NMDS, PERMANOVA, betadisper)
# ============================================================

# ---- Config (edit as needed) ----
min_lib_beta <- 1                     # remove samples with library < this (1 removes only zeros)
beta_methods <- c("bray","jaccard","aitchison")
# Explicit plot sizes per method/plot (no global defaults)
beta_sizes <- list(
  bray = list(
    pcoa_png = c(18, 12, "cm"),  pcoa_svg = c(18, 12, "cm"),
    nmds_png = c(18, 12, "cm"),  nmds_svg = c(18, 12, "cm")
  ),
  jaccard = list(
    pcoa_png = c(18, 12, "cm"),  pcoa_svg = c(18, 12, "cm"),
    nmds_png = c(18, 12, "cm"),  nmds_svg = c(18, 12, "cm")
  ),
  aitchison = list(
    pcoa_png = c(18, 12, "cm"),  pcoa_svg = c(18, 12, "cm"),
    nmds_png = c(18, 12, "cm"),  nmds_svg = c(18, 12, "cm")
  )
)

# ---- Helpers ----
percent_from_eig <- function(eig) {
  p <- eig / sum(abs(eig))
  round(100 * p[1:2], 1)
}
.clr <- function(x, pseudo = 1e-6) {
  x <- sweep(x, 1, pmax(rowSums(x), 1e-12), "/")
  x[x <= 0] <- pseudo
  l <- log(x)
  sweep(l, 1, rowMeans(l), "-")
}
pairwise_permanova <- function(dist_obj, grouping, permutations = 999) {
  gr <- factor(grouping)
  lev <- levels(gr)
  res <- data.frame(group1=character(), group2=character(), F=numeric(), R2=numeric(), p=numeric())
  if (length(lev) < 2) return(res)
  for (i in 1:(length(lev)-1)) {
    for (j in (i+1):length(lev)) {
      idx <- gr %in% c(lev[i], lev[j])
      sub_meta <- data.frame(g = droplevels(gr[idx]))
      ad <- vegan::adonis2(as.dist(as.matrix(dist_obj)[idx, idx]) ~ g,
                           data = sub_meta, permutations = permutations, by = "terms")
      res <- rbind(res,
                   data.frame(group1=lev[i], group2=lev[j],
                              F=ad$F[1], R2=ad$R2[1], p=ad$`Pr(>F)`[1]))
    }
  }
  res$p_adj <- p.adjust(res$p, method = "holm")
  res[order(res$p_adj), ]
}

# Diagnostic saver that prints paths, tries, and verifies files
save_plot_checked <- function(plot, file_png, file_svg,
                              png_width, png_height, png_units = "in", png_dpi = 300,
                              svg_width, svg_height, svg_units = "cm") {
  dir.create(dirname(file_png), showWarnings = FALSE, recursive = TRUE)
  
  message("[Saving PNG] ", normalizePath(file_png, mustWork = FALSE))
  ok_png <- tryCatch({
    ggplot2::ggsave(file_png, plot,
                    width = as.numeric(png_width),
                    height = as.numeric(png_height),
                    units = png_units, dpi = png_dpi)
    TRUE
  }, error = function(e) { message("PNG save error: ", e$message); FALSE })
  if (ok_png && !file.exists(file_png)) message("WARN: PNG not found after save: ", file_png)
  
  message("[Saving SVG] ", normalizePath(file_svg, mustWork = FALSE))
  ok_svg <- tryCatch({
    ggplot2::ggsave(file_svg, plot,
                    device = svglite::svglite,
                    width = as.numeric(svg_width),
                    height = as.numeric(svg_height),
                    units = svg_units)
    TRUE
  }, error = function(e) { message("SVG save error: ", e$message); FALSE })
  if (ok_svg && !file.exists(file_svg)) message("WARN: SVG not found after save: ", file_svg)
  
  invisible(ok_png && ok_svg)
}

# Quick write-permission test 
cat("write-test\n", file = file.path(EXPORT_DIR, ".__write_test__.txt"))

# ---- Prepare base matrices for beta ----
lib_beta <- rowSums(counts_mat, na.rm = TRUE)
keep_beta <- lib_beta >= min_lib_beta
counts_mat_beta <- counts_mat[keep_beta, , drop = FALSE]
meta_beta <- meta[meta[[sample_id_col]] %in% rownames(counts_mat_beta), , drop = FALSE]
stopifnot(all(rownames(counts_mat_beta) == meta_beta[[sample_id_col]]))
if (!group_col %in% names(meta_beta)) {
  meta_beta[[group_col]] <- "All"
  message("`", group_col, "` not found in metadata; set to 'All' for beta.")
}

# ---- Core pipeline (no saving inside; returns plots/frames) ----
run_beta_pipeline <- function(metric = c("bray","jaccard","aitchison"), show_plots = TRUE) {
  metric <- match.arg(metric)
  
  Xraw <- counts_mat_beta
  grp  <- factor(meta_beta[[group_col]])
  n_groups <- nlevels(grp)
  if (nrow(Xraw) < 2) stop("Beta: need >= 2 samples after filtering.")
  
  # Distance
  if (metric == "bray") {
    X <- Xraw / pmax(rowSums(Xraw), 1)                # relative
    dist_obj <- vegan::vegdist(X, method = "bray")
  } else if (metric == "jaccard") {
    X <- (Xraw > 0) * 1                               # presence/absence
    dist_obj <- vegan::vegdist(X, method = "jaccard", binary = TRUE)
  } else if (metric == "aitchison") {
    X <- .clr(Xraw)                                   # CLR
    dist_obj <- stats::dist(X, method = "euclidean")
  }
  if (anyNA(as.matrix(dist_obj))) stop("Beta: distance matrix contains NA. Check empty samples/IDs.")
  
  # PCoA (cmdscale) with safe fallback (ape::pcoa + Cailliez)
  pc <- tryCatch(cmdscale(dist_obj, eig = TRUE, k = 2), error = function(e) NULL)
  if (is.null(pc) || anyNA(pc$points)) {
    pc2 <- ape::pcoa(dist_obj, correction = "cailliez")
    pvar <- round(100 * pc2$values$Relative_eig[1:2], 1)
    pcoa_df <- data.frame(SampleID = rownames(Xraw),
                          Axis1 = pc2$vectors[,1], Axis2 = pc2$vectors[,2]) %>%
      dplyr::left_join(meta_beta, by = "SampleID")
  } else {
    pvar <- percent_from_eig(pc$eig)
    pcoa_df <- data.frame(SampleID = rownames(Xraw),
                          Axis1 = pc$points[,1], Axis2 = pc$points[,2]) %>%
      dplyr::left_join(meta_beta, by = "SampleID")
  }
  
  # NMDS
  nmds <- tryCatch(
    suppressWarnings(vegan::metaMDS(dist_obj, k = 2, trymax = 100,
                                    autotransform = FALSE, trace = FALSE)),
    error = function(e) { message("NMDS(", metric, ") failed: ", e$message); NULL }
  )
  nmds_df <- if (!is.null(nmds)) {
    data.frame(SampleID = rownames(Xraw),
               NMDS1 = nmds$points[,1], NMDS2 = nmds$points[,2]) %>%
      dplyr::left_join(meta_beta, by = "SampleID")
  } else NULL
  
  # Stats init
  bd <- bd_perm <- NULL
  ad_glob <- NULL
  pw <- NULL
  if (n_groups >= 2) {
    bd <- vegan::betadisper(dist_obj, group = grp)
    bd_perm <- vegan::permutest(bd, permutations = 999)
    ad_glob <- vegan::adonis2(dist_obj ~ grp, permutations = 999)
    pw <- pairwise_permanova(dist_obj, grp, permutations = 999)
    sub_perma <- paste0("PERMANOVA: F=", round(ad_glob$F[1],3),
                        " | R²=", round(ad_glob$R2[1],3),
                        " | p=", formatC(ad_glob$`Pr(>F)`[1], format="g", digits=3))
  } else {
    sub_perma <- "PERMANOVA: n/a (< 2 groups)"
  }
  
  # Plots
  can_ellipse <- all(table(grp) >= 3)
  p_pcoa <- ggplot2::ggplot(pcoa_df, ggplot2::aes(x = Axis1, y = Axis2, color = .data[[group_col]])) +
    ggplot2::geom_point(size = 3, alpha = 0.9, na.rm = TRUE) +
    ggrepel::geom_text_repel(ggplot2::aes(label = SampleID), size = 3, max.overlaps = 20, show.legend = FALSE, na.rm = TRUE) +
    ggplot2::labs(x = paste0("PCoA1 (", pvar[1], "%)"), y = paste0("PCoA2 (", pvar[2], "%)"),
                  color = group_col,
                  title = paste0("PCoA - ", tools::toTitleCase(metric)),
                  subtitle = sub_perma) +
    ggplot2::theme_bw()
  if (can_ellipse) {
    p_pcoa <- p_pcoa + ggplot2::stat_ellipse(ggplot2::aes(group = .data[[group_col]]),
                                             linewidth = 0.6, alpha = 0.6, show.legend = FALSE)
  }
  
  if (!is.null(nmds_df)) {
    p_nmds <- ggplot2::ggplot(nmds_df, ggplot2::aes(x = NMDS1, y = NMDS2, color = .data[[group_col]])) +
      ggplot2::geom_point(size = 3, alpha = 0.9, na.rm = TRUE) +
      ggrepel::geom_text_repel(ggplot2::aes(label = SampleID), size = 3, max.overlaps = 20, show.legend = FALSE, na.rm = TRUE) +
      ggplot2::labs(x = "NMDS1", y = "NMDS2", color = group_col,
                    title = paste0("NMDS - ", tools::toTitleCase(metric)),
                    subtitle = paste0("stress = ", round(nmds$stress, 3), " | ", sub_perma)) +
      ggplot2::theme_bw()
    if (can_ellipse) {
      p_nmds <- p_nmds + ggplot2::stat_ellipse(ggplot2::aes(group = .data[[group_col]]),
                                               linewidth = 0.6, alpha = 0.6, show.legend = FALSE)
    }
  } else {
    p_nmds <- NULL
  }
  
  # Preview
  if (isTRUE(show_plots)) {
    message("[Preview] PCoA - ", metric); print(p_pcoa)
    if (!is.null(p_nmds)) { message("[Preview] NMDS - ", metric); print(p_nmds) }
  }
  
  # Export tables 
  utils::write.csv(as.matrix(dist_obj), file.path(EXPORT_DIR, paste0("distance_", metric, ".csv")))
  data.table::fwrite(pcoa_df, file.path(EXPORT_DIR, paste0("pcoa_scores_", metric, ".csv")))
  if (!is.null(nmds_df)) data.table::fwrite(nmds_df, file.path(EXPORT_DIR, paste0("nmds_scores_", metric, ".csv")))
  
  # betadisper outputs (robust to version differences)
  if (!is.null(bd_perm) && !is.null(bd_perm$tab)) {
    bd_perm_tab <- as.data.frame(bd_perm$tab)
    nperm <- tryCatch({
      if (!is.null(bd_perm$permutations)) {
        if (is.numeric(bd_perm$permutations)) as.integer(bd_perm$permutations)
        else as.integer(attr(bd_perm$permutations, "nperm"))
      } else if (!is.null(attr(bd_perm, "permutations"))) {
        as.integer(attr(bd_perm, "permutations"))
      } else if (!is.null(bd_perm$perm)) {
        np <- attr(bd_perm$perm, "nperm")
        if (is.null(np)) as.integer(length(bd_perm$perm)) else as.integer(np)
      } else NA_integer_
    }, error = function(e) NA_integer_)
    if (length(nperm) == 1 && is.finite(nperm)) {
      bd_perm_tab$N.Perm <- rep(nperm, nrow(bd_perm_tab))
    }
    utils::write.csv(bd_perm_tab, file.path(EXPORT_DIR, paste0("betadisper_global_", metric, ".csv")), row.names = TRUE)
    
    tuk <- tryCatch(TukeyHSD(bd), error = function(e) NULL)
    if (!is.null(tuk)) {
      for (nm in names(tuk)) {
        tk <- as.data.frame(tuk[[nm]]); tk$contrast <- rownames(tk)
        utils::write.csv(tk, file.path(EXPORT_DIR, paste0("betadisper_tukey_", metric, "_", nm, ".csv")), row.names = FALSE)
      }
    }
  } else if (!is.null(bd_perm)) {
    message("betadisper: permutest returned no table; skipping CSV.")
  }
  
  ad_glob_df <- tryCatch({
    df <- as.data.frame(ad_glob); df$Term <- rownames(df); rownames(df) <- NULL; df
  }, error = function(e) data.frame(Term="grp", Note="permanova skipped (< 2 groups)"))
  utils::write.csv(ad_glob_df, file.path(EXPORT_DIR, paste0("permanova_global_", metric, ".csv")), row.names = FALSE)
  if (!is.null(pw) && nrow(pw) > 0) {
    utils::write.csv(pw, file.path(EXPORT_DIR, paste0("permanova_pairwise_", metric, ".csv")), row.names = FALSE)
  }
  
  # Return plots and frames (saving happens outside)
  list(metric = metric, p_pcoa = p_pcoa, p_nmds = p_nmds,
       pcoa_df = pcoa_df, nmds_df = nmds_df)
}

# ---- Run all requested methods (robust to individual failures) ----
beta_results <- lapply(beta_methods, function(m)
  tryCatch(run_beta_pipeline(m, show_plots = TRUE),
           error = function(e){ message("Method ", m, " failed: ", e$message); NULL })
)
names(beta_results) <- beta_methods
beta_results <- Filter(Negate(is.null), beta_results)  # drop failed methods (if any)
message("beta_results available: ", paste(names(beta_results), collapse = ", "))

# ---- Save figures with explicit sizes per method (diagnostic saver) ----
for (m in names(beta_results)) {
  r <- beta_results[[m]]
  s <- beta_sizes[[m]]
  if (is.null(s)) { message("No sizes declared for method: ", m); next }
  
  # Coerce sizes safely
  pcoa_png_w <- as.numeric(s$pcoa_png[[1]]); pcoa_png_h <- as.numeric(s$pcoa_png[[2]]); pcoa_png_u <- as.character(s$pcoa_png[[3]])
  pcoa_svg_w <- as.numeric(s$pcoa_svg[[1]]); pcoa_svg_h <- as.numeric(s$pcoa_svg[[2]]); pcoa_svg_u <- as.character(s$pcoa_svg[[3]])
  
  # PCoA
  save_plot_checked(
    r$p_pcoa,
    file_png = file.path(EXPORT_DIR, paste0("pcoa_", m, ".png")),
    file_svg = file.path(EXPORT_DIR, paste0("pcoa_", m, ".svg")),
    png_width = pcoa_png_w, png_height = pcoa_png_h, png_units = pcoa_png_u,
    svg_width = pcoa_svg_w, svg_height = pcoa_svg_h, svg_units = pcoa_svg_u
  )
  
  # NMDS 
  if (!is.null(r$p_nmds)) {
    nmds_png_w <- as.numeric(s$nmds_png[[1]]); nmds_png_h <- as.numeric(s$nmds_png[[2]]); nmds_png_u <- as.character(s$nmds_png[[3]])
    nmds_svg_w <- as.numeric(s$nmds_svg[[1]]); nmds_svg_h <- as.numeric(s$nmds_svg[[2]]); nmds_svg_u <- as.character(s$nmds_svg[[3]])
    
    save_plot_checked(
      r$p_nmds,
      file_png = file.path(EXPORT_DIR, paste0("nmds_", m, ".png")),
      file_svg = file.path(EXPORT_DIR, paste0("nmds_", m, ".svg")),
      png_width = nmds_png_w, png_height = nmds_png_h, png_units = nmds_png_u,
      svg_width = nmds_svg_w, svg_height = nmds_svg_h, svg_units = nmds_svg_u
    )
  }
}

# ---- List saved files ----
print(list.files(EXPORT_DIR, pattern = "pcoa|nmds|distance|permanova|betadisper", full.names = TRUE))
message("Beta analysis done. Outputs in: ", normalizePath(EXPORT_DIR, mustWork = FALSE))

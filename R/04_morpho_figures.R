# ==============================================================================
# R/04_morpho_figures.R
# Step 4: Publication-Grade Figures (Okabe-Ito Palette, Italicized Binomials,
# Non-Overlapping Labels, Embedded Statistical Summaries, Vector PDF + 600 DPI PNG)
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(ggtree)
  library(ape)
  library(dendextend)
})
reset_project_seed()

message("=== [04] Generating Publication Figures ===")

consensus_qc  <- read.csv(file.path(DIRS$processed, "morpho_consensus_qc.csv"), stringsAsFactors = FALSE)
pca_objs      <- readRDS(file.path(DIRS$rds, "morpho_pca_objects.rds"))
clust_objs    <- readRDS(file.path(DIRS$rds, "morpho_clustering_objects.rds"))

cont_cols          <- pca_objs$cont_cols
pretty_trait_names <- pca_objs$pretty_trait_names

# Prepare cultivar metadata annotation columns
consensus_qc <- consensus_qc %>%
  mutate(
    Short_Label = paste0("Musa '", Display_Name, "' (", sub(" \\(Local\\)", "", Genome_Group), ")")
  )

# ==============================================================================
# FIGURE 1: Trait Collinearity (Spearman), PCA Biplot & Z-Score Heatmap
# ==============================================================================

# --- Panel A: Spearman Correlation Matrix -------------------------------------
X_qc <- pca_objs$X_qc_scaled
cor_mat_s <- cor(X_qc, method = "spearman")
colnames(cor_mat_s) <- unname(pretty_trait_names[colnames(cor_mat_s)])
rownames(cor_mat_s) <- unname(pretty_trait_names[rownames(cor_mat_s)])

cor_long <- as.data.frame(as.table(cor_mat_s), stringsAsFactors = FALSE)
names(cor_long) <- c("Trait_1", "Trait_2", "Spearman_rho")
cor_long$Trait_1 <- factor(cor_long$Trait_1, levels = unname(pretty_trait_names))
cor_long$Trait_2 <- factor(cor_long$Trait_2, levels = rev(unname(pretty_trait_names)))

p1_corr <- ggplot(cor_long, aes(x = Trait_1, y = Trait_2, fill = Spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Spearman_rho),
                color = abs(Spearman_rho) > 0.65),
            size = 3.0, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "#111111", "TRUE" = "white")) +
  scale_fill_viridis_c(option = "cividis", limits = c(-1, 1), name = "Spearman \u03c1") +
  labs(x = NULL, y = NULL) +
  theme_thesis(base_size = 9.5) +
  theme(
    axis.text.x     = element_text(angle = 28, hjust = 1, vjust = 1, size = 8.5),
    axis.text.y     = element_text(size = 8.5),
    legend.position = "right"
  )

# --- Panel B: PCA Biplot (PC1 vs PC2) on Standardized Morphometrics -----------
pca_mod <- pca_objs$pca_qc_full
eig_pct <- pca_mod$eig[, 2]

ind_df <- as.data.frame(pca_mod$ind$coord[, 1:2])
colnames(ind_df) <- c("PC1", "PC2")
ind_df$Cultivar_Key <- rownames(ind_df)
ind_df <- ind_df %>% left_join(consensus_qc, by = "Cultivar_Key")

var_df <- as.data.frame(pca_mod$var$coord[, 1:2])
colnames(var_df) <- c("PC1", "PC2")
var_df$Trait <- rownames(var_df)
var_df$Trait_Label <- unname(pretty_trait_names[var_df$Trait])

# Scale loading arrows to fit cultivar coordinate space cleanly
arrow_scale <- 2.35
var_df <- var_df %>%
  mutate(
    xend = PC1 * arrow_scale,
    yend = PC2 * arrow_scale
  )

# Unified repulsion dataframe so cultivar labels and trait labels mutually repel
# Labels are black; only the cultivar points are colored by Genome Group
repel_all_df <- bind_rows(
  ind_df %>%
    transmute(
      x          = PC1,
      y          = PC2,
      label      = Short_Label,
      Label_Type = "Cultivar",
      Color_Hex  = "#111111",
      Font_Face  = "bold.italic",
      Fill_Col   = alpha("white", 0.85)
    ),
  var_df %>%
    transmute(
      x          = xend,
      y          = yend,
      label      = Trait_Label,
      Label_Type = "Trait",
      Color_Hex  = "#444444",
      Font_Face  = "plain",
      Fill_Col   = alpha("#F4F6F8", 0.92)
    )
)

p1_pca <- ggplot(ind_df, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  # Loading vectors
  geom_segment(
    data = var_df,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
    color = "#555555", linewidth = 0.55, alpha = 0.75,
    inherit.aes = FALSE
  ) +
  # Cultivar points colored by Genome Group
  geom_point(aes(color = Genome_Group), shape = 19, size = 4.2, alpha = 0.95) +
  # Single unified repulsion layer so NO cultivar label and trait label ever overlap
  ggrepel::geom_label_repel(
    data = repel_all_df,
    aes(x = x, y = y, label = label),
    color        = repel_all_df$Color_Hex,
    fill         = repel_all_df$Fill_Col,
    fontface     = repel_all_df$Font_Face,
    size         = ifelse(repel_all_df$Label_Type == "Cultivar", 2.95, 2.50),
    linewidth    = ifelse(repel_all_df$Label_Type == "Cultivar", 0.15, 0.25),
    label.padding= unit(0.18, "lines"),
    box.padding  = 0.45,
    point.padding= 0.35,
    force        = 5.0,
    max.overlaps = 50,
    segment.color= "#666666",
    segment.size = 0.35,
    seed         = PROJECT_SEED,
    inherit.aes  = FALSE
  ) +
  scale_color_manual(values = GENOME_PALETTE, name = "*Musa* Genome Group") +
  scale_x_continuous(expand = expansion(mult = 0.22)) +
  scale_y_continuous(expand = expansion(mult = 0.22)) +
  labs(
    x = sprintf("PC1 (%.1f%% variance explained)", eig_pct[1]),
    y = sprintf("PC2 (%.1f%% variance explained)", eig_pct[2])
  ) +
  theme_thesis(base_size = 9.5) +
  theme(legend.position = "right")

# --- Panel C: Standardized Phenotypic Z-Score Profile across Cultivars --------
z_long <- as.data.frame(X_qc) %>%
  mutate(Cultivar_Key = rownames(X_qc)) %>%
  pivot_longer(cols = all_of(cont_cols), names_to = "Trait", values_to = "Z_Score") %>%
  left_join(consensus_qc %>% select(Cultivar_Key, Display_Name, Genome_Group), by = "Cultivar_Key") %>%
  mutate(
    Trait_Label    = factor(unname(pretty_trait_names[Trait]), levels = unname(pretty_trait_names)),
    Cultivar_Label = paste0("*Musa* '", Display_Name, "' (", sub(" \\(Local\\)", "", Genome_Group), ")")
  )

# Order cultivars by Track A dendrogram tip order
hc_order_keys <- clust_objs$pv_trackA$hclust$labels[clust_objs$pv_trackA$hclust$order]
label_levels  <- unique(z_long$Cultivar_Label[match(hc_order_keys, z_long$Cultivar_Key)])
z_long$Cultivar_Label <- factor(z_long$Cultivar_Label, levels = label_levels)

p1_zheat <- ggplot(z_long, aes(x = Trait_Label, y = Cultivar_Label, fill = Z_Score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.2f", Z_Score),
                color = abs(Z_Score) > 1.3),
            size = 2.85, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "#111111", "TRUE" = "white")) +
  scale_fill_viridis_c(option = "viridis", limits = c(-2.3, 2.6), name = "Trait Z-Score\n(\u03bc=0, \u03c3=1)") +
  labs(x = NULL, y = NULL) +
  theme_thesis(base_size = 9.5) +
  theme(
    axis.text.y     = ggtext::element_markdown(size = 8.8, color = "#111111"),
    axis.text.x     = element_text(size = 8.8, face = "bold"),
    legend.position = "right"
  )

fig1_combined <- ((p1_corr | p1_pca) + plot_layout(widths = c(1, 1.35))) / p1_zheat +
  plot_layout(heights = c(1.22, 0.95)) +
  plot_annotation(tag_levels = "A")

save_thesis_plot(fig1_combined, "Fig1_Morpho_QC_PCA_Biplot", width = 14.2, height = 9.8, dpi = 600)

# ==============================================================================
# FIGURE 2: Cophenetic Correlation Benchmark & Shephard Goodness-of-Fit
# ==============================================================================

bench_qc <- clust_objs$bench_df %>%
  filter(Dataset == "QC_Corrected") %>%
  mutate(
    Model_Label = paste0(Distance_Metric, "  x  ", Linkage_Method),
    Is_Winner   = (Distance_Metric == clust_objs$best_trackA$Distance_Metric &
                     Linkage_Code == clust_objs$best_trackA$Linkage_Code) |
                  (Distance_Metric == clust_objs$best_trackB$Distance_Metric &
                     Linkage_Code == clust_objs$best_trackB$Linkage_Code)
  ) %>%
  arrange(r_coph_Pearson) %>%
  mutate(Model_Label = factor(Model_Label, levels = unique(Model_Label)))

p2_forest <- ggplot(bench_qc, aes(x = r_coph_Pearson, y = Model_Label,
                                  color = Linkage_Method, shape = Track)) +
  geom_vline(xintercept = 0.80, linetype = "dashed", color = "#D55E00", linewidth = 0.7) +
  geom_errorbar(aes(xmin = CI95_Low, xmax = CI95_High), width = 0.28, linewidth = 0.55, alpha = 0.8, orientation = "y") +
  geom_point(aes(size = Is_Winner), alpha = 0.95) +
  geom_text(aes(label = sprintf("%.3f", r_coph_Pearson)),
            nudge_y = 0.36, size = 2.55, fontface = "bold", show.legend = FALSE) +
  annotate("text", x = 0.81, y = 4.5, label = "Target Threshold\n(r_coph >= 0.80)",
           hjust = 0, size = 2.85, fontface = "bold", color = "#D55E00", lineheight = 1.05) +
  scale_size_manual(values = c("FALSE" = 2.8, "TRUE" = 4.5), guide = "none") +
  scale_color_manual(
    values = c(
      "UPGMA (Average)"  = "#0072B2",
      "Complete Linkage" = "#009E73",
      "Single Linkage"   = "#E69F00",
      "Ward.D2"          = "#CC79A7"
    ),
    name = "Linkage Algorithm"
  ) +
  scale_x_continuous(breaks = seq(0.3, 1.0, by = 0.1)) +
  coord_cartesian(xlim = c(0.28, 1.02)) +
  labs(
    x = "Cophenetic Correlation Coefficient (r_coph)",
    y = "Distance Metric  x  Linkage Method"
  ) +
  theme_thesis(base_size = 9.0) +
  theme(legend.position = "bottom", legend.box = "vertical")

# --- Panel B: Shephard Cophenetic Calibration Scatterplot ---------------------
build_shephard_df <- function(d_raw, hc_obj, model_tag) {
  d_coph <- cophenetic(hc_obj)
  tibble(
    Model      = model_tag,
    Raw_Dist   = as.numeric(d_raw),
    Coph_Dist  = as.numeric(d_coph)
  )
}

shep_df <- bind_rows(
  build_shephard_df(
    clust_objs$D_A_N10,
    clust_objs$pv_trackA$hclust,
    sprintf("Track A: Euclidean_VIF + UPGMA\nr_coph = %.3f (95%% CI [%.3f, %.3f])",
            clust_objs$best_trackA$r_coph_Pearson,
            clust_objs$best_trackA$CI95_Low,
            clust_objs$best_trackA$CI95_High)
  ),
  build_shephard_df(
    clust_objs$D_B_N10,
    clust_objs$pv_trackB$hclust,
    sprintf("Track B: Gower Unweighted + UPGMA\nr_coph = %.3f (95%% CI [%.3f, %.3f])",
            clust_objs$best_trackB$r_coph_Pearson,
            clust_objs$best_trackB$CI95_Low,
            clust_objs$best_trackB$CI95_High)
  )
)

p2_shephard <- ggplot(shep_df, aes(x = Raw_Dist, y = Coph_Dist)) +
  geom_point(color = "#0072B2", size = 2.4, alpha = 0.75) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#D55E00", fill = "#F5C6A5", linewidth = 0.8) +
  facet_wrap(~ Model, scales = "free", ncol = 1) +
  labs(
    x = "Observed Pairwise Trait Distance (D_morph)",
    y = "Ultrametric Dendrogram Cophenetic Distance (D_coph)"
  ) +
  theme_thesis(base_size = 9.0)

fig2_combined <- (p2_forest | p2_shephard) +
  plot_layout(widths = c(1.35, 1.0)) +
  plot_annotation(tag_levels = "A")

save_thesis_plot(fig2_combined, "Fig2_Cophenetic_Benchmark", width = 13.8, height = 7.2, dpi = 600)

# ==============================================================================
# FIGURE 3: Multiscale Bootstrapped Dendrograms (Track A & Track B, nboot=10,000)
# ==============================================================================

build_ggtree_pvclust_panel <- function(pv_obj, x_axis_label) {
  hc <- pv_obj$hclust
  tr <- ape::as.phylo(hc)

  n_tips  <- length(tr$tip.label)
  n_nodes <- tr$Nnode

  get_hc_clade <- function(idx) {
    if (idx < 0) return(hc$labels[-idx])
    c(Recall(hc$merge[idx, 1]), Recall(hc$merge[idx, 2]))
  }
  hc_clades <- sapply(seq_len(nrow(hc$merge)), function(i) {
    paste(sort(c(get_hc_clade(hc$merge[i, 1]), get_hc_clade(hc$merge[i, 2]))), collapse = "|")
  })

  ape_clades <- sapply((n_tips + 1):(n_tips + n_nodes), function(nd) {
    sub_tr <- ape::extract.clade(tr, nd)
    paste(sort(sub_tr$tip.label), collapse = "|")
  })

  match_idx <- match(ape_clades, hc_clades)
  au_vals   <- round(pv_obj$edges[match_idx, "au"] * 100, 0)
  bp_vals   <- round(pv_obj$edges[match_idx, "bp"] * 100, 0)

  node_df <- tibble(
    node        = (n_tips + 1):(n_tips + n_nodes),
    AU          = au_vals,
    BP          = bp_vals,
    Support_Lbl = sprintf("%d / %d", au_vals, bp_vals),
    Sig_Class   = ifelse(au_vals >= 95, "Significant (AU \u2265 95%)", "Moderate / Weak (AU < 95%)")
  )

  tip_df <- consensus_qc %>%
    select(label = Cultivar_Key, Display_Name, Genome_Group) %>%
    mutate(
      Tip_Display = paste0("*Musa* '", Display_Name, "' (", sub(" \\(Local\\)", "", Genome_Group), ")")
    )

  max_x <- max(ape::node.depth.edgelength(tr))

  p <- ggtree(tr, linewidth = 0.75, color = "#222222") %<+% tip_df %<+% node_df +
    geom_nodepoint(aes(fill = Sig_Class), shape = 21, size = 3.0, color = "#111111", na.rm = TRUE) +
    geom_label(
      aes(label = Support_Lbl, fill = Sig_Class),
      size = 2.55, fontface = "bold", color = "black",
      label.padding = unit(0.14, "lines"), linewidth = 0.2,
      hjust = 1.12, vjust = -0.35, na.rm = TRUE, show.legend = FALSE
    ) +
    geom_tippoint(aes(color = Genome_Group), shape = 19, size = 3.8, na.rm = TRUE) +
    ggtext::geom_richtext(
      aes(label = Tip_Display),
      color = "#111111",
      hjust = 0, nudge_x = max_x * 0.03,
      size = 3.15, fill = NA, label.color = NA, fontface = "bold",
      na.rm = TRUE, show.legend = FALSE
    ) +
    scale_color_manual(values = GENOME_PALETTE, name = "*Musa* Genome Group") +
    scale_fill_manual(
      values = c(
        "Significant (AU \u2265 95%)"     = "#A8E6CF",
        "Moderate / Weak (AU < 95%)" = "#FFD3B6"
      ),
      name = "Multiscale Bootstrap (AU / BP %)"
    ) +
    coord_cartesian(clip = "off") +
    hexpand(0.58, direction = 1) +
    hexpand(0.08, direction = -1) +
    theme_tree2() +
    labs(x = x_axis_label) +
    theme(
      axis.title.x    = element_text(size = 9.5, face = "bold"),
      legend.position = "bottom",
      legend.box      = "horizontal",
      legend.title    = ggtext::element_markdown(size = 8.8, face = "bold"),
      legend.text     = ggtext::element_markdown(size = 8.2),
      plot.margin     = margin(10, 20, 10, 10)
    )

  p
}

p3_trackA <- build_ggtree_pvclust_panel(
  clust_objs$pv_trackA,
  x_axis_label = "Euclidean Branch Height (UPGMA)"
)

p3_trackB <- build_ggtree_pvclust_panel(
  clust_objs$pv_trackB,
  x_axis_label = "Gower Dissimilarity Branch Height (UPGMA)"
)

fig3_combined <- (p3_trackA | p3_trackB) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(legend.position = "bottom")
  )

save_thesis_plot(fig3_combined, "Fig3_Morpho_Dendrogram_Bootstrapped", width = 14.0, height = 7.6, dpi = 600)

message("[04] All Phase 1 publication figures generated successfully.")


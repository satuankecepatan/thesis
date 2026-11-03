# ==============================================================================
# R/11_concordance_figures.R
# Phase 3 - Step 3: Publication-Grade Cophylogenetic & Concordance Figures
# (Title/Subtitle-Free, Black Italicized Binomial Labels, Locked Okabe-Ito
# Genome Group Point Palette, Vector PDF + 600 DPI PNG)
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(dendextend)
  library(ape)
  library(ggrepel)
})
reset_project_seed()

message("=== [11] Generating Phase 3 Comparative Concordance Publication Figures ===")

topo_objs <- readRDS(file.path(DIRS$rds, "concordance_topological_objects.rds"))
conc_objs <- readRDS(file.path(DIRS$rds, "concordance_analysis_objects.rds"))
meta_df   <- read.csv(file.path(DIRS$processed, "cultivar_metadata.csv"), stringsAsFactors = FALSE)

shared_meta <- topo_objs$shared_meta %>%
  mutate(
    Clean_Label_Markdown = paste0("*Musa* '", Display_Name, "' (", sub(" \\(Local\\)", "", Genome_Group), ")"),
    Clean_Label_Plain    = paste0("Musa '", Display_Name, "' (", sub(" \\(Local\\)", "", Genome_Group), ")"),
    Genome_Group         = factor(Genome_Group, levels = names(GENOME_PALETTE))
  )

# ==============================================================================
# FIGURE 7: Cophylogenetic Tanglegram (Panel A) & Baker's Gamma Null (Panel B)
# ==============================================================================

# Extract step2side-untangled dendrograms for the primary pair (Track A vs Phred-Weighted)
dl_prim   <- topo_objs$untangled_pairs[["TrackA_vs_Phred"]]$dend_list
dend_left <- dl_prim[[1]] # Morphological Track A
dend_right <- dl_prim[[2]] # Molecular Phred-Weighted ML

# Helper to extract segment and leaf coordinates from a dendrogram object
# with branch heights strictly normalized to [0, 1]
extract_dend_geom <- function(dend) {
  ggd  <- dendextend::as.ggdend(dend)
  segs <- as_tibble(ggd$segments)
  max_y <- max(c(segs$y, segs$yend), na.rm = TRUE)
  if (is.finite(max_y) && max_y > 0) {
    segs$y    <- pmax(0, segs$y / max_y)
    segs$yend <- pmax(0, segs$yend / max_y)
  }
  lbls <- as_tibble(ggd$labels)
  list(segments = segs, labels = lbls)
}

g_left  <- extract_dend_geom(dend_left)
g_right <- extract_dend_geom(dend_right)

# Map horizontal coordinates with generous label & spline corridors:
# - Left tree (Morphology): root at x = -1.90, tips at x = -1.18 (width = 0.72)
# - Left tip labels: x = -1.14 to -0.32 (hjust = 0)
# - Center connecting splines: x = -0.26 to +0.26
# - Right tip labels: x = +0.32 to +1.14 (hjust = 1)
# - Right tree (Molecular ITS): tips at x = +1.18, root at x = +1.90 (width = 0.72)
map_left_x  <- function(h) -1.18 - 0.72 * h
map_right_x <- function(h)  1.18 + 0.72 * h

segs_left <- tibble(
  x    = map_left_x(g_left$segments$y),
  xend = map_left_x(g_left$segments$yend),
  y    = g_left$segments$x,
  yend = g_left$segments$xend
)

segs_right <- tibble(
  x    = map_right_x(g_right$segments$y),
  xend = map_right_x(g_right$segments$yend),
  y    = g_right$segments$x,
  yend = g_right$segments$xend
)

tips_left <- tibble(
  Cultivar_Key = g_left$labels$label,
  x_tip        = map_left_x(0),
  y_tip        = g_left$labels$x
) %>%
  left_join(shared_meta, by = "Cultivar_Key")

tips_right <- tibble(
  Cultivar_Key = g_right$labels$label,
  x_tip        = map_right_x(0),
  y_tip        = g_right$labels$x
) %>%
  left_join(shared_meta, by = "Cultivar_Key")

links_df <- tips_left %>%
  select(Cultivar_Key, Genome_Group, y_left = y_tip) %>%
  inner_join(tips_right %>% select(Cultivar_Key, y_right = y_tip), by = "Cultivar_Key")

sigmoid_curves <- bind_rows(lapply(seq_len(nrow(links_df)), function(r) {
  row <- links_df[r, ]
  t_seq <- seq(0, 1, length.out = 60)
  x_seq <- -0.26 + 0.52 * t_seq
  s_val <- 0.5 * (1 - cos(pi * t_seq))
  y_seq <- row$y_left + (row$y_right - row$y_left) * s_val
  tibble(
    Cultivar_Key = row$Cultivar_Key,
    Genome_Group = row$Genome_Group,
    x            = x_seq,
    y            = y_seq
  )
}))

tp_row <- topo_objs$topo_df %>% filter(Pair_ID == "TrackA_vs_Phred")
stat_badge_text <- sprintf(
  "Entanglement E = %.2f   |   Baker's \u03b3 = %.3f (exact p = %.3f)   |   Normalized RF = %.2f",
  tp_row$Entanglement_Step2Side, tp_row$Bakers_Gamma, tp_row$Bakers_Exact_P_TwoSided, tp_row$RF_Unrooted_Normalized
)

p7_tangle <- ggplot() +
  # Left and Right Dendrogram Branches
  geom_segment(data = segs_left,  aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.80, color = "#222222") +
  geom_segment(data = segs_right, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.80, color = "#222222") +
  # Dotted leader lines from label ends (-0.32 / +0.32) to spline anchors (-0.26 / +0.26)
  geom_segment(data = tips_left,  aes(x = -0.32, xend = -0.26, y = y_tip, yend = y_tip), linetype = "dotted", color = "#888888", linewidth = 0.45) +
  geom_segment(data = tips_right, aes(x =  0.26, xend =  0.32, y = y_tip, yend = y_tip), linetype = "dotted", color = "#888888", linewidth = 0.45) +
  # Smooth Cophylogenetic Connecting Splines
  geom_path(
    data = sigmoid_curves,
    aes(x = x, y = y, group = Cultivar_Key, color = Genome_Group),
    linewidth = 1.15, alpha = 0.90, show.legend = FALSE
  ) +
  # Anchor dots at spline endpoints
  geom_point(data = links_df, aes(x = -0.26, y = y_left,  color = Genome_Group), size = 2.2, show.legend = FALSE) +
  geom_point(data = links_df, aes(x =  0.26, y = y_right, color = Genome_Group), size = 2.2, show.legend = FALSE) +
  # Dendrogram Tip Points Colored by Genome Group
  geom_point(data = tips_left,  aes(x = x_tip, y = y_tip, color = Genome_Group), shape = 19, size = 3.8) +
  geom_point(data = tips_right, aes(x = x_tip, y = y_tip, color = Genome_Group), shape = 19, size = 3.8, show.legend = FALSE) +
  # Black Italicized Cultivar Labels on Left and Right
  ggtext::geom_richtext(
    data = tips_left,
    aes(x = x_tip + 0.03, y = y_tip, label = Clean_Label_Markdown),
    color = "#111111", hjust = 0, size = 2.75, fontface = "bold",
    fill = NA, label.color = NA, show.legend = FALSE
  ) +
  ggtext::geom_richtext(
    data = tips_right,
    aes(x = x_tip - 0.03, y = y_tip, label = Clean_Label_Markdown),
    color = "#111111", hjust = 1, size = 2.75, fontface = "bold",
    fill = NA, label.color = NA, show.legend = FALSE
  ) +
  # Directional Tree Headers at Top of Corridor
  annotate("text", x = -1.54, y = 6.55, label = "Morphological Dendrogram\n(VIF-Pruned Euclidean + UPGMA)", size = 2.95, fontface = "bold", color = "#111111") +
  annotate("text", x =  1.54, y = 6.55, label = "Molecular ITS Dendrogram\n(Phred-Weighted TrN+G4 ML)", size = 2.95, fontface = "bold", color = "#111111") +
  annotate("label", x = 0, y = 0.35, label = stat_badge_text, size = 2.75, fontface = "bold", fill = "#F8F9FA", color = "#111111", linewidth = 0.3) +
  scale_color_manual(
    values = GENOME_PALETTE,
    drop   = TRUE,
    name   = "*Musa* Genome Group"
  ) +
  scale_x_continuous(
    limits = c(-1.98, 1.98),
    breaks = c(-1.90, -1.54, -1.18, 1.18, 1.54, 1.90),
    labels = c("1.0", "0.5", "0.0", "0.0", "0.5", "1.0")
  ) +
  scale_y_continuous(limits = c(0.15, 6.80), breaks = NULL) +
  labs(
    x = "Normalized Ultrametric Branch Height (Left: Morphological Distance | Right: ITS Genetic Distance)",
    y = NULL
  ) +
  theme_thesis(base_size = 9.5) +
  theme(
    panel.grid.major = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.position  = "bottom"
  )

# --- Panel B: Exact 6! = 720 Permutation Null Distribution of Baker's Gamma ---
null_prim <- topo_objs$null_df %>% filter(Pair_ID == "TrackA_vs_Phred")
bg_obs_val <- tp_row$Bakers_Gamma
ci_null    <- quantile(null_prim$Bakers_Gamma, probs = c(0.025, 0.975))

p7_null <- ggplot(null_prim, aes(x = Bakers_Gamma)) +
  annotate(
    "rect",
    xmin = ci_null[1], xmax = ci_null[2], ymin = 0, ymax = Inf,
    fill = "#E0E0E0", alpha = 0.45
  ) +
  geom_histogram(aes(y = after_stat(density)), bins = 18, fill = "#56B4E9", color = "#111111", linewidth = 0.35, alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "#444444", linewidth = 0.6) +
  geom_vline(xintercept = bg_obs_val, linetype = "solid", color = "#D55E00", linewidth = 1.1) +
  annotate(
    "label",
    x = ifelse(bg_obs_val >= 0, min(0.65, bg_obs_val + 0.05), max(-0.65, bg_obs_val - 0.05)),
    y = 1.85,
    label = sprintf("Observed \u03b3 = %.3f\nExact p = %.3f\n(6! = 720 Permutations)", bg_obs_val, tp_row$Bakers_Exact_P_TwoSided),
    size = 2.8, fontface = "bold", fill = "#FFF5EE", color = "#D55E00", linewidth = 0.35
  ) +
  scale_x_continuous(limits = c(-1.05, 1.05), breaks = seq(-1.0, 1.0, by = 0.5)) +
  labs(
    x = "Baker's Gamma Cophenetic Rank Correlation (\u03b3)",
    y = "Exact Permutation Null Density"
  ) +
  theme_thesis(base_size = 9.5)

fig7_combined <- (p7_tangle | p7_null) +
  plot_layout(widths = c(1.65, 1.0)) +
  plot_annotation(tag_levels = "A")

save_thesis_plot(fig7_combined, "Fig7_Cophylogenetic_Tanglegram_Topological_Benchmark", width = 14.5, height = 6.5, dpi = 600)

# ==============================================================================
# FIGURE 8: Mantel Distance Concordance (A) & Single-Trait Decomposition (B)
# ==============================================================================

pw_df <- conc_objs$pairwise_df
m_prim <- conc_objs$mantel_specs %>%
  filter(
    Morphological_Distance == "Track A: Continuous VIF-Pruned Euclidean (6 Traits)",
    Molecular_Distance     == "Phred-Weighted ML Genetic Distance (TrN+G4)",
    Test_Type              == "Standard Mantel"
  )
m_part <- conc_objs$mantel_specs %>%
  filter(
    Morphological_Distance == "Track A: Continuous VIF-Pruned Euclidean (6 Traits)",
    Molecular_Distance     == "Phred-Weighted ML Genetic Distance (TrN+G4)",
    Test_Type              == "Partial Mantel (| Genome Group)"
  )

mantel_annot_txt <- sprintf(
  "Standard Mantel: r_M = %.3f (95%% CI [%.2f, %.2f], exact p = %.3f)\nSpearman Rank: \u03c1_M = %.3f (exact p = %.3f)\nPartial Mantel (| Genome Group): r_M = %.3f (exact p = %.3f)",
  m_prim$Mantel_r_Pearson, m_prim$Pearson_CI_Low, m_prim$Pearson_CI_High, m_prim$Pearson_P_TwoSided,
  m_prim$Mantel_rho_Spearman, m_prim$Spearman_P_TwoSided,
  m_part$Mantel_r_Pearson, m_part$Pearson_P_TwoSided
)

cross_palette <- c(
  "Intra-Acuminata Diploid Cross (AA x AA/AAB)"              = "#E69F00",
  "Intra-Acuminata Ploidy Bridge (Diploid AA x Triploid AAA)" = "#0072B2",
  "Interspecific (Acuminata AA/AAA x Balbisiana ABB)"        = "#CC79A7"
)

p8_mantel <- ggplot(pw_df, aes(x = ITS_ML_Dist, y = Morph_Euclidean_VIF)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#333333", fill = "#D9D9D9", linetype = "dashed", linewidth = 0.75, alpha = 0.45) +
  geom_point(aes(color = Cross_Category, shape = Cross_Category), size = 3.8, alpha = 0.92) +
  ggrepel::geom_text_repel(
    aes(label = Pair_Label),
    color = "#111111", size = 2.75, fontface = "bold",
    box.padding = 0.35, point.padding = 0.25, max.overlaps = 20, seed = 20261102
  ) +
  annotate(
    "label",
    x = min(pw_df$ITS_ML_Dist) + 0.02 * diff(range(pw_df$ITS_ML_Dist)),
    y = max(pw_df$Morph_Euclidean_VIF) + 0.12 * diff(range(pw_df$Morph_Euclidean_VIF)),
    label = mantel_annot_txt,
    hjust = 0, vjust = 1, size = 2.65, fontface = "bold",
    fill = "#F8F9FA", color = "#111111", linewidth = 0.3
  ) +
  scale_color_manual(values = cross_palette, name = "Parental Cross Architecture") +
  scale_shape_manual(
    values = c(
      "Intra-Acuminata Diploid Cross (AA x AA/AAB)"              = 19,
      "Intra-Acuminata Ploidy Bridge (Diploid AA x Triploid AAA)" = 17,
      "Interspecific (Acuminata AA/AAA x Balbisiana ABB)"        = 15
    ),
    name = "Parental Cross Architecture"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  labs(
    x = "Pairwise ITS Maximum Likelihood Genetic Distance (Substitutions / Site)",
    y = "Pairwise Morphological Euclidean Distance (Z-Standardized, VIF < 5)"
  ) +
  theme_thesis(base_size = 9.2) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical"
  )

# --- Panel B: Character-by-Character Mantel Decomposition ---------------------
tm_df <- conc_objs$trait_mantel_df %>%
  arrange(Mantel_r_Pearson) %>%
  mutate(
    Trait_Label = factor(Trait_Label, levels = unique(Trait_Label)),
    Label_Txt   = sprintf("r = %+.2f (p = %.2f)", Mantel_r_Pearson, Pearson_P_TwoSided)
  )

trait_cat_pal <- c(
  "Continuous (Retained, VIF < 5)"       = "#0072B2",
  "Continuous (Collinear, VIF >= 5)"     = "#E69F00",
  "Qualitative IPGRI Descriptor (Gower)" = "#009E73"
)

p8_traits <- ggplot(tm_df, aes(x = Mantel_r_Pearson, y = Trait_Label, color = Trait_Category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#444444", linewidth = 0.6) +
  geom_segment(aes(x = 0, xend = Mantel_r_Pearson, y = Trait_Label, yend = Trait_Label), linewidth = 0.85) +
  geom_point(size = 3.5) +
  geom_text(
    aes(
      x     = ifelse(Mantel_r_Pearson >= 0, Mantel_r_Pearson + 0.03, Mantel_r_Pearson - 0.03),
      hjust = ifelse(Mantel_r_Pearson >= 0, 0, 1),
      label = Label_Txt
    ),
    color = "#111111", size = 2.60, fontface = "bold", show.legend = FALSE
  ) +
  scale_color_manual(values = trait_cat_pal, name = "Morphological Character Class") +
  scale_x_continuous(limits = c(-0.65, 1.15), breaks = seq(-0.6, 0.9, by = 0.3)) +
  labs(
    x = "Single-Trait Mantel Correlation with ITS Genetic Distance (r_M)",
    y = NULL
  ) +
  theme_thesis(base_size = 9.2) +
  theme(
    legend.position  = "bottom",
    legend.direction = "vertical"
  )

fig8_combined <- (p8_mantel | p8_traits) +
  plot_layout(widths = c(1.15, 1.05)) +
  plot_annotation(tag_levels = "A")

save_thesis_plot(fig8_combined, "Fig8_Mantel_Regression_Trait_Decomposition", width = 14.5, height = 6.8, dpi = 600)

# ==============================================================================
# FIGURE 9: Symmetric Procrustes Superimposition (A) & Residual Audit (B)
# ==============================================================================

proc_df <- conc_objs$proc_cultivar_df %>%
  left_join(shared_meta %>% select(Cultivar_Key, Clean_Label_Markdown), by = "Cultivar_Key") %>%
  mutate(Genome_Group = factor(Genome_Group, levels = names(GENOME_PALETTE)))

proc_pts_long <- bind_rows(
  proc_df %>% transmute(Cultivar_Key, Genome_Group, Dim1 = ITS_Target_X1,    Dim2 = ITS_Target_X2,    Space = "Molecular ITS PCoA (Target)"),
  proc_df %>% transmute(Cultivar_Key, Genome_Group, Dim1 = Morph_Rotated_Y1, Dim2 = Morph_Rotated_Y2, Space = "Morphological PCoA (Procrustes-Rotated)")
)

proc_badge_txt <- sprintf(
  "Symmetric Procrustes: M\u00b2 = %.3f   |   r_P = %.3f (exact p = %.3f)",
  conc_objs$m2_obs, conc_objs$rp_obs, conc_objs$p_proc_exact
)

p9_proc <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "#888888", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "#888888", linewidth = 0.5) +
  # Displacement vectors from Molecular ITS Target (X) to Rotated Morphology (Yrot)
  geom_segment(
    data = proc_df,
    aes(x = ITS_Target_X1, y = ITS_Target_X2, xend = Morph_Rotated_Y1, yend = Morph_Rotated_Y2, color = Genome_Group),
    arrow = grid::arrow(length = unit(0.20, "cm"), type = "closed"),
    linewidth = 0.85, alpha = 0.88, show.legend = FALSE
  ) +
  geom_point(
    data = proc_pts_long,
    aes(x = Dim1, y = Dim2, color = Genome_Group, shape = Space),
    size = 3.8, alpha = 0.95
  ) +
  # Place cultivar label near each cultivar's Procrustes vector with clean vertical offsets
  ggtext::geom_richtext(
    data = proc_df %>% mutate(
      lbl_x = case_when(
        Cultivar_Key == "Mas"   ~ -0.21,
        Cultivar_Key == "Mauli" ~ -0.25,
        Cultivar_Key == "Talas" ~ -0.09,
        Cultivar_Key == "Awak"  ~  0.02,
        TRUE                    ~  0.5 * (ITS_Target_X1 + Morph_Rotated_Y1)
      ),
      lbl_y = case_when(
        Cultivar_Key == "Mas"       ~  0.235,
        Cultivar_Key == "Mauli"     ~ -0.085,
        Cultivar_Key == "Talas"     ~ -0.175,
        Cultivar_Key == "Awak"      ~ -0.355,
        Cultivar_Key == "Kalimbang" ~  0.390,
        Cultivar_Key == "Ambon"     ~  0.065,
        TRUE                        ~  0.5 * (ITS_Target_X2 + Morph_Rotated_Y2)
      )
    ),
    aes(x = lbl_x, y = lbl_y, label = Clean_Label_Markdown),
    color = "#111111", size = 2.95, fontface = "bold",
    fill = alpha("white", 0.78), label.color = NA, show.legend = FALSE
  ) +
  annotate(
    "label",
    x = min(c(proc_df$ITS_Target_X1, proc_df$Morph_Rotated_Y1)) - 0.02,
    y = max(c(proc_df$ITS_Target_X2, proc_df$Morph_Rotated_Y2)) + 0.09,
    label = proc_badge_txt,
    hjust = 0, vjust = 1, size = 2.75, fontface = "bold",
    fill = "#F8F9FA", color = "#111111", linewidth = 0.3
  ) +
  scale_color_manual(
    values = GENOME_PALETTE,
    drop   = TRUE,
    name   = "*Musa* Genome Group"
  ) +
  scale_shape_manual(
    values = c(
      "Molecular ITS PCoA (Target)"             = 19,
      "Morphological PCoA (Procrustes-Rotated)" = 17
    ),
    name = "Ordination Configuration"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.18))) +
  guides(
    color = guide_legend(order = 1, nrow = 1, override.aes = list(size = 3.8)),
    shape = guide_legend(order = 2, nrow = 1, override.aes = list(size = 3.6, color = "#111111"))
  ) +
  labs(
    x = "Procrustes Superimposition Dimension 1",
    y = "Procrustes Superimposition Dimension 2"
  ) +
  theme_thesis(base_size = 9.2) +
  theme(
    legend.position = "bottom",
    legend.box      = "vertical"
  )

# --- Panel B: Cultivar-Specific Procrustes Residual Magnitude Ranking ---------
q50_val <- median(proc_df$Procrustes_Residual)
q75_val <- quantile(proc_df$Procrustes_Residual, 0.75)

proc_rank_df <- proc_df %>%
  arrange(Procrustes_Residual) %>%
  mutate(
    Clean_Label_Markdown = factor(Clean_Label_Markdown, levels = unique(Clean_Label_Markdown)),
    Res_Label            = sprintf("\u03b5 = %.3f (%.1f%% of M\u00b2)", Procrustes_Residual, Pct_Contribution_M2)
  )

p9_res <- ggplot(proc_rank_df, aes(x = Procrustes_Residual, y = Clean_Label_Markdown)) +
  geom_vline(xintercept = q50_val, linetype = "dotted", color = "#555555", linewidth = 0.6) +
  geom_vline(xintercept = q75_val, linetype = "dashed", color = "#D55E00", linewidth = 0.65) +
  geom_segment(aes(x = 0, xend = Procrustes_Residual, y = Clean_Label_Markdown, yend = Clean_Label_Markdown), color = "#444444", linewidth = 0.85) +
  geom_point(aes(color = Genome_Group), size = 4.0, show.legend = FALSE) +
  geom_text(
    aes(x = Procrustes_Residual + 0.012, label = Res_Label),
    hjust = 0, size = 2.75, fontface = "bold", color = "#111111"
  ) +
  annotate("text", x = q50_val - 0.005, y = 6.45, label = "Median (Q50)", hjust = 1, size = 2.5, fontface = "italic", color = "#444444") +
  annotate("text", x = q75_val + 0.005, y = 6.45, label = "Outlier (Q75)", hjust = 0, size = 2.5, fontface = "bold", color = "#D55E00") +
  scale_color_manual(values = GENOME_PALETTE, drop = TRUE) +
  scale_x_continuous(limits = c(0, max(proc_rank_df$Procrustes_Residual) * 1.45)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.85))) +
  labs(
    x = "Cultivar Procrustes Residual Magnitude (\u03b5_i)",
    y = NULL
  ) +
  theme_thesis(base_size = 9.2) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 9.0, color = "#111111")
  )

fig9_combined <- (p9_proc | p9_res) +
  plot_layout(widths = c(1.25, 1.0)) +
  plot_annotation(tag_levels = "A")

save_thesis_plot(fig9_combined, "Fig9_Procrustes_Superimposition_Residual_Audit", width = 14.5, height = 6.6, dpi = 600)

message("[11] All Phase 3 comparative concordance figures generated successfully.")


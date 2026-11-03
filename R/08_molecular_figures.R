# ==============================================================================
# R/08_molecular_figures.R
# Phase 2 - Step 4: Report-Ready Molecular Phylogenetics Figures
# (Title/Subtitle-Free, Black Italicized Labels, Locked Okabe-Ito Genome Group
# Tip Points, SH-aLRT / ML-BS Node Badges, Vector PDF + 600 DPI PNG)
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(ggtree)
  library(ape)
})
reset_project_seed()

message("=== [08] Generating Phase 2 Molecular Publication Figures ===")

cons_objs  <- readRDS(file.path(DIRS$rds, "molecular_consensus_objects.rds"))
aln_objs   <- readRDS(file.path(DIRS$rds, "molecular_alignment_objects.rds"))
phylo_objs <- readRDS(file.path(DIRS$rds, "molecular_phylogeny_objects.rds"))
meta_df    <- read.csv(file.path(DIRS$processed, "cultivar_metadata.csv"), stringsAsFactors = FALSE)

# Build clean display labels without redundant ITS mapping clutter
its_meta <- meta_df %>%
  filter(In_ITS) %>%
  mutate(
    Clean_Label_Markdown = case_when(
      ITS_Accession_ID == "Ensete_ventricosum" ~ "*Ensete ventricosum* (Outgroup)",
      ITS_Accession_ID == "Hutan"              ~ "*Musa acuminata* 'Hutan' (Wild)",
      ITS_Accession_ID == "II_Kalimbang"       ~ "*Musa* 'Kalimbang II' (AA/AAB)",
      ITS_Accession_ID == "Awang_Bangkal"      ~ "*Musa* sp. 'Awang Bangkal' (AA/AAB)",
      TRUE ~ paste0("*Musa* '", Display_Name, "' (", sub(" \\(Local\\)", "", Genome_Group), ")")
    ),
    Clean_Label_Plain = case_when(
      ITS_Accession_ID == "Ensete_ventricosum" ~ "Ensete ventricosum (Outgroup)",
      ITS_Accession_ID == "Hutan"              ~ "Musa acuminata 'Hutan' (Wild)",
      ITS_Accession_ID == "II_Kalimbang"       ~ "Musa 'Kalimbang II' (AA/AAB)",
      ITS_Accession_ID == "Awang_Bangkal"      ~ "Musa sp. 'Awang Bangkal' (AA/AAB)",
      TRUE ~ paste0("Musa '", Display_Name, "' (", sub(" \\(Local\\)", "", Genome_Group), ")")
    )
  )

# ==============================================================================
# FIGURE 4: Consensus Ambiguity Reduction (Panel A) & Alignment Domain Map (B)
# ==============================================================================

qc_df <- cons_objs$qc_df %>%
  left_join(its_meta %>% select(ITS_Accession_ID, Clean_Label_Markdown, Genome_Group), by = "ITS_Accession_ID") %>%
  arrange(Galaxy_Ambig_Pct) %>%
  mutate(Clean_Label_Markdown = factor(Clean_Label_Markdown, levels = unique(Clean_Label_Markdown)))

qc_long <- qc_df %>%
  select(ITS_Accession_ID, Clean_Label_Markdown, Genome_Group,
         `Galaxy Unweighted Consensus` = Galaxy_Ambig_Pct,
         `Phred-Weighted Consensus`    = Phred_Ambig_Pct) %>%
  pivot_longer(cols = c("Galaxy Unweighted Consensus", "Phred-Weighted Consensus"),
               names_to = "Assembly_Method", values_to = "Ambig_Pct")

p4_ambig <- ggplot(qc_long, aes(x = Ambig_Pct, y = Clean_Label_Markdown)) +
  geom_line(aes(group = Clean_Label_Markdown), color = "#666666", linewidth = 0.75) +
  geom_point(aes(color = Assembly_Method, shape = Assembly_Method), size = 3.6, alpha = 0.95) +
  geom_text(
    data = qc_df %>% filter(Galaxy_Ambig_Pct > 1.0),
    aes(x = Galaxy_Ambig_Pct, y = Clean_Label_Markdown,
        label = sprintf("%.1f%% -> %.1f%%", Galaxy_Ambig_Pct, Phred_Ambig_Pct)),
    hjust = -0.15, size = 2.65, fontface = "bold", color = "#222222",
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Galaxy Unweighted Consensus" = "#D55E00",
      "Phred-Weighted Consensus"    = "#0072B2"
    ),
    name = "Consensus Assembly"
  ) +
  scale_shape_manual(
    values = c(
      "Galaxy Unweighted Consensus" = 17,
      "Phred-Weighted Consensus"    = 19
    ),
    name = "Consensus Assembly"
  ) +
  scale_x_continuous(limits = c(0, 30.5), breaks = seq(0, 25, by = 5)) +
  labs(
    x = "IUPAC Ambiguous Bases (% of Amplicon Length)",
    y = NULL
  ) +
  theme_thesis(base_size = 9.5) +
  theme(
    axis.text.y     = ggtext::element_markdown(size = 8.8, color = "#111111"),
    legend.position = "bottom"
  )

# --- Panel B: Alignment Architecture & Trimming Profile across rRNA Domains ---
col_audit <- aln_objs$res_phred$col_audit %>%
  mutate(
    # Smooth Shannon information across a 9-bp sliding window for clean visualization
    Info_Smooth = as.numeric(stats::filter(Info_Content_Bits, rep(1 / 9, 9), sides = 2)),
    Info_Smooth = ifelse(is.na(Info_Smooth), Info_Content_Bits, Info_Smooth),
    Status      = ifelse(Retained, "Retained Column", "Trimmed (Overhang / Gappy)")
  )

domain_bounds <- col_audit %>%
  group_by(rRNA_Domain) %>%
  summarise(
    xmin = min(Column_Index),
    xmax = max(Column_Index),
    xmid = 0.5 * (min(Column_Index) + max(Column_Index)),
    .groups = "drop"
  ) %>%
  mutate(
    Short_Label = case_when(
      rRNA_Domain == "18S rRNA (3')" ~ "18S",
      rRNA_Domain == "26S rRNA (5')" ~ "26S",
      TRUE                           ~ rRNA_Domain
    ),
    Domain_Fill = case_when(
      rRNA_Domain == "18S rRNA (3')" ~ "#E69F00",
      rRNA_Domain == "ITS1"          ~ "#56B4E9",
      rRNA_Domain == "5.8S rRNA"     ~ "#009E73",
      rRNA_Domain == "ITS2"          ~ "#CC79A7",
      TRUE                           ~ "#0072B2"
    )
  )

p4_aln <- ggplot(col_audit, aes(x = Column_Index)) +
  # Information content ribbon colored and filled by Retained vs Trimmed status
  geom_col(aes(y = Info_Content_Bits, fill = Status, color = Status), width = 1.0, alpha = 0.75) +
  geom_line(aes(y = Info_Smooth), color = "#111111", linewidth = 0.65) +
  # Structural domain bar at top (y = 2.06 to 2.24) using explicit fill outside aes()
  geom_rect(
    data = domain_bounds,
    aes(xmin = xmin, xmax = xmax, ymin = 2.06, ymax = 2.24),
    fill = domain_bounds$Domain_Fill,
    color = "#222222", linewidth = 0.35, alpha = 0.88, inherit.aes = FALSE
  ) +
  geom_text(
    data = domain_bounds,
    aes(x = xmid, y = 2.15, label = Short_Label),
    size = 2.65, fontface = "bold", color = "#111111", inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Retained Column"            = "#0072B2",
      "Trimmed (Overhang / Gappy)" = "#D55E00"
    ),
    name = "Alignment Column Filter"
  ) +
  scale_color_manual(
    values = c(
      "Retained Column"            = "#0072B2",
      "Trimmed (Overhang / Gappy)" = "#D55E00"
    ),
    name = "Alignment Column Filter"
  ) +
  scale_y_continuous(limits = c(0, 2.30), breaks = seq(0, 2.0, by = 0.5)) +
  labs(
    x = "Multiple Sequence Alignment Position (bp)",
    y = "Sequence Conservation (Information Bits)"
  ) +
  theme_thesis(base_size = 9.5) +
  theme(legend.position = "bottom")

fig4_combined <- (p4_ambig | p4_aln) +
  plot_layout(widths = c(1.05, 1.30)) +
  plot_annotation(tag_levels = "A")

save_thesis_plot(fig4_combined, "Fig4_ITS_Consensus_Alignment_QC", width = 14.0, height = 6.2, dpi = 600)

# ==============================================================================
# FIGURE 5: Substitution ModelTest (BIC/AICc) & Cophenetic / Support Calibration
# ==============================================================================

mt_top <- phylo_objs$phylo_phred$mt_df %>%
  slice_head(n = 20) %>%
  mutate(
    Is_Winner = row_number() == 1,
    Model     = factor(Model, levels = rev(Model))
  )

mt_long <- mt_top %>%
  select(Model, Is_Winner, Delta_BIC, Delta_AICc) %>%
  pivot_longer(cols = c("Delta_BIC", "Delta_AICc"), names_to = "Criterion", values_to = "Delta_Val") %>%
  mutate(Criterion = ifelse(Criterion == "Delta_BIC", "\u0394BIC", "\u0394AICc"))

p5_model <- ggplot(mt_long, aes(x = Delta_Val, y = Model, color = Criterion, shape = Criterion)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#0072B2", linewidth = 0.6) +
  geom_vline(xintercept = 2, linetype = "dotted", color = "#666666", linewidth = 0.5) +
  geom_line(aes(group = Model), color = "#999999", linewidth = 0.45) +
  geom_point(aes(size = Is_Winner), alpha = 0.92) +
  scale_size_manual(values = c("FALSE" = 2.6, "TRUE" = 4.2), guide = "none") +
  scale_color_manual(
    values = c("\u0394BIC" = "#0072B2", "\u0394AICc" = "#E69F00"),
    name   = "Information Criterion"
  ) +
  scale_shape_manual(
    values = c("\u0394BIC" = 19, "\u0394AICc" = 17),
    name   = "Information Criterion"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3.5)),
    shape = guide_legend(override.aes = list(size = 3.5))
  ) +
  labs(
    x = "Information Criterion Difference from Optimal Model (\u0394BIC / \u0394AICc)",
    y = "Candidate Nucleotide Substitution Model (Top 20 of 92)"
  ) +
  theme_thesis(base_size = 9.0) +
  theme(legend.position = "bottom")

# --- Panel B: ML Pairwise Distance vs. Tree Patristic Distance Calibration ----
cal_df <- bind_rows(
  tibble(
    Assembly = sprintf("Phred-Weighted (%s)\nIngroup r_coph = %.3f",
                       phylo_objs$phylo_phred$best_model_name,
                       phylo_objs$phylo_phred$ingroup_coph_r),
    ML_Dist  = as.numeric(phylo_objs$phylo_phred$d_ml_in),
    Pat_Dist = as.numeric(phylo_objs$phylo_phred$d_pat_in)
  ),
  tibble(
    Assembly = sprintf("Galaxy Unweighted (%s)\nIngroup r_coph = %.3f",
                       phylo_objs$phylo_galaxy$best_model_name,
                       phylo_objs$phylo_galaxy$ingroup_coph_r),
    ML_Dist  = as.numeric(phylo_objs$phylo_galaxy$d_ml_in),
    Pat_Dist = as.numeric(phylo_objs$phylo_galaxy$d_pat_in)
  )
)

p5_cal <- ggplot(cal_df, aes(x = ML_Dist, y = Pat_Dist)) +
  geom_point(color = "#0072B2", size = 2.4, alpha = 0.78) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#D55E00", fill = "#F5C6A5", linewidth = 0.8) +
  facet_wrap(~ Assembly, scales = "free", ncol = 1) +
  labs(
    x = "Pairwise Maximum Likelihood Genetic Distance (Substitutions / Site)",
    y = "ML Tree Patristic Distance (Substitutions / Site)"
  ) +
  theme_thesis(base_size = 9.0)

fig5_combined <- (p5_model | p5_cal) +
  plot_layout(widths = c(1.25, 1.0)) +
  plot_annotation(tag_levels = "A")

save_thesis_plot(fig5_combined, "Fig5_Substitution_ModelTest_Benchmark", width = 13.5, height = 6.8, dpi = 600)

# ==============================================================================
# FIGURE 6: Rooted Maximum Likelihood ITS Phylogram (Full N=11 & Shared N=6)
# ==============================================================================

build_its_phylogram_panel <- function(tr_obj, ns_df, x_label, hexp_right = 0.65, show_leg = TRUE) {
  n_tips  <- length(tr_obj$tip.label)
  n_nodes <- tr_obj$Nnode

  # Match internal nodes by clade membership
  ape_clades <- sapply((n_tips + 1):(n_tips + n_nodes), function(nd) {
    paste(sort(ape::extract.clade(tr_obj, nd)$tip.label), collapse = " + ")
  })

  # Parse node labels directly if present ("SH/BS")
  if (!is.null(tr_obj$node.label) && length(tr_obj$node.label) == n_nodes) {
    parts <- strsplit(tr_obj$node.label, "/")
    sh_v  <- as.numeric(sapply(parts, `[`, 1))
    bs_v  <- as.numeric(sapply(parts, `[`, 2))
  } else {
    m_idx <- match(ape_clades, ns_df$Clade_Taxa)
    sh_v  <- ns_df$SH_aLRT_Pct[m_idx]
    bs_v  <- ns_df$ML_BS_Pct[m_idx]
  }

  support_levels <- c(
    "Strong (SH-aLRT \u2265 80% & BS \u2265 95%)",
    "Moderate (SH-aLRT \u2265 80% or BS \u2265 70%)",
    "Weak / Unresolved (< 80% / < 70%)"
  )

  node_ids <- (n_tips + 1):(n_tips + n_nodes)
  clade_sizes <- sapply(node_ids, function(nd) length(ape::extract.clade(tr_obj, nd)$tip.label))
  # Suppress duplicate root badge at node (n_tips + 1)
  is_root_split <- (node_ids == (n_tips + 1))

  # When Hutan is pruned in the N=6 tree, the (Mas + Mauli + Talas) clade inherits the full
  # AA clade support (Node 17: 91.8% SH-aLRT / 87.2% ML-BS) from the parent branch
  for (k in seq_len(n_nodes)) {
    if (ape_clades[k] %in% c("Mas + Mauli + Talas") && n_tips == 6) {
      aa_row <- ns_df %>% filter(Clade_Taxa == "Hutan + Mas + Mauli + Talas")
      if (nrow(aa_row) == 1) {
        sh_v[k] <- aa_row$SH_aLRT_Pct
        bs_v[k] <- aa_row$ML_BS_Pct
      }
    }
  }

  node_annot <- tibble(
    node        = node_ids,
    Clade_Size  = clade_sizes,
    SH          = ifelse(is_root_split, NA_real_, sh_v),
    BS          = ifelse(is_root_split, NA_real_, bs_v),
    Support_Lbl = ifelse(is_root_split, NA_character_, sprintf("%.0f / %.0f", sh_v, bs_v)),
    Support_Cat = factor(
      ifelse(
        is_root_split, NA_character_,
        case_when(
          sh_v >= 80 & bs_v >= 95 ~ "Strong (SH-aLRT \u2265 80% & BS \u2265 95%)",
          sh_v >= 80 | bs_v >= 70 ~ "Moderate (SH-aLRT \u2265 80% or BS \u2265 70%)",
          TRUE                    ~ "Weak / Unresolved (< 80% / < 70%)"
        )
      ),
      levels = support_levels
    )
  )

  # Map tip metadata (handle both ITS_Accession_ID tips in N=11 and Cultivar_Key tips in N=6)
  tip_annot <- tibble(label = tr_obj$tip.label) %>%
    left_join(
      its_meta %>% select(ITS_Accession_ID, Cultivar_Key, Clean_Label_Markdown, Genome_Group),
      by = c("label" = "ITS_Accession_ID")
    )

  # Fill any tips matched by Cultivar_Key (e.g. Awak in N=6)
  na_idx <- which(is.na(tip_annot$Genome_Group))
  if (length(na_idx) > 0) {
    m_ck <- match(tip_annot$label[na_idx], its_meta$Cultivar_Key)
    tip_annot$Clean_Label_Markdown[na_idx] <- its_meta$Clean_Label_Markdown[m_ck]
    tip_annot$Genome_Group[na_idx]         <- its_meta$Genome_Group[m_ck]
  }
  tip_annot$Genome_Group <- factor(tip_annot$Genome_Group, levels = names(GENOME_PALETTE))

  max_x <- max(ape::node.depth.edgelength(tr_obj))

  p <- ggtree(tr_obj, linewidth = 0.75, color = "#222222") %<+% tip_annot %<+% node_annot +
    geom_nodepoint(
      data = ~ dplyr::filter(.x, !is.na(Support_Cat) & (Clade_Size == 10 | is.na(branch.length) | branch.length > 1e-4)),
      aes(fill = Support_Cat),
      shape = 21, size = 3.0, color = "#111111",
      na.rm = TRUE, show.legend = show_leg
    ) +
    geom_label(
      data = ~ dplyr::filter(.x, !is.na(Support_Cat) & (Clade_Size == 10 | is.na(branch.length) | branch.length > 1e-4)),
      aes(label = Support_Lbl, fill = Support_Cat),
      size = 2.40, fontface = "bold", color = "black",
      label.padding = unit(0.12, "lines"), linewidth = 0.2,
      hjust = 1.12, vjust = -0.35, na.rm = TRUE, show.legend = FALSE
    ) +
    geom_tippoint(
      aes(color = Genome_Group),
      shape = 19, size = 3.8,
      na.rm = TRUE, show.legend = show_leg
    ) +
    ggtext::geom_richtext(
      aes(label = Clean_Label_Markdown),
      color = "#111111",
      hjust = 0, nudge_x = max_x * 0.025,
      size = 3.1, fill = NA, label.color = NA, fontface = "bold",
      na.rm = TRUE, show.legend = FALSE
    ) +
    scale_color_manual(
      values = GENOME_PALETTE,
      limits = names(GENOME_PALETTE),
      drop   = FALSE,
      na.translate = FALSE,
      name   = "*Musa* / Outgroup Clade"
    ) +
    scale_fill_manual(
      values = c(
        "Strong (SH-aLRT \u2265 80% & BS \u2265 95%)"    = "#A8E6CF",
        "Moderate (SH-aLRT \u2265 80% or BS \u2265 70%)" = "#FDFFAB",
        "Weak / Unresolved (< 80% / < 70%)"              = "#FFD3B6"
      ),
      limits = support_levels,
      drop   = FALSE,
      na.translate = FALSE,
      name   = "Branch Support (SH-aLRT / ML-BS %)"
    ) +
    guides(
      color = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 19, size = 4.0)),
      fill  = guide_legend(order = 2, nrow = 1, override.aes = list(shape = 21, size = 3.8, color = "#111111"))
    ) +
    coord_cartesian(clip = "off") +
    hexpand(hexp_right, direction = 1) +
    hexpand(0.08, direction = -1) +
    theme_tree2() +
    labs(x = x_label) +
    theme(
      axis.title.x    = element_text(size = 9.5, face = "bold"),
      legend.position = if (show_leg) "bottom" else "none",
      legend.box      = "vertical",
      legend.title    = ggtext::element_markdown(size = 8.8, face = "bold"),
      legend.text     = element_text(size = 8.2),
      plot.margin     = margin(10, 25, 10, 10)
    )

  p
}

p6_full <- build_its_phylogram_panel(
  phylo_objs$tr_N11,
  phylo_objs$phylo_phred$node_support_df,
  x_label    = sprintf("ML Substitutions / Site (%s, Rooted on Ensete ventricosum)", phylo_objs$phylo_phred$best_model_name),
  hexp_right = 0.65,
  show_leg   = TRUE
)

p6_shared <- build_its_phylogram_panel(
  phylo_objs$tr_N6_synced,
  phylo_objs$phylo_phred$node_support_df,
  x_label    = sprintf("ML Substitutions / Site (Pruned N = 6 Shared Cultivars)"),
  hexp_right = 0.62,
  show_leg   = FALSE
)

fig6_combined <- (p6_full | p6_shared) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(
      legend.position = "bottom",
      legend.box      = "vertical"
    )
  )

save_thesis_plot(fig6_combined, "Fig6_ITS_ML_Phylogram_Rooted", width = 14.5, height = 7.8, dpi = 600)

message("[08] All Phase 2 molecular publication figures generated successfully.")


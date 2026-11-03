# ==============================================================================
# R/03_morpho_clustering_validation.R
# Step 3: Dual-Track Distance & Clustering Benchmark (Cophenetic Correlation
# r_coph + 95% CI + Mantel Test), Multiscale Bootstrap Validation (pvclust,
# nboot = 10,000), and Two-Tier Export (N = 10 Full & N = 6 Shared ITS)
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(cluster)
  library(pvclust)
  library(ape)
  library(dendextend)
  library(vegan)
})
reset_project_seed()

message("=== [03] Starting Dual-Track Distance & Clustering Validation ===")

consensus_qc  <- read.csv(file.path(DIRS$processed, "morpho_consensus_qc.csv"),  stringsAsFactors = FALSE)
consensus_raw <- read.csv(file.path(DIRS$processed, "morpho_consensus_raw.csv"), stringsAsFactors = FALSE)
pca_objs      <- readRDS(file.path(DIRS$rds, "morpho_pca_objects.rds"))

cont_cols     <- pca_objs$cont_cols
retained_cols <- pca_objs$retained_cols

harm_cat_cols <- c(
  "Fruit_Position_Cat",
  "Fruit_Shape_Cat",
  "Fruit_Tip_Cat",
  "Stalk_Surface_Cat",
  "Skin_Young_Cat"
)

all_cat_cols <- c(
  harm_cat_cols,
  "Fruits_Per_Hand_Cat", "Bunch_Position_T1", "Bunch_Shape_T1",
  "Relic_Tip_T2", "Skin_Ripe_T2", "Skin_Adhesion_T2", "Skin_Cracks_T2",
  "Pulp_Presence_T2", "Pulp_Unripe_T2", "Flesh_Color_T2", "Flesh_Texture_T2",
  "Dominant_Taste_T2", "Seed_Surface_T2", "Seed_Shape_T2",
  "Fruits_Bunch_Cat_T2", "Cross_Section_T2"
)

# 1. Helper to Build Candidate Distance Matrices -------------------------------
prepare_mixed_df <- function(df, cont_vars, cat_vars) {
  sub_df <- df[, c(cont_vars, cat_vars), drop = FALSE]
  rownames(sub_df) <- df$Cultivar_Key
  for (cv in cont_vars) {
    sub_df[[cv]] <- as.numeric(scale(sub_df[[cv]]))
  }
  for (kv in cat_vars) {
    sub_df[[kv]] <- as.factor(sub_df[[kv]])
  }
  sub_df
}

build_distance_candidates <- function(df) {
  X_full <- scale(as.matrix(df[, cont_cols]))
  rownames(X_full) <- df$Cultivar_Key

  X_vif  <- X_full[, retained_cols, drop = FALSE]

  mixed_12 <- prepare_mixed_df(df, cont_cols, harm_cat_cols)
  mixed_vif11 <- prepare_mixed_df(df, retained_cols, harm_cat_cols)
  # Block-balanced weights: 70% continuous morphometrics (6 VIF-pruned traits), 30% categorical IPGRI (5 traits)
  w_block <- c(rep(0.70 / length(retained_cols), length(retained_cols)),
               rep(0.30 / length(harm_cat_cols), length(harm_cat_cols)))
  mixed_28 <- prepare_mixed_df(df, c(cont_cols, "Fruit_Weight_g", "Hands_Per_Bunch"), all_cat_cols)

  list(
    "Euclidean_Full (7 Cont)" = list(
      track   = "Track A: Continuous Morphometrics",
      dist    = dist(X_full, method = "euclidean"),
      mat     = X_full,
      weights = rep(1, ncol(X_full))
    ),
    "Euclidean_VIF (6 Cont, VIF<5)" = list(
      track   = "Track A: Continuous Morphometrics",
      dist    = dist(X_vif, method = "euclidean"),
      mat     = X_vif,
      weights = rep(1, ncol(X_vif))
    ),
    "Manhattan_Full (7 Cont)" = list(
      track   = "Track A: Continuous Morphometrics",
      dist    = dist(X_full, method = "manhattan"),
      mat     = X_full,
      weights = rep(1, ncol(X_full))
    ),
    "Gower_BlockBalanced (11 VIF+IPGRI)" = list(
      track   = "Track B: Mixed Continuous + IPGRI Categorical",
      dist    = cluster::daisy(mixed_vif11, metric = "gower", weights = w_block),
      mat     = mixed_vif11,
      weights = w_block
    ),
    "Gower_Unweighted (12 Mixed, 0% NA)" = list(
      track   = "Track B: Mixed Continuous + IPGRI Categorical",
      dist    = cluster::daisy(mixed_12, metric = "gower"),
      mat     = mixed_12,
      weights = rep(1, ncol(mixed_12))
    ),
    "Gower_Extended (28 Mixed, Partial T1/T2)" = list(
      track   = "Track B: Mixed Continuous + IPGRI Categorical",
      dist    = cluster::daisy(mixed_28, metric = "gower"),
      mat     = mixed_28,
      weights = rep(1, ncol(mixed_28))
    )
  )
}

dist_qc  <- build_distance_candidates(consensus_qc)
dist_raw <- build_distance_candidates(consensus_raw)

# 2. Benchmark Linkage Methods via Cophenetic Correlation (r_coph) -------------
linkage_methods <- c("average", "ward.D2", "complete", "single")
linkage_labels  <- c(
  "average"  = "UPGMA (Average)",
  "ward.D2"  = "Ward.D2",
  "complete" = "Complete Linkage",
  "single"   = "Single Linkage"
)

evaluate_clustering_grid <- function(dist_list, dataset_name) {
  reset_project_seed()
  results <- list()

  for (d_name in names(dist_list)) {
    d_obj   <- dist_list[[d_name]]$dist
    t_label <- dist_list[[d_name]]$track
    d_vec   <- as.numeric(d_obj)

    for (m_code in linkage_methods) {
      hc   <- hclust(d_obj, method = m_code)
      coph <- cophenetic(hc)
      c_vec <- as.numeric(coph)

      # Parametric Pearson & Non-Parametric Spearman Cophenetic Correlation
      ct_p <- cor.test(d_vec, c_vec, method = "pearson")
      ct_s <- suppressWarnings(cor.test(d_vec, c_vec, method = "spearman", exact = FALSE))

      # Non-parametric bootstrap 95% CI on pairwise distance pairs (B = 2,000)
      boot_r <- replicate(2000, {
        idx <- sample.int(length(d_vec), replace = TRUE)
        cor(d_vec[idx], c_vec[idx], method = "pearson")
      })
      ci_boot <- quantile(boot_r, probs = c(0.025, 0.975), na.rm = TRUE)

      # Two-sided Mantel permutation test (720 to 9,999 permutations)
      mt <- quiet_mantel <- suppressMessages(
        vegan::mantel(d_obj, coph, method = "pearson", permutations = 9999)
      )

      # Normalized 2-norm stress after optimal scalar projection
      scale_fac <- sum(d_vec * c_vec) / sum(c_vec^2)
      stress_2  <- sqrt(sum((d_vec - scale_fac * c_vec)^2) / sum(d_vec^2))

      results[[length(results) + 1]] <- tibble(
        Dataset          = dataset_name,
        Track            = t_label,
        Distance_Metric  = d_name,
        Linkage_Code     = m_code,
        Linkage_Method   = linkage_labels[[m_code]],
        r_coph_Pearson   = round(unname(ct_p$estimate), 4),
        CI95_Low         = round(unname(ci_boot[1]), 4),
        CI95_High        = round(unname(ci_boot[2]), 4),
        r_coph_Spearman  = round(unname(ct_s$estimate), 4),
        Stress_2Norm     = round(stress_2, 4),
        Mantel_p         = mt$signif,
        Meets_080_Target = unname(ct_p$estimate) >= 0.80
      )
    }
  }
  bind_rows(results)
}

bench_df <- bind_rows(
  evaluate_clustering_grid(dist_qc,  "QC_Corrected"),
  evaluate_clustering_grid(dist_raw, "Raw_Uncorrected")
) %>%
  mutate(
    Mantel_p_adj_BH = p.adjust(Mantel_p, method = "BH")
  ) %>%
  arrange(Dataset, Track, desc(r_coph_Pearson))

write.csv(bench_df, file.path(DIRS$tables, "clustering_cophenetic_benchmark.csv"), row.names = FALSE)

# 3. Select Optimal Track A & Track B Models Programmatically ------------------
best_trackA <- bench_df %>%
  filter(Dataset == "QC_Corrected", grepl("^Track A", Track)) %>%
  slice_max(order_by = r_coph_Pearson, n = 1, with_ties = FALSE)

best_trackB <- bench_df %>%
  filter(Dataset == "QC_Corrected", grepl("^Track B", Track)) %>%
  slice_max(order_by = r_coph_Pearson, n = 1, with_ties = FALSE)

message(sprintf("[Best Track A] %s + %s: r_coph = %.4f (95%% CI [%.4f, %.4f], Mantel p = %.4f)",
                best_trackA$Distance_Metric, best_trackA$Linkage_Method,
                best_trackA$r_coph_Pearson, best_trackA$CI95_Low, best_trackA$CI95_High, best_trackA$Mantel_p))

message(sprintf("[Best Track B] %s + %s: r_coph = %.4f (95%% CI [%.4f, %.4f], Mantel p = %.4f)",
                best_trackB$Distance_Metric, best_trackB$Linkage_Method,
                best_trackB$r_coph_Pearson, best_trackB$CI95_Low, best_trackB$CI95_High, best_trackB$Mantel_p))

stopifnot(best_trackA$r_coph_Pearson >= 0.80)

# Also record QC sensitivity summary comparing QC_Corrected vs Raw_Uncorrected
qc_sens_summary <- bench_df %>%
  filter(Linkage_Code == "average") %>%
  select(Dataset, Track, Distance_Metric, Linkage_Method, r_coph_Pearson, CI95_Low, CI95_High, Stress_2Norm) %>%
  pivot_wider(
    names_from  = Dataset,
    values_from = c(r_coph_Pearson, CI95_Low, CI95_High, Stress_2Norm)
  ) %>%
  mutate(
    Delta_r_coph = round(r_coph_Pearson_QC_Corrected - r_coph_Pearson_Raw_Uncorrected, 4),
    Delta_Stress = round(Stress_2Norm_QC_Corrected - Stress_2Norm_Raw_Uncorrected, 4)
  )
write.csv(qc_sens_summary, file.path(DIRS$tables, "qc_sensitivity_comparison.csv"), row.names = FALSE)

# 4. Multiscale Bootstrap Resampling via pvclust (nboot = 10,000) --------------
# Note: pvclust clusters COLUMNS, so we transpose (traits = rows, cultivars = cols)
reset_project_seed()
N_BOOT <- 10000L

# Track A: Continuous traits (standardized)
mat_trackA <- dist_qc[[best_trackA$Distance_Metric]]$mat
dist_method_A <- ifelse(grepl("^Manhattan", best_trackA$Distance_Metric), "manhattan", "euclidean")

message(sprintf("[pvclust] Running %d multiscale bootstrap replicates for Track A (%s + %s)...",
                N_BOOT, best_trackA$Distance_Metric, best_trackA$Linkage_Method))

pv_trackA <- pvclust::pvclust(
  data          = t(mat_trackA),
  method.hclust = best_trackA$Linkage_Code,
  method.dist   = dist_method_A,
  nboot         = N_BOOT,
  r             = seq(0.5, 1.4, by = 0.1),
  quiet         = TRUE
)

# Track B: Mixed traits via custom Gower distance function on resampled rows (traits)
# Encode factor categories as integer codes and record which rows are categorical vs continuous
mixed_df_B <- dist_qc[[best_trackB$Distance_Metric]]$mat
weights_B  <- dist_qc[[best_trackB$Distance_Metric]]$weights
names(weights_B) <- colnames(mixed_df_B)
is_cat_trait <- sapply(mixed_df_B, is.factor)
mat_B_num <- sapply(mixed_df_B, function(col) as.numeric(col))
rownames(mat_B_num) <- rownames(mixed_df_B)

# Custom Gower distance across columns (cultivars) when rows (traits) are resampled by pvclust
gower_pvclust_dist <- function(x_traits_by_cultivars) {
  # x_traits_by_cultivars has resampled traits in rows, cultivars in columns
  n_cult <- ncol(x_traits_by_cultivars)
  n_tr   <- nrow(x_traits_by_cultivars)
  rn     <- rownames(x_traits_by_cultivars)
  # Identify which resampled rows came from categorical traits by matching row names
  base_rn   <- sub("\\.[0-9]+$", "", rn)
  cat_flags <- is_cat_trait[base_rn]
  row_wts   <- weights_B[base_rn]
  row_wts[is.na(row_wts)] <- 1.0

  d_accum <- matrix(0, n_cult, n_cult)
  w_accum <- matrix(0, n_cult, n_cult)

  for (r_idx in seq_len(n_tr)) {
    v <- x_traits_by_cultivars[r_idx, ]
    wt <- row_wts[r_idx]
    valid <- !is.na(v)
    w_pair <- outer(valid, valid, `&`)
    if (sum(w_pair) == 0) next

    if (isTRUE(cat_flags[r_idx])) {
      # Nominal/Categorical mismatch (0 if equal, 1 if different)
      diff_mat <- outer(v, v, function(a, b) ifelse(a == b, 0, 1))
    } else {
      # Continuous range-normalized Manhattan distance
      rng <- max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
      if (is.na(rng) || rng < 1e-12) {
        diff_mat <- matrix(0, n_cult, n_cult)
      } else {
        diff_mat <- abs(outer(v, v, `-`)) / rng
      }
    }
    diff_mat[!w_pair] <- 0
    d_accum <- d_accum + (diff_mat * wt)
    w_accum <- w_accum + (w_pair * wt)
  }

  out_mat <- ifelse(w_accum > 0, d_accum / w_accum, 0)
  colnames(out_mat) <- colnames(x_traits_by_cultivars)
  rownames(out_mat) <- colnames(x_traits_by_cultivars)
  as.dist(out_mat)
}

reset_project_seed()
message(sprintf("[pvclust] Running %d multiscale bootstrap replicates for Track B (%s + %s)...",
                N_BOOT, best_trackB$Distance_Metric, best_trackB$Linkage_Method))

pv_trackB <- pvclust::pvclust(
  data          = t(mat_B_num),
  method.hclust = best_trackB$Linkage_Code,
  method.dist   = gower_pvclust_dist,
  nboot         = N_BOOT,
  r             = seq(0.5, 1.4, by = 0.1),
  quiet         = TRUE
)

# Extract node-by-node bootstrap support tables
extract_pv_table <- function(pv_obj, track_label) {
  edges <- pv_obj$edges
  hc    <- pv_obj$hclust
  n_nodes <- nrow(edges)

  # Map each internal node to its descendant tip labels
  phylo_tr <- ape::as.phylo(hc)
  clade_members <- sapply(seq_len(n_nodes), function(i) {
    # Internal node in ape is Ntip + i (ordered by hc$merge)
    m_row <- hc$merge[i, ]
    get_tips <- function(idx) {
      if (idx < 0) return(hc$labels[-idx])
      c(Recall(hc$merge[idx, 1]), Recall(hc$merge[idx, 2]))
    }
    paste(sort(c(get_tips(m_row[1]), get_tips(m_row[2]))), collapse = " + ")
  })

  tibble(
    Track          = track_label,
    Node_ID        = seq_len(n_nodes),
    Branch_Height  = round(hc$height, 4),
    Clade_Taxa     = clade_members,
    AU_Pct         = round(edges[, "au"] * 100, 2),
    BP_Pct         = round(edges[, "bp"] * 100, 2),
    SE_AU_Pct      = round(edges[, "se.au"] * 100, 2),
    SE_BP_Pct      = round(edges[, "se.bp"] * 100, 2),
    AU_p_value     = round(edges[, "pchi"], 4),
    Significant_AU95 = (edges[, "au"] >= 0.95)
  )
}

pv_support_df <- bind_rows(
  extract_pv_table(pv_trackA, paste0("Track A (", best_trackA$Distance_Metric, " + ", best_trackA$Linkage_Method, ")")),
  extract_pv_table(pv_trackB, paste0("Track B (", best_trackB$Distance_Metric, " + ", best_trackB$Linkage_Method, ")"))
)
write.csv(pv_support_df, file.path(DIRS$tables, "pvclust_node_support.csv"), row.names = FALSE)

# 5. Compare Track A vs Track B Tree Concordance (Mantel & Baker's Gamma) ------
D_A <- dist_qc[[best_trackA$Distance_Metric]]$dist
D_B <- dist_qc[[best_trackB$Distance_Metric]]$dist

mantel_AB <- suppressMessages(vegan::mantel(D_A, D_B, method = "spearman", permutations = 9999))
dend_A    <- as.dendrogram(pv_trackA$hclust)
dend_B    <- as.dendrogram(pv_trackB$hclust)
bakers_AB <- dendextend::cor_bakers_gamma(dend_A, dend_B)

message(sprintf("[Track A vs Track B Concordance] Mantel Spearman r_M = %.4f (p = %.4f); Baker's Gamma = %.4f",
                mantel_AB$statistic, mantel_AB$signif, bakers_AB))

# 6. Export Tier 1 (N=10 Full) and Tier 2 (N=6 Shared ITS) Trees & Matrices ----
shared_cultivars <- consensus_qc %>%
  filter(In_ITS) %>%
  pull(Cultivar_Key)

stopifnot(length(shared_cultivars) == 6)

subset_dist_to_shared <- function(d_obj, shared_keys, meta_sub) {
  m <- as.matrix(d_obj)[shared_keys, shared_keys]
  # Also provide an ITS-accession-named copy (e.g. Awak -> Awa) for direct 1:1 matching
  its_names <- meta_sub$ITS_Accession_ID[match(shared_keys, meta_sub$Cultivar_Key)]
  m_its <- m
  rownames(m_its) <- its_names
  colnames(m_its) <- its_names
  list(by_cultivar_key = as.dist(m), by_its_accession = as.dist(m_its))
}

D_A_shared <- subset_dist_to_shared(D_A, shared_cultivars, consensus_qc)
D_B_shared <- subset_dist_to_shared(D_B, shared_cultivars, consensus_qc)

phylo_A_N10 <- ape::as.phylo(pv_trackA$hclust)
phylo_B_N10 <- ape::as.phylo(pv_trackB$hclust)

phylo_A_N6  <- ape::keep.tip(phylo_A_N10, shared_cultivars)
phylo_B_N6  <- ape::keep.tip(phylo_B_N10, shared_cultivars)

ape::write.tree(phylo_A_N10, file = file.path(DIRS$trees, "morpho_trackA_continuous_N10.nwk"))
ape::write.tree(phylo_B_N10, file = file.path(DIRS$trees, "morpho_trackB_gower_N10.nwk"))
ape::write.tree(phylo_A_N6,  file = file.path(DIRS$trees, "morpho_trackA_continuous_N6_shared.nwk"))
ape::write.tree(phylo_B_N6,  file = file.path(DIRS$trees, "morpho_trackB_gower_N6_shared.nwk"))

saveRDS(
  list(
    bench_df       = bench_df,
    best_trackA    = best_trackA,
    best_trackB    = best_trackB,
    qc_sens_summary= qc_sens_summary,
    pv_trackA      = pv_trackA,
    pv_trackB      = pv_trackB,
    pv_support_df  = pv_support_df,
    D_A_N10        = D_A,
    D_B_N10        = D_B,
    D_A_N6         = D_A_shared,
    D_B_N6         = D_B_shared,
    phylo_A_N10    = phylo_A_N10,
    phylo_B_N10    = phylo_B_N10,
    phylo_A_N6     = phylo_A_N6,
    phylo_B_N6     = phylo_B_N6,
    mantel_AB      = mantel_AB,
    bakers_AB      = bakers_AB,
    dist_qc        = dist_qc,
    dist_raw       = dist_raw
  ),
  file = file.path(DIRS$rds, "morpho_clustering_objects.rds")
)
message("[03] Saved clustering benchmark tables, pvclust bootstrap objects, and Newick trees (N=10 & N=6).")


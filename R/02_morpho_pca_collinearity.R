# ==============================================================================
# R/02_morpho_pca_collinearity.R
# Step 2: Parametric Assumption Checks (Shapiro-Wilk), Collinearity Audit (VIF),
# Pairwise Correlations (Spearman/Pearson + 95% CI + BH-FDR), and FactoMineR PCA
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(FactoMineR)
  library(factoextra)
})
reset_project_seed()

message("=== [02] Starting Normality, Collinearity (VIF) & PCA Analysis ===")

consensus_qc  <- read.csv(file.path(DIRS$processed, "morpho_consensus_qc.csv"),  stringsAsFactors = FALSE)
consensus_raw <- read.csv(file.path(DIRS$processed, "morpho_consensus_raw.csv"), stringsAsFactors = FALSE)

cont_cols <- c(
  "Fruit_Length_cm",
  "Fruit_Width_cm",
  "Stalk_Length_cm",
  "Stalk_Width_cm",
  "Diam_With_Skin_cm",
  "Diam_No_Skin_cm",
  "Skin_Thickness_cm"
)

pretty_trait_names <- c(
  "Fruit_Length_cm"   = "Fruit Length (cm)",
  "Fruit_Width_cm"    = "Fruit Width (cm)",
  "Stalk_Length_cm"   = "Stalk Length (cm)",
  "Stalk_Width_cm"    = "Stalk Width (cm)",
  "Diam_With_Skin_cm" = "Flesh Diam. w/ Skin (cm)",
  "Diam_No_Skin_cm"   = "Flesh Diam. w/o Skin (cm)",
  "Skin_Thickness_cm" = "Skin Thickness (cm)"
)

# 1. Shapiro-Wilk Normality Tests (AGENT.md Rule #2) ---------------------------
run_normality_audit <- function(df, dataset_label) {
  res <- lapply(cont_cols, function(col) {
    x <- df[[col]]
    sw <- shapiro.test(x)
    tibble(
      Dataset      = dataset_label,
      Trait        = col,
      Trait_Label  = pretty_trait_names[[col]],
      N            = length(x),
      Mean         = round(mean(x), 3),
      SD           = round(sd(x), 3),
      Median       = round(median(x), 3),
      IQR          = round(IQR(x), 3),
      Skewness     = round(mean(((x - mean(x)) / sd(x))^3), 3),
      Shapiro_W    = round(unname(sw$statistic), 4),
      p_value      = sw$p.value
    )
  }) %>% bind_rows()

  res %>%
    mutate(
      p_adj_BH       = p.adjust(p_value, method = "BH"),
      Normality_Pass = p_adj_BH >= 0.05
    )
}

norm_audit <- bind_rows(
  run_normality_audit(consensus_qc,  "QC_Corrected (N=10)"),
  run_normality_audit(consensus_raw, "Raw_Uncorrected (N=10)")
)
write.csv(norm_audit, file.path(DIRS$tables, "normality_shapiro_tests.csv"), row.names = FALSE)

# 2. Pairwise Spearman & Pearson Correlation Matrix with 95% CIs & BH-FDR ------
compute_pairwise_correlations <- function(df, dataset_label) {
  pairs <- combn(cont_cols, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pr) {
    x <- df[[pr[1]]]
    y <- df[[pr[2]]]
    ct_p <- suppressWarnings(cor.test(x, y, method = "pearson"))
    ct_s <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    # Fisher z 95% CI for Spearman rho
    r_s <- unname(ct_s$estimate)
    z_s <- atanh(pmin(pmax(r_s, -0.999), 0.999))
    se_z <- 1 / sqrt(max(length(x) - 3, 1))
    ci_s <- tanh(z_s + c(-1, 1) * 1.96 * se_z)

    tibble(
      Dataset          = dataset_label,
      Trait_1          = pr[1],
      Trait_2          = pr[2],
      Pearson_r        = round(unname(ct_p$estimate), 4),
      Pearson_CI_Low   = round(ct_p$conf.int[1], 4),
      Pearson_CI_High  = round(ct_p$conf.int[2], 4),
      Pearson_p        = ct_p$p.value,
      Spearman_rho     = round(r_s, 4),
      Spearman_CI_Low  = round(ci_s[1], 4),
      Spearman_CI_High = round(ci_s[2], 4),
      Spearman_p       = ct_s$p.value
    )
  }) %>% bind_rows()

  rows %>%
    mutate(
      Pearson_p_adj_BH  = p.adjust(Pearson_p, method = "BH"),
      Spearman_p_adj_BH = p.adjust(Spearman_p, method = "BH")
    )
}

corr_audit <- bind_rows(
  compute_pairwise_correlations(consensus_qc,  "QC_Corrected"),
  compute_pairwise_correlations(consensus_raw, "Raw_Uncorrected")
)
write.csv(corr_audit, file.path(DIRS$tables, "trait_correlations_pairwise.csv"), row.names = FALSE)

# 3. Variance Inflation Factor (VIF < 5) Collinearity Audit & Pruning ----------
compute_vif_vector <- function(mat_scaled) {
  R <- cor(mat_scaled, use = "complete.obs")
  # Ridge regularization only if near-singular
  vif_vals <- tryCatch(
    diag(solve(R)),
    error = function(e) diag(solve(R + diag(1e-6, ncol(R))))
  )
   vif_vals
}

X_qc_scaled  <- scale(as.matrix(consensus_qc[, cont_cols]))
rownames(X_qc_scaled) <- consensus_qc$Cultivar_Key

X_raw_scaled <- scale(as.matrix(consensus_raw[, cont_cols]))
rownames(X_raw_scaled) <- consensus_raw$Cultivar_Key

vif_initial <- compute_vif_vector(X_qc_scaled)

# Stepwise VIF pruning until all retained variables have VIF < 5.0
retained_cols <- cont_cols
step_id <- 1L
vif_history <- list()

repeat {
  curr_vifs <- compute_vif_vector(X_qc_scaled[, retained_cols, drop = FALSE])
  vif_history[[step_id]] <- tibble(
    Step      = step_id,
    Trait     = names(curr_vifs),
    VIF       = round(as.numeric(curr_vifs), 3),
    Max_VIF   = round(max(curr_vifs), 3),
    Status    = ifelse(curr_vifs < 5.0, "Pass (VIF < 5)", "Collinear (VIF >= 5)")
  )
  if (max(curr_vifs) < 5.0 || length(retained_cols) <= 3) break
  drop_var <- names(which.max(curr_vifs))
  retained_cols <- setdiff(retained_cols, drop_var)
  step_id <- step_id + 1L
}

vif_audit_df <- bind_rows(vif_history) %>%
  mutate(Retained_In_Final_VIF_Set = Trait %in% retained_cols)

write.csv(vif_audit_df, file.path(DIRS$tables, "vif_collinearity_audit.csv"), row.names = FALSE)
message(sprintf("[VIF Audit] Initial max VIF = %.2f; after pruning -> %d traits retained with max VIF = %.2f (%s)",
                max(vif_initial), length(retained_cols),
                max(compute_vif_vector(X_qc_scaled[, retained_cols, drop = FALSE])),
                paste(retained_cols, collapse = ", ")))

# 4. Principal Component Analysis (FactoMineR::PCA) ----------------------------
pca_qc_full <- FactoMineR::PCA(X_qc_scaled, scale.unit = FALSE, ncp = 5, graph = FALSE)
pca_qc_vif  <- FactoMineR::PCA(X_qc_scaled[, retained_cols, drop = FALSE], scale.unit = FALSE, ncp = 5, graph = FALSE)
pca_raw_full <- FactoMineR::PCA(X_raw_scaled, scale.unit = FALSE, ncp = 5, graph = FALSE)

extract_eig_table <- function(pca_obj, model_label) {
  eig <- as.data.frame(pca_obj$eig)
  tibble(
    PCA_Model              = model_label,
    Axis                   = rownames(eig),
    Eigenvalue             = round(eig[, 1], 4),
    Variance_Pct           = round(eig[, 2], 2),
    Cumulative_Variance_Pct= round(eig[, 3], 2),
    Kaiser_Criterion_GT1   = eig[, 1] > 1.0
  )
}

eig_table <- bind_rows(
  extract_eig_table(pca_qc_full,  "QC_Corrected_Full (7 traits)"),
  extract_eig_table(pca_qc_vif,   paste0("QC_Corrected_VIF_Pruned (", length(retained_cols), " traits)")),
  extract_eig_table(pca_raw_full, "Raw_Uncorrected_Full (7 traits)")
)
write.csv(eig_table, file.path(DIRS$tables, "pca_eigenvalues.csv"), row.names = FALSE)

# Extract loadings, cos2, and % contribution for QC Full and QC VIF models
extract_var_table <- function(pca_obj, model_label) {
  coord <- as.data.frame(pca_obj$var$coord[, 1:3, drop = FALSE])
  cos2  <- as.data.frame(pca_obj$var$cos2[, 1:3, drop = FALSE])
  contrib <- as.data.frame(pca_obj$var$contrib[, 1:3, drop = FALSE])
  tibble(
    PCA_Model       = model_label,
    Trait           = rownames(coord),
    Trait_Label     = unname(pretty_trait_names[rownames(coord)]),
    Loading_PC1     = round(coord[, 1], 4),
    Loading_PC2     = round(coord[, 2], 4),
    Loading_PC3     = round(coord[, 3], 4),
    Cos2_PC1        = round(cos2[, 1], 4),
    Cos2_PC2        = round(cos2[, 2], 4),
    Contrib_PC1_Pct = round(contrib[, 1], 2),
    Contrib_PC2_Pct = round(contrib[, 2], 2)
  )
}

var_table <- bind_rows(
  extract_var_table(pca_qc_full, "QC_Corrected_Full"),
  extract_var_table(pca_qc_vif,  "QC_Corrected_VIF_Pruned"),
  extract_var_table(pca_raw_full, "Raw_Uncorrected_Full")
)
write.csv(var_table, file.path(DIRS$tables, "pca_loadings_contributions.csv"), row.names = FALSE)

# Extract cultivar coordinates merged with metadata
extract_ind_table <- function(pca_obj, df_meta, model_label) {
  ind_coord <- as.data.frame(pca_obj$ind$coord[, 1:3, drop = FALSE])
  ind_cos2  <- as.data.frame(pca_obj$ind$cos2[, 1:3, drop = FALSE])
  ind_contrib <- as.data.frame(pca_obj$ind$contrib[, 1:3, drop = FALSE])
  tibble(
    PCA_Model           = model_label,
    Cultivar_Key        = rownames(ind_coord),
    PC1                 = round(ind_coord[, 1], 4),
    PC2                 = round(ind_coord[, 2], 4),
    PC3                 = round(ind_coord[, 3], 4),
    Cos2_PC1_PC2        = round(ind_cos2[, 1] + ind_cos2[, 2], 4),
    Contrib_PC1_Pct     = round(ind_contrib[, 1], 2),
    Contrib_PC2_Pct     = round(ind_contrib[, 2], 2)
  ) %>%
    left_join(df_meta %>% select(Cultivar_Key, Display_Name, Plot_Label_Markdown,
                                 Genome_Group, In_ITS, ITS_Accession_ID),
              by = "Cultivar_Key")
}

ind_table <- bind_rows(
  extract_ind_table(pca_qc_full,  consensus_qc,  "QC_Corrected_Full"),
  extract_ind_table(pca_qc_vif,   consensus_qc,  "QC_Corrected_VIF_Pruned"),
  extract_ind_table(pca_raw_full, consensus_raw, "Raw_Uncorrected_Full")
)
write.csv(ind_table, file.path(DIRS$tables, "pca_cultivar_coordinates.csv"), row.names = FALSE)

saveRDS(
  list(
    pca_qc_full   = pca_qc_full,
    pca_qc_vif    = pca_qc_vif,
    pca_raw_full  = pca_raw_full,
    X_qc_scaled   = X_qc_scaled,
    X_raw_scaled  = X_raw_scaled,
    retained_cols = retained_cols,
    cont_cols     = cont_cols,
    pretty_trait_names = pretty_trait_names,
    norm_audit    = norm_audit,
    corr_audit    = corr_audit,
    vif_audit_df  = vif_audit_df
  ),
  file = file.path(DIRS$rds, "morpho_pca_objects.rds")
)
message("[02] Completed Normality, Correlation, VIF, and PCA analyses.")


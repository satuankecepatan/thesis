# ==============================================================================
# R/10_concordance_mantel_procrustes.R
# Phase 3 - Step 2: Distance-Matrix Concordance (Exact 6! = 720 Permutation
# Mantel & Partial Mantel Tests), Single-Trait Evolutionary Signal Decomposition,
# Symmetric Procrustes Superimposition (M^2, r_P, Cultivar Residuals), and
# Quantitative Germplasm / Heterotic Breeding Pair Matrix
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(vegan)
  library(ape)
  library(cluster)
  library(permute)
})
reset_project_seed()

message("=== [10] Starting Mantel Permutation, Trait Decomposition & Procrustes Analysis ===")

morpho_cl  <- readRDS(file.path(DIRS$rds, "morpho_clustering_objects.rds"))
morpho_pca <- readRDS(file.path(DIRS$rds, "morpho_pca_objects.rds"))
phylo_objs <- readRDS(file.path(DIRS$rds, "molecular_phylogeny_objects.rds"))
topo_objs  <- readRDS(file.path(DIRS$rds, "concordance_topological_objects.rds"))
morpho_qc  <- read.csv(file.path(DIRS$processed, "morpho_consensus_qc.csv"), stringsAsFactors = FALSE)
meta_df    <- read.csv(file.path(DIRS$processed, "cultivar_metadata.csv"), stringsAsFactors = FALSE)

shared_meta <- topo_objs$shared_meta
shared_keys <- topo_objs$shared_keys
n_shared    <- length(shared_keys)

# 1. Assemble Synchronized N = 6 Distance Matrices -----------------------------
sub_mat <- function(D_obj, keys) {
  M <- as.matrix(D_obj)
  as.dist(M[keys, keys])
}

# Morphological Distance Matrices (N = 6)
D_morph_vif  <- sub_mat(morpho_cl$D_A_N6$by_cultivar_key, shared_keys)
D_morph_full <- sub_mat(morpho_cl$dist_qc[["Euclidean_Full (7 Cont)"]]$dist, shared_keys)
D_morph_gow  <- sub_mat(morpho_cl$D_B_N6$by_cultivar_key, shared_keys)
D_morph_coph <- as.dist(as.matrix(cophenetic(topo_objs$tr_morph_A))[shared_keys, shared_keys])

# Molecular ITS Distance Matrices (N = 6)
D_its_ml_phred  <- sub_mat(phylo_objs$D_its_ml_N6,  shared_keys)
D_its_pat_phred <- sub_mat(phylo_objs$D_its_pat_N6, shared_keys)

# Galaxy Raw ML & Patristic Distance Matrices (N = 6, mapped from ITS_Accession_ID to Cultivar_Key)
gal_ml_mat  <- as.matrix(phylo_objs$phylo_galaxy$d_ml_all11)[shared_meta$ITS_Accession_ID, shared_meta$ITS_Accession_ID]
gal_pat_mat <- as.matrix(phylo_objs$phylo_galaxy$d_patristic_all11)[shared_meta$ITS_Accession_ID, shared_meta$ITS_Accession_ID]
rownames(gal_ml_mat)  <- colnames(gal_ml_mat)  <- shared_keys
rownames(gal_pat_mat) <- colnames(gal_pat_mat) <- shared_keys
D_its_ml_galaxy  <- as.dist(gal_ml_mat[shared_keys, shared_keys])
D_its_pat_galaxy <- as.dist(gal_pat_mat[shared_keys, shared_keys])

# Genomic Ploidy / B-Genome Dosage Distance Matrix for Partial Mantel Control
# Quantifies standardized divergence in ploidy level (2x vs 3x) and B-genome fraction
genome_traits <- shared_meta %>%
  mutate(
    Ploidy_Num = case_when(
      Genome_Group == "AA"             ~ 2.0,
      Genome_Group == "AAA"            ~ 3.0,
      Genome_Group == "ABB"            ~ 3.0,
      Genome_Group == "AA/AAB (Local)" ~ 2.5,
      TRUE                             ~ 2.5
    ),
    B_Genome_Proportion = case_when(
      Genome_Group == "AA"             ~ 0.00,
      Genome_Group == "AAA"            ~ 0.00,
      Genome_Group == "AA/AAB (Local)" ~ 0.15,
      Genome_Group == "ABB"            ~ 0.67,
      TRUE                             ~ 0.00
    )
  )
G_mat <- scale(cbind(genome_traits$Ploidy_Num, genome_traits$B_Genome_Proportion))
rownames(G_mat) <- shared_keys
D_genome <- dist(G_mat, method = "euclidean")

# 2. Exact 6! = 720 Permutation Mantel Engine ----------------------------------
perm_mat <- rbind(seq_len(n_shared), permute::allPerms(n_shared))
n_perms  <- nrow(perm_mat) # 720 complete permutations of S_6

run_exact_mantel <- function(D_y, D_x, D_z = NULL, label_y = "", label_x = "") {
  M_y <- as.matrix(D_y)[shared_keys, shared_keys]
  M_x <- as.matrix(D_x)[shared_keys, shared_keys]
  v_y <- M_y[lower.tri(M_y)]
  v_x <- M_x[lower.tri(M_x)]

  if (!is.null(D_z)) {
    M_z <- as.matrix(D_z)[shared_keys, shared_keys]
    v_z <- M_z[lower.tri(M_z)]
  }

  calc_stat <- function(vy_curr, vx_curr, method) {
    ok <- !is.na(vy_curr) & !is.na(vx_curr)
    if (sum(ok) < 3) return(0.0)
    s_y <- sd(vy_curr[ok]); s_x <- sd(vx_curr[ok])
    if (is.na(s_y) || is.na(s_x) || s_y == 0 || s_x == 0) return(0.0)
    if (is.null(D_z)) {
      cor(vy_curr[ok], vx_curr[ok], method = method)
    } else {
      # Partial correlation controlling for v_z
      ryz <- cor(vy_curr[ok], v_z[ok], method = method)
      rxz <- cor(vx_curr[ok], v_z[ok], method = method)
      rxy <- cor(vy_curr[ok], vx_curr[ok], method = method)
      denom <- sqrt(max(1e-12, (1 - ryz^2) * (1 - rxz^2)))
      (rxy - ryz * rxz) / denom
    }
  }

  r_pear_obs <- calc_stat(v_y, v_x, "pearson")
  r_spear_obs <- calc_stat(v_y, v_x, "spearman")

  r_pear_null  <- numeric(n_perms)
  r_spear_null <- numeric(n_perms)

  for (p in seq_len(n_perms)) {
    idx <- perm_mat[p, ]
    M_x_perm <- M_x[idx, idx]
    v_x_perm <- M_x_perm[lower.tri(M_x_perm)]
    r_pear_null[p]  <- calc_stat(v_y, v_x_perm, "pearson")
    r_spear_null[p] <- calc_stat(v_y, v_x_perm, "spearman")
  }

  # Exact two-sided and one-sided permutation p-values
  p_pear_two  <- mean(abs(r_pear_null) >= (abs(r_pear_obs) - 1e-12))
  p_pear_pos  <- mean(r_pear_null >= (r_pear_obs - 1e-12))
  p_pear_neg  <- mean(r_pear_null <= (r_pear_obs + 1e-12))
  sd_p        <- sd(r_pear_null)
  z_pear      <- ifelse(sd_p > 0, (r_pear_obs - mean(r_pear_null)) / sd_p, 0.0)

  p_spear_two <- mean(abs(r_spear_null) >= (abs(r_spear_obs) - 1e-12))
  sd_s        <- sd(r_spear_null)
  z_spear     <- ifelse(sd_s > 0, (r_spear_obs - mean(r_spear_null)) / sd_s, 0.0)

  # Bootstrap 95% CI across pairwise distances
  reset_project_seed()
  boot_pear <- replicate(2000, {
    b_idx <- sample(seq_along(v_y), replace = TRUE)
    vy_b <- v_y[b_idx]; vx_b <- v_x[b_idx]
    ok_b <- !is.na(vy_b) & !is.na(vx_b)
    if (sum(ok_b) < 3) return(NA_real_)
    sy_b <- sd(vy_b[ok_b]); sx_b <- sd(vx_b[ok_b])
    if (is.na(sy_b) || is.na(sx_b) || sy_b == 0 || sx_b == 0) return(NA_real_)
    cor(vy_b[ok_b], vx_b[ok_b], method = "pearson")
  })
  ci_pear <- if (all(is.na(boot_pear))) c(0, 0) else quantile(boot_pear, probs = c(0.025, 0.975), na.rm = TRUE)

  tibble(
    Morphological_Distance   = label_y,
    Molecular_Distance       = label_x,
    Test_Type                = ifelse(is.null(D_z), "Standard Mantel", "Partial Mantel (| Genome Group)"),
    N_Taxa                   = n_shared,
    N_Pairs                  = length(v_y),
    Exact_Permutations       = n_perms,
    Mantel_r_Pearson         = round(r_pear_obs, 4),
    Pearson_CI_Low           = round(unname(ci_pear[1]), 4),
    Pearson_CI_High          = round(unname(ci_pear[2]), 4),
    Pearson_Z_Score          = round(z_pear, 3),
    Pearson_P_TwoSided       = round(p_pear_two, 4),
    Pearson_P_OneSided_Pos   = round(p_pear_pos, 4),
    Pearson_P_OneSided_Neg   = round(p_pear_neg, 4),
    Mantel_rho_Spearman      = round(r_spear_obs, 4),
    Spearman_Z_Score         = round(z_spear, 3),
    Spearman_P_TwoSided      = round(p_spear_two, 4)
  )
}

mantel_specs <- bind_rows(
  run_exact_mantel(D_morph_vif,  D_its_ml_phred,  NULL,     "Track A: Continuous VIF-Pruned Euclidean (6 Traits)", "Phred-Weighted ML Genetic Distance (TrN+G4)"),
  run_exact_mantel(D_morph_vif,  D_its_pat_phred, NULL,     "Track A: Continuous VIF-Pruned Euclidean (6 Traits)", "Phred-Weighted ML Tree Patristic Distance"),
  run_exact_mantel(D_morph_full, D_its_ml_phred,  NULL,     "Full Continuous Euclidean (7 Traits)",                "Phred-Weighted ML Genetic Distance (TrN+G4)"),
  run_exact_mantel(D_morph_gow,  D_its_ml_phred,  NULL,     "Track B: Mixed 12-Trait Gower Distance",              "Phred-Weighted ML Genetic Distance (TrN+G4)"),
  run_exact_mantel(D_morph_gow,  D_its_pat_phred, NULL,     "Track B: Mixed 12-Trait Gower Distance",              "Phred-Weighted ML Tree Patristic Distance"),
  run_exact_mantel(D_morph_coph, D_its_pat_phred, NULL,     "Track A: Morphological UPGMA Cophenetic Distance",    "Phred-Weighted ML Tree Patristic Distance"),
  run_exact_mantel(D_morph_vif,  D_its_ml_galaxy, NULL,     "Track A: Continuous VIF-Pruned Euclidean (6 Traits)", "Galaxy Unweighted ML Genetic Distance (TrN+G4)"),
  run_exact_mantel(D_morph_gow,  D_its_ml_galaxy, NULL,     "Track B: Mixed 12-Trait Gower Distance",              "Galaxy Unweighted ML Genetic Distance (TrN+G4)"),
  run_exact_mantel(D_morph_vif,  D_its_ml_phred,  D_genome, "Track A: Continuous VIF-Pruned Euclidean (6 Traits)", "Phred-Weighted ML Genetic Distance (TrN+G4)"),
  run_exact_mantel(D_morph_gow,  D_its_ml_phred,  D_genome, "Track B: Mixed 12-Trait Gower Distance",              "Phred-Weighted ML Genetic Distance (TrN+G4)")
) %>%
  mutate(
    Pearson_P_BH_Adj  = round(p.adjust(Pearson_P_TwoSided, method = "BH"), 4),
    Spearman_P_BH_Adj = round(p.adjust(Spearman_P_TwoSided, method = "BH"), 4),
    Biological_Regime = case_when(
      Pearson_P_TwoSided < 0.05 & Mantel_r_Pearson > 0 ~ "Positive Concordance (Lineage-Tracked Divergence)",
      Pearson_P_TwoSided < 0.05 & Mantel_r_Pearson < 0 ~ "Negative Concordance (Convergent / Somatic Decoupling)",
      TRUE                                             ~ "Null Concordance (Polyploid / Somatic Decoupling, r_M ~ 0)"
    )
  )

write.csv(mantel_specs, file.path(DIRS$tables, "concordance_mantel_tests.csv"), row.names = FALSE)

message(sprintf(
  "[Mantel Concordance] Primary Test (Track A vs Phred ML Dist): Pearson r_M = %.4f (95%% CI [%.3f, %.3f], exact two-sided p = %.4f); Spearman rho_M = %.4f (p = %.4f)",
  mantel_specs$Mantel_r_Pearson[1], mantel_specs$Pearson_CI_Low[1], mantel_specs$Pearson_CI_High[1],
  mantel_specs$Pearson_P_TwoSided[1], mantel_specs$Mantel_rho_Spearman[1], mantel_specs$Spearman_P_TwoSided[1]
))

# 3. Single-Trait Mantel Decomposition (12 Morphological Descriptors) ----------
morpho_shared_df <- morpho_qc %>%
  filter(Cultivar_Key %in% shared_keys) %>%
  arrange(match(Cultivar_Key, shared_keys))
rownames(morpho_shared_df) <- shared_keys

cont_traits <- c(
  "Fruit_Length_cm"   = "Fruit Length (cm)",
  "Fruit_Width_cm"    = "Fruit Width (cm)",
  "Stalk_Length_cm"   = "Pedicel / Stalk Length (cm)",
  "Stalk_Width_cm"    = "Pedicel / Stalk Width (cm)",
  "Diam_With_Skin_cm" = "Fruit Diameter With Peel (cm)",
  "Diam_No_Skin_cm"   = "Pulp Diameter Without Peel (cm)",
  "Skin_Thickness_cm" = "Pericarp / Peel Thickness (cm)"
)

cat_traits <- c(
  "Fruits_Per_Hand_Cat" = "Fruits Per Hand Class (IPGRI)",
  "Fruit_Position_Cat"  = "Fruit Rachis Posture (IPGRI)",
  "Fruit_Shape_Cat"     = "Longitudinal Fruit Curvature (IPGRI)",
  "Fruit_Tip_Cat"       = "Fruit Apex / Tip Architecture (IPGRI)",
  "Stalk_Surface_Cat"   = "Pedicel Surface Hairiness (IPGRI)",
  "Skin_Young_Cat"      = "Immature Pericarp Color (IPGRI)"
)

trait_rows <- list()

for (col_nm in names(cont_traits)) {
  vec <- scale(morpho_shared_df[[col_nm]])
  rownames(vec) <- shared_keys
  D_tr <- dist(vec, method = "euclidean")
  res_m <- run_exact_mantel(D_tr, D_its_ml_phred, NULL, cont_traits[[col_nm]], "Phred-Weighted ML Genetic Distance")
  res_m$Trait_Code <- col_nm
  res_m$Trait_Category <- ifelse(col_nm == "Diam_No_Skin_cm", "Continuous (Collinear, VIF >= 5)", "Continuous (Retained, VIF < 5)")
  trait_rows[[length(trait_rows) + 1]] <- res_m
}

for (col_nm in names(cat_traits)) {
  if (col_nm %in% colnames(morpho_shared_df) && length(unique(na.omit(morpho_shared_df[[col_nm]]))) > 1) {
    df_cat <- data.frame(val = factor(morpho_shared_df[[col_nm]]))
    rownames(df_cat) <- shared_keys
    D_tr <- cluster::daisy(df_cat, metric = "gower")
    res_m <- run_exact_mantel(D_tr, D_its_ml_phred, NULL, cat_traits[[col_nm]], "Phred-Weighted ML Genetic Distance")
    res_m$Trait_Code <- col_nm
    res_m$Trait_Category <- "Qualitative IPGRI Descriptor (Gower)"
    trait_rows[[length(trait_rows) + 1]] <- res_m
  }
}

trait_mantel_df <- bind_rows(trait_rows) %>%
  select(Trait_Code, Trait_Label = Morphological_Distance, Trait_Category,
         Mantel_r_Pearson, Pearson_CI_Low, Pearson_CI_High, Pearson_Z_Score, Pearson_P_TwoSided,
         Mantel_rho_Spearman, Spearman_Z_Score, Spearman_P_TwoSided) %>%
  mutate(
    Pearson_P_BH_Adj  = round(p.adjust(Pearson_P_TwoSided, method = "BH"), 4),
    Spearman_P_BH_Adj = round(p.adjust(Spearman_P_TwoSided, method = "BH"), 4)
  ) %>%
  arrange(desc(Mantel_r_Pearson))

write.csv(trait_mantel_df, file.path(DIRS$tables, "concordance_trait_mantel_decomposition.csv"), row.names = FALSE)

# 4. Symmetric Procrustes Superimposition & Cultivar Residual Audit ------------
# Principal Coordinates Analysis (PCoA, k = 2) on Morphological Track A and Phred ML Genetic Distance
pcoa_morph_A <- cmdscale(D_morph_vif, k = 2, eig = TRUE)
pcoa_morph_B <- cmdscale(D_morph_gow, k = 2, eig = TRUE)
pcoa_its_ml  <- cmdscale(D_its_ml_phred, k = 2, eig = TRUE)

X_its   <- pcoa_its_ml$points[shared_keys, 1:2]
Y_morph <- pcoa_morph_A$points[shared_keys, 1:2]

proc_fit <- vegan::procrustes(X = X_its, Y = Y_morph, symmetric = TRUE)

# Exact 6! = 720 permutation test for Symmetric Procrustes M^2 and r_P = sqrt(1 - M^2)
m2_obs <- proc_fit$ss
rp_obs <- sqrt(max(0, 1 - m2_obs))

m2_null <- numeric(n_perms)
for (p in seq_len(n_perms)) {
  idx <- perm_mat[p, ]
  pf_p <- vegan::procrustes(X = X_its, Y = Y_morph[idx, , drop = FALSE], symmetric = TRUE)
  m2_null[p] <- pf_p$ss
}
rp_null <- sqrt(pmax(0, 1 - m2_null))
p_proc_exact <- mean(m2_null <= (m2_obs + 1e-12))
z_proc       <- (rp_obs - mean(rp_null)) / sd(rp_null)

# Also evaluate Track B (Mixed Gower) Procrustes superimposition
Y_morph_B <- pcoa_morph_B$points[shared_keys, 1:2]
proc_fit_B <- vegan::procrustes(X = X_its, Y = Y_morph_B, symmetric = TRUE)
m2_obs_B   <- proc_fit_B$ss
rp_obs_B   <- sqrt(max(0, 1 - m2_obs_B))
m2_null_B  <- sapply(seq_len(n_perms), function(p) {
  vegan::procrustes(X = X_its, Y = Y_morph_B[perm_mat[p, ], , drop = FALSE], symmetric = TRUE)$ss
})
p_proc_B_exact <- mean(m2_null_B <= (m2_obs_B + 1e-12))

# Extract Cultivar-Specific Procrustes Coordinates and Residual Vectors (Track A Primary)
res_vec <- residuals(proc_fit)
res_sq  <- res_vec^2
pct_m2  <- 100 * res_sq / sum(res_sq)

q50_res <- median(res_vec)
q75_res <- quantile(res_vec, 0.75)

proc_cultivar_df <- tibble(
  Cultivar_Key        = shared_keys,
  Display_Name        = shared_meta$Display_Name,
  Binomial_Label      = shared_meta$Plot_Label_Plain,
  Genome_Group        = shared_meta$Genome_Group,
  ITS_Target_X1       = round(proc_fit$X[, 1], 4),
  ITS_Target_X2       = round(proc_fit$X[, 2], 4),
  Morph_Rotated_Y1    = round(proc_fit$Yrot[, 1], 4),
  Morph_Rotated_Y2    = round(proc_fit$Yrot[, 2], 4),
  Delta_Dim1          = round(proc_fit$Yrot[, 1] - proc_fit$X[, 1], 4),
  Delta_Dim2          = round(proc_fit$Yrot[, 2] - proc_fit$X[, 2], 4),
  Procrustes_Residual = round(res_vec, 4),
  Pct_Contribution_M2 = round(pct_m2, 2),
  Residual_Rank       = rank(-res_vec, ties.method = "first"),
  Concordance_Status  = case_when(
    res_vec >= q75_res ~ "High Discordance Outlier (>= Q75)",
    res_vec >= q50_res ~ "Moderate Divergence (Q50 - Q75)",
    TRUE               ~ "High Phenotypic-Genetic Concordance (< Q50)"
  ),
  Biological_Mechanism = case_when(
    Cultivar_Key == "Ambon"     ~ "Triploid Gigas Effect (AAA autopolyploid fruit enlargement relative to acuminata ITS clade)",
    Cultivar_Key == "Awak"      ~ "Allopolyploid B-Genome Divergence (ABB hybrid vigor and thick-skinned culinary fruit architecture)",
    Cultivar_Key == "Kalimbang" ~ "Local Cultivar Selection (Medium-fruited South Kalimantan landrace within acuminata ITS subclade)",
    Cultivar_Key == "Talas"     ~ "Divergent Pericarp Architecture (Thick-skinned, broad-diameter local AA cultivar vs small AA sisters)",
    Cultivar_Key == "Mas"       ~ "Conserved Diploid AA Archetype (Small, thin-skinned dessert fruit matching AA ITS position)",
    Cultivar_Key == "Mauli"     ~ "Conserved Diploid AA Archetype (Small, thin-skinned Banjar dessert cultivar matching AA ITS sisterhood)",
    TRUE                        ~ "Standard Phenotypic-Genetic Alignment"
  )
) %>%
  arrange(Residual_Rank)

write.csv(proc_cultivar_df, file.path(DIRS$tables, "concordance_procrustes_cultivar_residuals.csv"), row.names = FALSE)

message(sprintf(
  "[Procrustes Superimposition] Track A vs Phred ITS: M^2 = %.4f, r_P = %.4f (Z = %.2f, exact p = %.4f); Top Outlier = %s (eps = %.4f, %.1f%% of M^2)",
  m2_obs, rp_obs, z_proc, p_proc_exact,
  proc_cultivar_df$Cultivar_Key[1], proc_cultivar_df$Procrustes_Residual[1], proc_cultivar_df$Pct_Contribution_M2[1]
))

# 5. Pairwise Germplasm Concordance & Heterotic Breeding Matrix (15 Pairs) -----
M_morph_A <- as.matrix(D_morph_vif)[shared_keys, shared_keys]
M_morph_B <- as.matrix(D_morph_gow)[shared_keys, shared_keys]
M_its_ml  <- as.matrix(D_its_ml_phred)[shared_keys, shared_keys]
M_its_pat <- as.matrix(D_its_pat_phred)[shared_keys, shared_keys]

pair_comb <- combn(shared_keys, 2, simplify = FALSE)

pair_rows <- lapply(pair_comb, function(pr) {
  k1 <- pr[1]; k2 <- pr[2]
  g1 <- shared_meta$Genome_Group[shared_meta$Cultivar_Key == k1]
  g2 <- shared_meta$Genome_Group[shared_meta$Cultivar_Key == k2]
  g1_short <- sub(" \\(Local\\)", "", g1)
  g2_short <- sub(" \\(Local\\)", "", g2)

  has_abb <- ("ABB" %in% c(g1_short, g2_short))
  has_aaa <- ("AAA" %in% c(g1_short, g2_short))

  cross_type <- case_when(
    has_abb ~ "Interspecific (Acuminata AA/AAA x Balbisiana ABB)",
    has_aaa ~ "Intra-Acuminata Ploidy Bridge (Diploid AA x Triploid AAA)",
    TRUE    ~ "Intra-Acuminata Diploid Cross (AA x AA/AAB)"
  )

  d_m <- M_morph_A[k1, k2]
  d_g <- M_its_ml[k1, k2]

  tibble(
    Parent_1                = k1,
    Parent_2                = k2,
    Pair_Label              = paste0(k1, " \u00d7 ", k2),
    Pair_Label_ASCII        = paste0(k1, " x ", k2),
    Genome_Cross            = paste0(g1_short, " x ", g2_short),
    Cross_Category          = cross_type,
    ITS_ML_Dist             = round(d_g, 5),
    ITS_Patristic_Dist      = round(M_its_pat[k1, k2], 5),
    Morph_Euclidean_VIF     = round(d_m, 4),
    Morph_Gower_12Trait     = round(M_morph_B[k1, k2], 4)
  )
})

pairwise_df <- bind_rows(pair_rows) %>%
  mutate(
    Z_ITS_Dist          = as.numeric(scale(ITS_ML_Dist)),
    Z_Morph_Dist        = as.numeric(scale(Morph_Euclidean_VIF)),
    Discordance_Diff_Z  = round(Z_Morph_Dist - Z_ITS_Dist, 3),
    Heterotic_Index     = round(0.5 * (Z_ITS_Dist + Z_Morph_Dist), 3),
    Breeding_Classification = case_when(
      Z_ITS_Dist >= 0.4 & Z_Morph_Dist >= 0.4 ~ "Transgressive Heterotic Cross (High Genetic & High Phenotypic Divergence)",
      Z_ITS_Dist >= 0.4 & Z_Morph_Dist <  0.4 ~ "Cryptic Genetic Diversity Cross (High ITS Divergence, Convergent Fruit Type)",
      Z_ITS_Dist <  0.4 & Z_Morph_Dist >= 0.4 ~ "Polyploid / Somatic Phenotypic Shift (Close ITS Lineage, Divergent Fruit Architecture)",
      TRUE                                    ~ "Conserved Within-Clade Cross (Close Genetic & Phenotypic Similarity)"
    )
  ) %>%
  arrange(desc(Heterotic_Index))

write.csv(pairwise_df, file.path(DIRS$tables, "concordance_pairwise_breeding_matrix.csv"), row.names = FALSE)

saveRDS(
  list(
    mantel_specs     = mantel_specs,
    trait_mantel_df  = trait_mantel_df,
    proc_fit         = proc_fit,
    proc_fit_B       = proc_fit_B,
    m2_obs           = m2_obs,
    rp_obs           = rp_obs,
    p_proc_exact     = p_proc_exact,
    z_proc           = z_proc,
    m2_null          = m2_null,
    rp_null          = rp_null,
    m2_obs_B         = m2_obs_B,
    rp_obs_B         = rp_obs_B,
    p_proc_B_exact   = p_proc_B_exact,
    proc_cultivar_df = proc_cultivar_df,
    pairwise_df      = pairwise_df,
    pcoa_morph_A     = pcoa_morph_A,
    pcoa_its_ml      = pcoa_its_ml
  ),
  file = file.path(DIRS$rds, "concordance_analysis_objects.rds")
)
message("[10] Saved Mantel tests, trait decomposition, Procrustes residuals, and pairwise breeding matrix.")


# ==============================================================================
# R/09_concordance_topological_metrics.R
# Phase 3 - Step 1: Cophylogenetic & Topological Concordance Analysis
# (Step2Side Untangling, Entanglement E, Baker's Gamma Rank Correlation with
# Exact 6! = 720 Complete Permutation Null Distribution, Robinson-Foulds nRF,
# Steel-Penny Path Distance, and Kuhner-Felsenstein Branch Score Distance)
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(dendextend)
  library(permute)
})
reset_project_seed()

message("=== [09] Starting Topological & Cophylogenetic Concordance Analysis ===")

morpho_cl  <- readRDS(file.path(DIRS$rds, "morpho_clustering_objects.rds"))
phylo_objs <- readRDS(file.path(DIRS$rds, "molecular_phylogeny_objects.rds"))
meta_df    <- read.csv(file.path(DIRS$processed, "cultivar_metadata.csv"), stringsAsFactors = FALSE)

shared_meta <- meta_df %>%
  filter(In_Morpho & In_ITS) %>%
  arrange(Cultivar_Key)

shared_keys <- shared_meta$Cultivar_Key
stopifnot(length(shared_keys) == 6L)

# 1. Construct Synchronized Ultrametric Dendrograms & Phylo Objects (N = 6) ----
# Morphological Track A (VIF-Pruned 6 Continuous Traits, Euclidean + UPGMA)
D_morph_A <- as.dist(as.matrix(morpho_cl$D_A_N6$by_cultivar_key)[shared_keys, shared_keys])
hc_morph_A <- hclust(D_morph_A, method = "average")
dend_morph_A <- as.dendrogram(hc_morph_A)
tr_morph_A   <- ape::as.phylo(hc_morph_A)

# Morphological Track B (Mixed 12-Trait Gower + UPGMA)
D_morph_B <- as.dist(as.matrix(morpho_cl$D_B_N6$by_cultivar_key)[shared_keys, shared_keys])
hc_morph_B <- hclust(D_morph_B, method = "average")
dend_morph_B <- as.dendrogram(hc_morph_B)
tr_morph_B   <- ape::as.phylo(hc_morph_B)

# Helper to convert a rooted ML phylogram / patristic distance matrix (N = 6)
# into a valid ultrametric dendrogram with root height normalized to 1.0
phylo_to_ultrametric_dend <- function(tr_in, target_keys, D_fallback = NULL) {
  stopifnot(setequal(tr_in$tip.label, target_keys))
  D_pat <- if (!is.null(D_fallback)) {
    as.dist(as.matrix(D_fallback)[target_keys, target_keys])
  } else {
    as.dist(ape::cophenetic.phylo(tr_in)[target_keys, target_keys])
  }
  hc_u <- hclust(D_pat, method = "average")
  hc_u$height <- hc_u$height / max(hc_u$height)
  as.dendrogram(hc_u)
}

# Normalize morphological hclust heights to [0, 1] before converting to dendrogram
hc_morph_A$height <- hc_morph_A$height / max(hc_morph_A$height)
dend_morph_A <- as.dendrogram(hc_morph_A)

hc_morph_B$height <- hc_morph_B$height / max(hc_morph_B$height)
dend_morph_B <- as.dendrogram(hc_morph_B)

# Phred-Weighted ML Tree (N = 6, synced to Cultivar_Key)
tr_its_phred <- phylo_objs$tr_N6_synced
dend_its_phred <- phylo_to_ultrametric_dend(tr_its_phred, shared_keys, phylo_objs$D_its_pat_N6)

# Galaxy Raw ML Tree (N = 6, pruned from Galaxy N = 11 and synced to Cultivar_Key)
tr_gal_full <- phylo_objs$phylo_galaxy$tr_rooted
tr_its_galaxy <- ape::keep.tip(tr_gal_full, shared_meta$ITS_Accession_ID)
m_gal <- match(tr_its_galaxy$tip.label, shared_meta$ITS_Accession_ID)
tr_its_galaxy$tip.label <- shared_meta$Cultivar_Key[m_gal]
dend_its_galaxy <- phylo_to_ultrametric_dend(tr_its_galaxy, shared_keys, NULL)

# 2. Evaluate All 4 Morphological x Molecular Tree Pairings --------------------
pair_specs <- list(
  list(
    Pair_ID      = "TrackA_vs_Phred",
    Morph_Model  = "Track A: Continuous VIF-Pruned (Euclidean + UPGMA)",
    Mol_Model    = "Phred-Weighted Consensus (TrN+G4 ML)",
    dend_morph   = dend_morph_A,
    dend_mol     = dend_its_phred,
    tr_morph     = tr_morph_A,
    tr_mol       = tr_its_phred
  ),
  list(
    Pair_ID      = "TrackB_vs_Phred",
    Morph_Model  = "Track B: Mixed 12-Trait (Gower + UPGMA)",
    Mol_Model    = "Phred-Weighted Consensus (TrN+G4 ML)",
    dend_morph   = dend_morph_B,
    dend_mol     = dend_its_phred,
    tr_morph     = tr_morph_B,
    tr_mol       = tr_its_phred
  ),
  list(
    Pair_ID      = "TrackA_vs_Galaxy",
    Morph_Model  = "Track A: Continuous VIF-Pruned (Euclidean + UPGMA)",
    Mol_Model    = "Galaxy Unweighted Consensus (TrN+G4 ML)",
    dend_morph   = dend_morph_A,
    dend_mol     = dend_its_galaxy,
    tr_morph     = tr_morph_A,
    tr_mol       = tr_its_galaxy
  ),
  list(
    Pair_ID      = "TrackB_vs_Galaxy",
    Morph_Model  = "Track B: Mixed 12-Trait (Gower + UPGMA)",
    Mol_Model    = "Galaxy Unweighted Consensus (TrN+G4 ML)",
    dend_morph   = dend_morph_B,
    dend_mol     = dend_its_galaxy,
    tr_morph     = tr_morph_B,
    tr_mol       = tr_its_galaxy
  )
)

# Generate all 6! = 720 exact permutations of the 6 cultivar labels for exact null testing
perm_mat <- rbind(seq_len(6), permute::allPerms(6))
n_exact_perms <- nrow(perm_mat) # 720 exact permutations

topo_rows <- list()
null_rows <- list()
untangled_pairs <- list()

for (i in seq_along(pair_specs)) {
  sp <- pair_specs[[i]]
  dl_init <- dendextend::dendlist(sp$dend_morph, sp$dend_mol)
  ent_pre <- dendextend::entanglement(dl_init)

  # Benchmark untangling algorithms to guarantee global minimum entanglement
  dl_s1 <- dendextend::untangle(dl_init, method = "step1side")
  dl_s2 <- dendextend::untangle(dl_init, method = "step2side")
  dl_lb <- dendextend::untangle(dl_init, method = "labels")

  cands <- list(step2side = dl_s2, step1side = dl_s1, labels = dl_lb, initial = dl_init)
  ents  <- sapply(cands, dendextend::entanglement)
  best_idx <- which.min(ents)
  dl_opt   <- cands[[best_idx]]
  ent_post <- ents[best_idx]

  untangled_pairs[[sp$Pair_ID]] <- list(
    dend_list = dl_opt,
    ent_pre   = ent_pre,
    ent_post  = ent_post,
    method    = names(cands)[best_idx]
  )

  # Observed Baker's Gamma rank correlation
  bg_obs <- dendextend::cor_bakers_gamma(dl_opt[[1]], dl_opt[[2]])

  # Exact 6! = 720 permutation null distribution for Baker's Gamma
  orig_labels <- labels(dl_opt[[2]])
  bg_null <- numeric(n_exact_perms)
  for (p_idx in seq_len(n_exact_perms)) {
    d2_perm <- dl_opt[[2]]
    labels(d2_perm) <- orig_labels[perm_mat[p_idx, ]]
    bg_null[p_idx]  <- dendextend::cor_bakers_gamma(dl_opt[[1]], d2_perm)
  }

  # Exact two-sided permutation p-value and standardized Z-score
  p_exact <- mean(abs(bg_null) >= (abs(bg_obs) - 1e-12))
  bg_sd   <- sd(bg_null)
  bg_z    <- ifelse(bg_sd > 0, (bg_obs - mean(bg_null)) / bg_sd, 0.0)

  # Bootstrap 95% CI for Baker's Gamma via pairwise cophenetic rank resampling
  coph_1 <- as.numeric(as.matrix(cophenetic(dl_opt[[1]]))[shared_keys, shared_keys][lower.tri(matrix(0, 6, 6))])
  coph_2 <- as.numeric(as.matrix(cophenetic(dl_opt[[2]]))[shared_keys, shared_keys][lower.tri(matrix(0, 6, 6))])
  reset_project_seed()
  bg_boot <- replicate(2000, {
    idx_b <- sample(seq_along(coph_1), replace = TRUE)
    if (sd(coph_1[idx_b]) == 0 || sd(coph_2[idx_b]) == 0) return(NA_real_)
    cor(coph_1[idx_b], coph_2[idx_b], method = "spearman")
  })
  bg_ci <- quantile(bg_boot, probs = c(0.025, 0.975), na.rm = TRUE)

  # Topological & Branch-Length Distances (phangorn)
  # Use unrooted trees for standard Robinson-Foulds bipartition comparison (max splits = 2*(6-3) = 6)
  tr1_u <- ape::unroot(ape::multi2di(sp$tr_morph))
  tr2_u <- ape::unroot(ape::multi2di(sp$tr_mol))
  rf_raw  <- as.numeric(phangorn::RF.dist(tr1_u, tr2_u, normalize = FALSE, check.labels = TRUE))
  rf_norm <- as.numeric(phangorn::RF.dist(tr1_u, tr2_u, normalize = TRUE,  check.labels = TRUE))
  max_splits <- 2 * (length(shared_keys) - 3)
  shared_splits <- max(0, (max_splits - rf_raw) / 2)

  # Rooted RF distance
  rf_rooted_norm <- as.numeric(phangorn::RF.dist(
    ape::multi2di(sp$tr_morph), ape::multi2di(sp$tr_mol),
    rooted = TRUE, normalize = TRUE, check.labels = TRUE
  ))

  # Steel-Penny Path Distance & Kuhner-Felsenstein Branch-Score Distance
  path_d <- as.numeric(phangorn::path.dist(tr1_u, tr2_u, check.labels = TRUE))
  kf_d   <- as.numeric(phangorn::KF.dist(tr1_u, tr2_u, check.labels = TRUE))

  topo_rows[[i]] <- tibble(
    Pair_ID                  = sp$Pair_ID,
    Morphological_Tree       = sp$Morph_Model,
    Molecular_Tree           = sp$Mol_Model,
    N_Shared_Taxa            = length(shared_keys),
    Entanglement_Pre         = round(ent_pre, 4),
    Entanglement_Step2Side   = round(ent_post, 4),
    Bakers_Gamma             = round(bg_obs, 4),
    Bakers_Gamma_CI_Low      = round(unname(bg_ci[1]), 4),
    Bakers_Gamma_CI_High     = round(unname(bg_ci[2]), 4),
    Bakers_Null_Mean         = round(mean(bg_null), 4),
    Bakers_Null_SD           = round(bg_sd, 4),
    Bakers_Z_Score           = round(bg_z, 3),
    Bakers_Exact_P_TwoSided  = round(p_exact, 4),
    RF_Unrooted_Raw          = rf_raw,
    RF_Unrooted_Normalized   = round(rf_norm, 4),
    Shared_Internal_Splits   = shared_splits,
    RF_Rooted_Normalized     = round(rf_rooted_norm, 4),
    Steel_Penny_Path_Dist    = round(path_d, 4),
    Kuhner_Felsenstein_Dist  = round(kf_d, 4)
  )

  null_rows[[i]] <- tibble(
    Pair_ID        = sp$Pair_ID,
    Perm_Index     = seq_len(n_exact_perms),
    Bakers_Gamma   = bg_null,
    Observed_Gamma = bg_obs
  )
}

topo_df <- bind_rows(topo_rows) %>%
  mutate(Bakers_P_BH_Adj = round(p.adjust(Bakers_Exact_P_TwoSided, method = "BH"), 4))

null_df <- bind_rows(null_rows)

write.csv(topo_df, file.path(DIRS$tables, "concordance_topological_metrics.csv"), row.names = FALSE)
write.csv(null_df, file.path(DIRS$tables, "bakers_gamma_null_distribution.csv"), row.names = FALSE)

message(sprintf(
  "[Topological Concordance] Primary Pair (Track A vs Phred-Weighted): Entanglement E = %.3f (pre: %.3f), Baker's Gamma = %.4f (Z = %.2f, exact p = %.4f), Normalized RF = %.3f",
  topo_df$Entanglement_Step2Side[1], topo_df$Entanglement_Pre[1],
  topo_df$Bakers_Gamma[1], topo_df$Bakers_Z_Score[1],
  topo_df$Bakers_Exact_P_TwoSided[1], topo_df$RF_Unrooted_Normalized[1]
))

saveRDS(
  list(
    topo_df         = topo_df,
    null_df         = null_df,
    untangled_pairs = untangled_pairs,
    shared_meta     = shared_meta,
    shared_keys     = shared_keys,
    dend_morph_A    = dend_morph_A,
    dend_morph_B    = dend_morph_B,
    dend_its_phred  = dend_its_phred,
    dend_its_galaxy = dend_its_galaxy,
    tr_morph_A      = tr_morph_A,
    tr_morph_B      = tr_morph_B,
    tr_its_phred    = tr_its_phred,
    tr_its_galaxy   = tr_its_galaxy
  ),
  file = file.path(DIRS$rds, "concordance_topological_objects.rds")
)
message("[09] Saved concordance_topological_metrics.csv, bakers_gamma_null_distribution.csv, and RDS objects.")


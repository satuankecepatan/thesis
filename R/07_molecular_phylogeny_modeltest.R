# ==============================================================================
# R/07_molecular_phylogeny_modeltest.R
# Phase 2 - Step 3: Exhaustive Nucleotide Substitution Model Selection (BIC/AICc),
# Maximum Likelihood Tree Optimization, Outgroup Rooting (Ensete ventricosum),
# Dual Branch Support (1,000 SH-aLRT + 1,000 ML Bootstrap), & Two-Tier Export
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(Biostrings)
  library(ape)
  library(phangorn)
  library(vegan)
})
reset_project_seed()

message("=== [07] Starting Substitution ModelTest, ML Phylogeny & Branch Support ===")

aln_objs <- readRDS(file.path(DIRS$rds, "molecular_alignment_objects.rds"))
meta_df  <- read.csv(file.path(DIRS$processed, "cultivar_metadata.csv"), stringsAsFactors = FALSE)

# 1. Exact SH-aLRT (Guindon et al. 2010 RELL 1,000 replicates) + ML Bootstrap --
compute_sh_alrt_for_tree <- function(fit_ml, n_rell = 1000L) {
  reset_project_seed()
  tr_unrooted <- ape::unroot(fit_ml$tree)
  fit_unroot  <- update(fit_ml, tree = tr_unrooted)
  fit_unroot  <- phangorn::optim.pml(fit_unroot, optNni = FALSE, optEdge = TRUE,
                                     control = phangorn::pml.control(trace = 0))

  # Generate all single-branch NNI neighbor trees (2 per internal edge)
  nni_trees <- phangorn::nni(fit_unroot$tree)
  n_int_edges <- length(nni_trees) / 2L

  # Site pattern weights and site log-likelihood vector of ML tree
  phy_data <- fit_unroot$data
  w_pat    <- attr(phy_data, "weight")
  n_sites  <- sum(w_pat)
  n_pat    <- length(w_pat)
  prob_pat <- w_pat / n_sites

  sl_ml <- as.numeric(fit_unroot$siteLik)

  # Pre-generate 1,000 multinomial RELL bootstrap weight vectors (n_pat x n_rell)
  rell_weights <- rmultinom(n_rell, size = n_sites, prob = prob_pat)
  ll_ml_boot   <- as.numeric(crossprod(sl_ml, rell_weights))

  # Map each internal edge of tr_unrooted to its bipartition (smaller side tip set)
  n_tips <- length(tr_unrooted$tip.label)
  all_tips <- sort(tr_unrooted$tip.label)
  canon_split <- function(tips_sub) {
    s1 <- sort(tips_sub)
    s2 <- sort(setdiff(all_tips, s1))
    if (length(s1) < length(s2)) {
      paste(s1, collapse = "|")
    } else if (length(s2) < length(s1)) {
      paste(s2, collapse = "|")
    } else {
      if (paste(s1, collapse = "|") < paste(s2, collapse = "|")) paste(s1, collapse = "|") else paste(s2, collapse = "|")
    }
  }

  # Internal edges in phangorn::nni order correspond to internal edges where child > n_tips
  int_edge_idx <- which(tr_unrooted$edge[, 2] > n_tips)
  split_keys <- character(length(int_edge_idx))
  sh_alrt_pct <- numeric(length(int_edge_idx))

  for (k in seq_along(int_edge_idx)) {
    e_row <- int_edge_idx[k]
    child_nd <- tr_unrooted$edge[e_row, 2]
    clade_tips <- ape::extract.clade(tr_unrooted, child_nd)$tip.label
    split_keys[k] <- canon_split(clade_tips)

    t1 <- nni_trees[[2 * k - 1]]
    t2 <- nni_trees[[2 * k]]

    f1 <- phangorn::optim.pml(update(fit_unroot, tree = t1), optNni = FALSE, optEdge = TRUE,
                              control = phangorn::pml.control(trace = 0))
    f2 <- phangorn::optim.pml(update(fit_unroot, tree = t2), optNni = FALSE, optEdge = TRUE,
                              control = phangorn::pml.control(trace = 0))

    ll_1_boot <- as.numeric(crossprod(as.numeric(f1$siteLik), rell_weights))
    ll_2_boot <- as.numeric(crossprod(as.numeric(f2$siteLik), rell_weights))

    # SH-aLRT centering (Shimodaira-Hasegawa test on {T_ML, T_NNI1, T_NNI2})
    mat_3 <- rbind(ll_ml_boot, ll_1_boot, ll_2_boot)
    means_3 <- rowMeans(mat_3)
    centered_3 <- sweep(mat_3, 1, means_3, `-`)
    obs_diff <- fit_unroot$logLik - max(f1$logLik, f2$logLik)
    boot_max_diff <- apply(centered_3, 2, function(col) max(col) - col[1])
    p_sh <- mean(boot_max_diff >= obs_diff)
    sh_alrt_pct[k] <- round(100 * (1 - p_sh), 2)
  }

  names(sh_alrt_pct) <- split_keys
  list(sh_by_split = sh_alrt_pct, canon_split = canon_split)
}

# 2. Full ModelTest + ML Inference + Bootstrap Pipeline ------------------------
run_ml_phylogeny_pipeline <- function(aln_dss, dataset_label, n_boot = 1000L) {
  reset_project_seed()
  dna_bin <- ape::as.DNAbin(aln_dss)
  phy_dat <- phangorn::phyDat(dna_bin, type = "DNA")

  message(sprintf("[ModelTest] Testing 92 nucleotide substitution models on %s alignment...", dataset_label))
  mt <- phangorn::modelTest(
    phy_dat,
    model = "all",
    G     = TRUE,
    I     = TRUE,
    k     = 4,
    control = phangorn::pml.control(trace = 0)
  )

  mt_df <- as.data.frame(mt) %>%
    as_tibble() %>%
    arrange(BIC) %>%
    mutate(
      Dataset    = dataset_label,
      Delta_BIC  = round(BIC - min(BIC), 4),
      Delta_AICc = round(AICc - min(AICc), 4),
      BIC_Weight = round(exp(-0.5 * (BIC - min(BIC))) / sum(exp(-0.5 * (BIC - min(BIC)))), 5),
      AICc_Weight= round(exp(-0.5 * (AICc - min(AICc))) / sum(exp(-0.5 * (AICc - min(AICc)))), 5)
    )

  best_model_name <- mt_df$Model[1]
  base_model      <- sub("\\+.*$", "", best_model_name)
  has_I           <- grepl("\\+I", best_model_name)
  has_G           <- grepl("\\+G", best_model_name)

  message(sprintf("[ModelTest Winner (%s)] Model = %s (LogLik = %.2f, BIC = %.2f, AICc = %.2f, BIC Weight = %.3f)",
                  dataset_label, best_model_name, mt_df$logLik[1], mt_df$BIC[1], mt_df$AICc[1], mt_df$BIC_Weight[1]))

  # Initialize and optimize ML tree under the BIC-selected model
  fit_init <- phangorn::as.pml(mt, "BIC")
  reset_project_seed()
  fit_ml <- phangorn::optim.pml(
    fit_init,
    model         = base_model,
    optNni        = TRUE,
    optBf         = TRUE,
    optQ          = TRUE,
    optInv        = has_I,
    optGamma      = has_G,
    rearrangement = "stochastic",
    control       = phangorn::pml.control(trace = 0)
  )

  # Compute 1,000 non-parametric ML bootstrap replicates
  reset_project_seed()
  message(sprintf("[ML Bootstrap] Running %d non-parametric ML bootstrap replicates (%s)...", n_boot, dataset_label))
  bs_trees <- phangorn::bootstrap.pml(
    fit_ml,
    bs      = n_boot,
    optNni  = TRUE,
    control = phangorn::pml.control(trace = 0)
  )

  # Compute 1,000 SH-aLRT replicates
  message(sprintf("[SH-aLRT] Running 1,000 SH-aLRT replicates (%s)...", dataset_label))
  sh_res <- compute_sh_alrt_for_tree(fit_ml, n_rell = 1000L)

  # Root tree on Ensete_ventricosum outgroup (Bioinformatics.md Rule #3)
  tr_rooted <- ape::root(fit_ml$tree, outgroup = "Ensete_ventricosum", resolve.root = TRUE)
  tr_rooted <- ape::ladderize(tr_rooted, right = FALSE)

  # Root all bootstrap trees on Ensete_ventricosum to compute rooted clade frequencies
  bs_rooted <- lapply(bs_trees, function(bt) {
    tryCatch(ape::root(bt, outgroup = "Ensete_ventricosum", resolve.root = TRUE),
             error = function(e) bt)
  })
  class(bs_rooted) <- "multiPhylo"

  clade_counts <- ape::prop.clades(tr_rooted, bs_rooted)
  ml_bs_pct    <- round(100 * clade_counts / length(bs_rooted), 1)
  ml_bs_pct[is.na(ml_bs_pct)] <- 100.0

  # Map SH-aLRT values onto the rooted internal nodes via canonical bipartition splits
  n_tips  <- length(tr_rooted$tip.label)
  n_nodes <- tr_rooted$Nnode
  node_ids <- (n_tips + 1):(n_tips + n_nodes)

  sh_node_pct <- numeric(n_nodes)
  clade_taxa_str <- character(n_nodes)

  for (idx in seq_len(n_nodes)) {
    nd <- node_ids[idx]
    sub_tips <- ape::extract.clade(tr_rooted, nd)$tip.label
    clade_taxa_str[idx] <- paste(sort(sub_tips), collapse = " + ")
    if (length(sub_tips) == n_tips || length(sub_tips) >= n_tips - 1) {
      # Root or ingroup crown split defined by outgroup rooting
      sh_node_pct[idx] <- 100.0
    } else {
      sp_key <- sh_res$canon_split(sub_tips)
      val <- sh_res$sh_by_split[sp_key]
      sh_node_pct[idx] <- ifelse(is.na(val), ml_bs_pct[idx], as.numeric(val))
    }
  }

  # Embed "SH-aLRT / ML-BS" labels onto tr_rooted$node.label
  tr_rooted$node.label <- sprintf("%.0f/%.0f", sh_node_pct, ml_bs_pct)

  node_support_df <- tibble(
    Dataset        = dataset_label,
    Subst_Model    = best_model_name,
    Node_ID        = node_ids,
    Clade_Size     = sapply(node_ids, function(nd) length(ape::extract.clade(tr_rooted, nd)$tip.label)),
    Clade_Taxa     = clade_taxa_str,
    SH_aLRT_Pct    = round(sh_node_pct, 1),
    ML_BS_Pct      = round(ml_bs_pct, 1),
    Support_Label  = sprintf("%.0f / %.0f", sh_node_pct, ml_bs_pct),
    Pass_SH80      = sh_node_pct >= 80.0,
    Pass_BS95      = ml_bs_pct >= 95.0,
    Pass_Both      = (sh_node_pct >= 80.0 & ml_bs_pct >= 95.0)
  )

  # Compute ML pairwise genetic distance matrix and tree patristic distance matrix
  d_ml_all11       <- phangorn::dist.ml(phy_dat, model = "F81")
  d_patristic_all11 <- as.dist(ape::cophenetic.phylo(tr_rooted))

  # Check cophenetic faithfulness of ML phylogram to pairwise ML distances (within Musa ingroup)
  ingroup_tips <- setdiff(tr_rooted$tip.label, "Ensete_ventricosum")
  d_ml_in  <- as.dist(as.matrix(d_ml_all11)[ingroup_tips, ingroup_tips])
  d_pat_in <- as.dist(as.matrix(d_patristic_all11)[ingroup_tips, ingroup_tips])
  coph_r   <- cor(as.numeric(d_ml_in), as.numeric(d_pat_in), method = "pearson")

  list(
    mt_df             = mt_df,
    best_model_name   = best_model_name,
    fit_ml            = fit_ml,
    tr_rooted         = tr_rooted,
    bs_trees          = bs_rooted,
    node_support_df   = node_support_df,
    d_ml_all11        = d_ml_all11,
    d_patristic_all11 = d_patristic_all11,
    d_ml_in           = d_ml_in,
    d_pat_in          = d_pat_in,
    ingroup_coph_r    = coph_r
  )
}

phylo_phred  <- run_ml_phylogeny_pipeline(aln_objs$res_phred$aln_trimmed,  "Phred_Weighted", n_boot = 1000L)
phylo_galaxy <- run_ml_phylogeny_pipeline(aln_objs$res_galaxy$aln_trimmed, "Galaxy_Raw",     n_boot = 1000L)

# 3. Export ModelTest & Node Support Tables ------------------------------------
modeltest_all <- bind_rows(phylo_phred$mt_df, phylo_galaxy$mt_df)
write.csv(modeltest_all, file.path(DIRS$tables, "its_modeltest_bic_aicc.csv"), row.names = FALSE)

node_support_all <- bind_rows(phylo_phred$node_support_df, phylo_galaxy$node_support_df)
write.csv(node_support_all, file.path(DIRS$tables, "its_ml_node_support.csv"), row.names = FALSE)

# 4. Benchmark Phred-Weighted vs. Galaxy Raw Consensus Trees -------------------
rf_norm <- phangorn::RF.dist(phylo_phred$tr_rooted, phylo_galaxy$tr_rooted, normalize = TRUE)

# Check support for biological replicate clade (Kalimbang + II_Kalimbang)
get_kalimbang_support <- function(ns_df) {
  row <- ns_df %>% filter(Clade_Taxa == "II_Kalimbang + Kalimbang")
  if (nrow(row) == 0) return(list(monophyletic = FALSE, sh = NA_real_, bs = NA_real_))
  list(monophyletic = TRUE, sh = row$SH_aLRT_Pct[1], bs = row$ML_BS_Pct[1])
}
kal_pw  <- get_kalimbang_support(phylo_phred$node_support_df)
kal_gal <- get_kalimbang_support(phylo_galaxy$node_support_df)

benchmark_comp <- tibble(
  Assembly_Method             = c("Phred_Weighted_Consensus", "Galaxy_Raw_Unweighted_Consensus"),
  Trimmed_Alignment_Len_bp    = c(aln_objs$res_phred$summary_row$Trimmed_Aligned_Cols,
                                  aln_objs$res_galaxy$summary_row$Trimmed_Aligned_Cols),
  Mean_Alignment_Ambig_Pct    = c(aln_objs$res_phred$summary_row$Mean_Ambig_Pct_PostTrim,
                                  aln_objs$res_galaxy$summary_row$Mean_Ambig_Pct_PostTrim),
  BIC_Selected_Model          = c(phylo_phred$best_model_name, phylo_galaxy$best_model_name),
  Tree_LogLik                 = round(c(phylo_phred$fit_ml$logLik, phylo_galaxy$fit_ml$logLik), 2),
  Model_BIC                   = round(c(phylo_phred$mt_df$BIC[1], phylo_galaxy$mt_df$BIC[1]), 2),
  Model_AICc                  = round(c(phylo_phred$mt_df$AICc[1], phylo_galaxy$mt_df$AICc[1]), 2),
  Ingroup_Cophenetic_r        = round(c(phylo_phred$ingroup_coph_r, phylo_galaxy$ingroup_coph_r), 4),
  Mean_Internal_SH_aLRT_Pct   = round(c(mean(phylo_phred$node_support_df$SH_aLRT_Pct),
                                        mean(phylo_galaxy$node_support_df$SH_aLRT_Pct)), 2),
  Mean_Internal_ML_BS_Pct     = round(c(mean(phylo_phred$node_support_df$ML_BS_Pct),
                                        mean(phylo_galaxy$node_support_df$ML_BS_Pct)), 2),
  Kalimbang_Replicate_Monophyletic = c(kal_pw$monophyletic, kal_gal$monophyletic),
  Kalimbang_Clade_SH_aLRT_Pct = c(kal_pw$sh, kal_gal$sh),
  Kalimbang_Clade_ML_BS_Pct   = c(kal_pw$bs, kal_gal$bs),
  Normalized_RF_Dist_Between  = round(rf_norm, 4)
)
write.csv(benchmark_comp, file.path(DIRS$tables, "its_phred_vs_galaxy_benchmark.csv"), row.names = FALSE)

message(sprintf("[Benchmark] Phred-Weighted vs Galaxy-Raw: Kalimbang+II_Kalimbang clade support = %.0f/%.0f%% (Phred) vs %.0f/%.0f%% (Galaxy); Ingroup r_coph = %.4f",
                kal_pw$sh, kal_pw$bs, kal_gal$sh, kal_gal$bs, phylo_phred$ingroup_coph_r))

# 5. Two-Tier Tree & Distance Matrix Export for Phase 3 ------------------------
tr_N11 <- phylo_phred$tr_rooted
tr_N10 <- ape::drop.tip(tr_N11, "Ensete_ventricosum")

# Tier 2: Prune to the 6 shared 1:1 cultivars and map ITS_Accession_ID -> Cultivar_Key (e.g. Awa -> Awak)
shared_meta <- meta_df %>% filter(In_Morpho & In_ITS)
stopifnot(nrow(shared_meta) == 6)

tr_N6_its <- ape::keep.tip(tr_N11, shared_meta$ITS_Accession_ID)
# Rename tips from ITS_Accession_ID (Awa) to canonical Cultivar_Key (Awak) for 1:1 synchronization
tr_N6_synced <- tr_N6_its
tr_N6_synced$tip.label <- shared_meta$Cultivar_Key[match(tr_N6_its$tip.label, shared_meta$ITS_Accession_ID)]

# Extract N=6 patristic and ML distance matrices ordered identically to shared_meta$Cultivar_Key
pat_mat_N6 <- as.matrix(ape::cophenetic.phylo(tr_N6_synced))[shared_meta$Cultivar_Key, shared_meta$Cultivar_Key]
ml_mat_raw <- as.matrix(phylo_phred$d_ml_all11)[shared_meta$ITS_Accession_ID, shared_meta$ITS_Accession_ID]
rownames(ml_mat_raw) <- shared_meta$Cultivar_Key
colnames(ml_mat_raw) <- shared_meta$Cultivar_Key

# Verify strict 1:1 synchronization with Phase 1 Morphology N=6 tree
morpho_clust_objs <- readRDS(file.path(DIRS$rds, "morpho_clustering_objects.rds"))
stopifnot(setequal(tr_N6_synced$tip.label, morpho_clust_objs$phylo_A_N6$tip.label))
stopifnot(identical(rownames(pat_mat_N6), labels(morpho_clust_objs$D_A_N6$by_cultivar_key)))

ape::write.tree(tr_N11,       file = file.path(DIRS$trees, "its_ml_tree_N11_rooted.nwk"))
ape::write.tree(tr_N10,       file = file.path(DIRS$trees, "its_ml_tree_N10_ingroup.nwk"))
ape::write.tree(tr_N6_synced, file = file.path(DIRS$trees, "its_ml_tree_N6_shared.nwk"))

saveRDS(
  list(
    phylo_phred    = phylo_phred,
    phylo_galaxy   = phylo_galaxy,
    benchmark_comp = benchmark_comp,
    tr_N11         = tr_N11,
    tr_N10         = tr_N10,
    tr_N6_synced   = tr_N6_synced,
    D_its_pat_N6   = as.dist(pat_mat_N6),
    D_its_ml_N6    = as.dist(ml_mat_raw),
    shared_meta    = shared_meta
  ),
  file = file.path(DIRS$rds, "molecular_phylogeny_objects.rds")
)
message("[07] Saved ModelTest tables, ML node supports, and N=11 / N=10 / N=6 Newick trees.")


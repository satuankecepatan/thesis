# ==============================================================================
# R/run_concordance_pipeline.R
# Master Runner for Phase 3: Comparative Phenetic-Phylogenetic Concordance
# ==============================================================================

t0 <- Sys.time()
message("======================================================================")
message("  STARTING PHASE 3 PIPELINE: PHENETIC-PHYLOGENETIC CONCORDANCE")
message("======================================================================")

source("R/09_concordance_topological_metrics.R")
source("R/10_concordance_mantel_procrustes.R")
source("R/11_concordance_figures.R")

t1 <- Sys.time()
elapsed <- round(as.numeric(difftime(t1, t0, units = "secs")), 1)
message("======================================================================")
message(sprintf("  PHASE 3 PIPELINE COMPLETED SUCCESSFULLY IN %.1f SECONDS", elapsed))
message("======================================================================")


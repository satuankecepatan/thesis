# ==============================================================================
# R/run_morphology_pipeline.R
# Master Runner for Phase 1: Fruit Morphology Tree & Ordination Pipeline
# ==============================================================================

t_start <- Sys.time()
message("======================================================================")
message("  STARTING PHASE 1 PIPELINE: BANANA FRUIT MORPHOLOGY TREE ANALYSIS")
message("======================================================================")

source("R/00_config_aesthetics.R")
source("R/01_morpho_preprocess_qc.R")
source("R/02_morpho_pca_collinearity.R")
source("R/03_morpho_clustering_validation.R")
source("R/04_morpho_figures.R")

# Write deterministic session info log
writeLines(capture.output(sessionInfo()), file.path(DIRS$tables, "R_sessionInfo_phase1.txt"))

t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
message("======================================================================")
message(sprintf("  PHASE 1 PIPELINE COMPLETED SUCCESSFULLY IN %.1f SECONDS", t_elapsed))
message("======================================================================")


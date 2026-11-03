# ==============================================================================
# R/run_molecular_pipeline.R
# Master Runner for Phase 2: ITS Molecular Phylogenetics Pipeline
# ==============================================================================

t_start <- Sys.time()
message("======================================================================")
message("  STARTING PHASE 2 PIPELINE: BANANA ITS MOLECULAR PHYLOGENETICS")
message("======================================================================")

source("R/00_config_aesthetics.R")
source("R/05_molecular_consensus_qc.R")
source("R/06_molecular_alignment_trimming.R")
source("R/07_molecular_phylogeny_modeltest.R")
source("R/08_molecular_figures.R")

writeLines(capture.output(sessionInfo()), file.path(DIRS$tables, "R_sessionInfo_phase2.txt"))

t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
message("======================================================================")
message(sprintf("  PHASE 2 PIPELINE COMPLETED SUCCESSFULLY IN %.1f SECONDS", t_elapsed))
message("======================================================================")


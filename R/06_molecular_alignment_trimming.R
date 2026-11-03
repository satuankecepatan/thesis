# ==============================================================================
# R/06_molecular_alignment_trimming.R
# Phase 2 - Step 2: Orientation Verification, Iterative Local-Pair/Structure
# Alignment (DECIPHER::AlignSeqs), rRNA Domain Mapping, and Automated Trimming
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(ape)
  library(phangorn)
})
reset_project_seed()

message("=== [06] Starting Iterative Multiple Sequence Alignment & Automated Trimming ===")

dss_phred  <- Biostrings::readDNAStringSet(file.path(DIRS$processed, "its_consensus_phred_weighted.fasta"))
dss_galaxy <- Biostrings::readDNAStringSet(file.path(DIRS$processed, "its_consensus_galaxy_raw.fasta"))

# 1. Verify Orientation & Run Iterative Refinement Alignment -------------------
align_and_trim_pipeline <- function(dss_in, dataset_label) {
  reset_project_seed()
  # Verify 5' -> 3' orientation (Bioinformatics.md Rule #3)
  dss_oriented <- DECIPHER::OrientNucleotides(dss_in, verbose = FALSE)

  # Iterative refinement alignment (100 iterations, 100 refinements; Bioinformatics.md Rule #4)
  aln_full <- DECIPHER::AlignSeqs(
    dss_oriented,
    iterations  = 100,
    refinements = 100,
    processors  = 1,
    verbose     = FALSE
  )

  aln_mat <- as.matrix(aln_full)
  n_seq   <- nrow(aln_mat)
  n_col   <- ncol(aln_mat)

  # Locate conserved 18S, 5.8S, and 26S boundaries using consensus motifs
  cons_str <- paste(apply(aln_mat, 2, function(col) {
    v <- col[col %in% c("A", "C", "G", "T")]
    if (length(v) == 0) return("N")
    names(sort(table(v), decreasing = TRUE))[1]
  }), collapse = "")

  # 18S 3' end motif: GATCATTG (end of 18S rRNA)
  m_18s <- regexpr("GATCATTG", cons_str)
  end_18s <- if (m_18s[1] > 0) as.integer(m_18s[1] + attr(m_18s, "match.length") - 1L) else 45L

  # 5.8S rRNA start motif: CAACGGATATCT or CGATGAAGAACG
  m_58s_start <- regexpr("CAACGGATAT|ACTCTCGGCAAC", cons_str)
  start_58s <- if (m_58s_start[1] > 0) as.integer(m_58s_start[1]) else 255L

  # 5.8S rRNA 3' end motif: GAGGCCATCCGG or GCACGCCTGCCT
  m_58s_end <- regexpr("GCCATCCGG[A-Z]{0,15}AAGGGCACGCCTGC", cons_str)
  end_58s <- if (m_58s_end[1] > 0) as.integer(m_58s_end[1] + 12L) else (start_58s + 160L)

  # 26S rRNA 5' start motif: GACCCCAGGTCAGG
  m_26s <- regexpr("GACCCC[A-Z]{1,3}GTCAGG", cons_str)
  start_26s <- if (m_26s[1] > 0) as.integer(m_26s[1]) else (n_col - 55L)

  # Compute per-column gap fraction, ambiguity fraction, and Shannon information (bits)
  gap_frac   <- numeric(n_col)
  ambig_frac <- numeric(n_col)
  entropy_h  <- numeric(n_col)
  info_bits  <- numeric(n_col)
  is_var     <- logical(n_col)
  is_pars_inf<- logical(n_col)

  for (j in seq_len(n_col)) {
    col <- aln_mat[, j]
    gap_frac[j]   <- mean(col == "-")
    non_gap       <- col[col != "-"]
    ambig_frac[j] <- if (length(non_gap) > 0) mean(!non_gap %in% c("A", "C", "G", "T")) else 1.0
    acgt          <- non_gap[non_gap %in% c("A", "C", "G", "T")]

    if (length(acgt) >= 2) {
      freqs <- table(acgt) / length(acgt)
      h_val <- -sum(freqs * log2(freqs))
      entropy_h[j] <- h_val
      info_bits[j] <- max(0, 2.0 - h_val) * (1 - gap_frac[j])
      is_var[j]    <- length(freqs) >= 2
      is_pars_inf[j] <- sum(table(acgt) >= 2) >= 2
    } else {
      entropy_h[j] <- 0
      info_bits[j] <- 0
    }
  }

  # Programmatic Trimming (trimal -automated1 / ClipKIT logic; Bioinformatics.md Rule #5)
  # 1. Trim ragged 5' and 3' terminal overhangs until coverage >= 70% (gap_frac <= 0.30)
  covered_cols <- which(gap_frac <= 0.30)
  left_bound   <- min(covered_cols)
  right_bound  <- max(covered_cols)

  # 2. Within internal core [left_bound..right_bound], retain columns with gap_frac <= 0.35
  #    (allows at most 3/11 gaps so variable indels shared by clades are preserved while gappy artifacts are dropped)
  keep_col <- rep(FALSE, n_col)
  for (j in seq_len(n_col)) {
    if (j >= left_bound && j <= right_bound && gap_frac[j] <= 0.35) {
      keep_col[j] <- TRUE
    }
  }

  domain_vec <- case_when(
    seq_len(n_col) <= end_18s   ~ "18S rRNA (3')",
    seq_len(n_col) <  start_58s ~ "ITS1",
    seq_len(n_col) <= end_58s   ~ "5.8S rRNA",
    seq_len(n_col) <  start_26s ~ "ITS2",
    TRUE                        ~ "26S rRNA (5')"
  )

  col_audit <- tibble(
    Dataset           = dataset_label,
    Column_Index      = seq_len(n_col),
    rRNA_Domain       = domain_vec,
    Gap_Fraction      = round(gap_frac, 4),
    Ambig_Fraction    = round(ambig_frac, 4),
    Shannon_Entropy   = round(entropy_h, 4),
    Info_Content_Bits = round(info_bits, 4),
    Variable_Site     = is_var,
    Parsimony_Inform  = is_pars_inf,
    Retained          = keep_col,
    Trim_Reason       = case_when(
      keep_col ~ "Retained",
      seq_len(n_col) < left_bound  ~ "Trimmed (5' Terminal Overhang)",
      seq_len(n_col) > right_bound ~ "Trimmed (3' Terminal Overhang)",
      TRUE                         ~ "Trimmed (Internal Gappy Column > 35%)"
    )
  )

  # Construct trimmed DNAStringSet
  trimmed_mat <- aln_mat[, keep_col, drop = FALSE]
  trimmed_seqs <- apply(trimmed_mat, 1, paste, collapse = "")
  aln_trimmed  <- DNAStringSet(trimmed_seqs)

  summary_row <- tibble(
    Dataset                   = dataset_label,
    N_Taxa                    = n_seq,
    Raw_Aligned_Cols          = n_col,
    Trimmed_Aligned_Cols      = sum(keep_col),
    Cols_Removed              = sum(!keep_col),
    Pct_Cols_Retained         = round(100 * mean(keep_col), 2),
    Variable_Sites_PostTrim   = sum(is_var[keep_col]),
    Parsimony_Inform_PostTrim = sum(is_pars_inf[keep_col]),
    Mean_Gap_Pct_PostTrim     = round(100 * mean(gap_frac[keep_col]), 2),
    Mean_Ambig_Pct_PostTrim   = round(100 * mean(ambig_frac[keep_col]), 2)
  )

  list(
    aln_full    = aln_full,
    aln_trimmed = aln_trimmed,
    col_audit   = col_audit,
    summary_row = summary_row
  )
}

res_phred  <- align_and_trim_pipeline(dss_phred,  "Phred_Weighted")
res_galaxy <- align_and_trim_pipeline(dss_galaxy, "Galaxy_Raw")

# Export FASTA alignments (full and trimmed)
Biostrings::writeXStringSet(res_phred$aln_full,     file.path(DIRS$processed, "its_aligned_full_phred.fasta"))
Biostrings::writeXStringSet(res_phred$aln_trimmed,  file.path(DIRS$processed, "its_aligned_trimmed_phred.fasta"))
Biostrings::writeXStringSet(res_galaxy$aln_full,    file.path(DIRS$processed, "its_aligned_full_galaxy.fasta"))
Biostrings::writeXStringSet(res_galaxy$aln_trimmed, file.path(DIRS$processed, "its_aligned_trimmed_galaxy.fasta"))

audit_all <- bind_rows(res_phred$col_audit, res_galaxy$col_audit)
summ_all  <- bind_rows(res_phred$summary_row, res_galaxy$summary_row)

write.csv(audit_all, file.path(DIRS$tables, "its_alignment_trimming_audit.csv"), row.names = FALSE)
write.csv(summ_all,  file.path(DIRS$tables, "its_alignment_summary_metrics.csv"), row.names = FALSE)

message(sprintf("[Alignment & Trimming] Phred-Weighted: %d raw cols -> %d trimmed cols (%d variable, %d parsimony-informative, mean ambig = %.2f%%)",
                res_phred$summary_row$Raw_Aligned_Cols,
                res_phred$summary_row$Trimmed_Aligned_Cols,
                res_phred$summary_row$Variable_Sites_PostTrim,
                res_phred$summary_row$Parsimony_Inform_PostTrim,
                res_phred$summary_row$Mean_Ambig_Pct_PostTrim))

message(sprintf("[Alignment & Trimming] Galaxy-Raw:     %d raw cols -> %d trimmed cols (%d variable, %d parsimony-informative, mean ambig = %.2f%%)",
                res_galaxy$summary_row$Raw_Aligned_Cols,
                res_galaxy$summary_row$Trimmed_Aligned_Cols,
                res_galaxy$summary_row$Variable_Sites_PostTrim,
                res_galaxy$summary_row$Parsimony_Inform_PostTrim,
                res_galaxy$summary_row$Mean_Ambig_Pct_PostTrim))

saveRDS(
  list(
    res_phred  = res_phred,
    res_galaxy = res_galaxy,
    summ_all   = summ_all
  ),
  file = file.path(DIRS$rds, "molecular_alignment_objects.rds")
)
message("[06] Saved full and trimmed ITS alignments and column audit tables.")


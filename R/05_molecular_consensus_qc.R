# ==============================================================================
# R/05_molecular_consensus_qc.R
# Phase 2 - Step 1: Phred-Quality-Weighted Sanger Consensus Assembly,
# Outgroup Integration (Ensete ventricosum FJ428078.1), and Galaxy Raw Benchmark
# ==============================================================================

source("R/00_config_aesthetics.R")
suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(ape)
})
reset_project_seed()

message("=== [05] Starting Phred-Quality-Weighted ITS Consensus Assembly & QC ===")

f_dir      <- "data/galaxy/539_ FASTQ Groomer on Banana ITS (forward)/FASTQ Groomer on Banana ITS (forward)"
r_dir      <- "data/galaxy/562_ Reverse-Complement on Banana ITS (reverse)/Reverse-Complement on Banana ITS (reverse)"
gal_file   <- "data/galaxy/743_ Consensus sequence of Aligned Banana ITS_ Aligned sequences/Consensus sequence of Aligned Banana ITS_ Aligned sequences/Concatenated.fasta"
outgrp_fa  <- "data/raw/outgroup_FJ428078_Ensete_ventricosum.fasta"

# 1. Verify / Fetch Cached Outgroup (Ensete ventricosum FJ428078.1) ------------
if (!file.exists(outgrp_fa)) {
  message("[Outgroup] Fetching FJ428078.1 (Ensete ventricosum) from NCBI...")
  gb_dna <- ape::read.GenBank("FJ428078.1", as.character = TRUE)
  seq_str <- toupper(paste(gb_dna[[1]], collapse = ""))
  writeLines(c(">Ensete_ventricosum", seq_str), outgrp_fa)
}
outgrp_lines <- readLines(outgrp_fa, warn = FALSE)
outgrp_seq   <- toupper(paste(outgrp_lines[-1], collapse = ""))
stopifnot(nchar(outgrp_seq) > 500)

# 2. Read Galaxy Raw Consensus Sequences ---------------------------------------
gal_dss <- Biostrings::readDNAStringSet(gal_file)
accessions <- names(gal_dss)

# IUPAC ambiguity lookup for two bases of comparable high quality
iupac_pair <- function(b1, b2) {
  pair <- paste(sort(unique(c(b1, b2))), collapse = "")
  switch(pair,
    "A"  = "A", "C" = "C", "G" = "G", "T" = "T",
    "AG" = "R", "CT" = "Y", "CG" = "S", "AT" = "W",
    "GT" = "K", "AC" = "M", "N"
  )
}

# 3. Scaffolded Two-Regime Phred-Quality-Weighted Consensus Builder ------------
hutan_ref <- toupper(as.character(gal_dss[["Hutan"]]))

assemble_phred_weighted_consensus <- function(seq_F, qual_F, seq_R, qual_R, ref_seq = hutan_ref) {
  # Anchor F and R on the high-quality Musa reference scaffold so short reads
  # (e.g., Awang_Bangkal 76-bp F) align at their exact homologous coordinates.
  # The reference sequence is used strictly for coordinate scaffolding; all
  # consensus bases are called exclusively from F and R.
  dss <- DNAStringSet(c(Ref = toupper(ref_seq), F = toupper(seq_F), R = toupper(seq_R)))
  aln <- DECIPHER::AlignSeqs(dss, iterations = 50, refinements = 50, verbose = FALSE)

  chars_F <- strsplit(as.character(aln["F"]), "")[[1]]
  chars_R <- strsplit(as.character(aln["R"]), "")[[1]]
  n_col   <- length(chars_F)

  # Map Phred quality vectors onto alignment columns (NA at gaps '-')
  q_F_col <- rep(NA_real_, n_col)
  q_R_col <- rep(NA_real_, n_col)
  non_gap_F <- which(chars_F != "-")
  non_gap_R <- which(chars_R != "-")
  q_F_col[non_gap_F] <- qual_F
  q_R_col[non_gap_R] <- qual_R

  # Compute smoothed local window quality (window = +/- 7 columns)
  local_q <- function(q_vec) {
    out <- rep(NA_real_, length(q_vec))
    valid_idx <- which(!is.na(q_vec))
    if (length(valid_idx) == 0) return(out)
    for (j in seq_along(q_vec)) {
      lo <- max(1L, j - 7L)
      hi <- min(length(q_vec), j + 7L)
      w_vals <- q_vec[lo:hi]
      if (any(!is.na(w_vals))) out[j] <- mean(w_vals, na.rm = TRUE)
    }
    out
  }

  win_F <- local_q(q_F_col)
  win_R <- local_q(q_R_col)

  ov_start <- max(min(non_gap_F), min(non_gap_R))
  ov_end   <- min(max(non_gap_F), max(non_gap_R))

  cons_chars <- character(0)
  cons_quals <- numeric(0)

  for (j in seq_len(n_col)) {
    bF <- chars_F[j]
    bR <- chars_R[j]
    qF <- q_F_col[j]
    qR <- q_R_col[j]
    wF <- win_F[j]
    wR <- win_R[j]

    in_overlap <- (j >= ov_start && j <= ov_end)

    if (!in_overlap) {
      # Single-strand overhang (5' or 3' extension)
      if (bF != "-") {
        cons_chars <- c(cons_chars, bF)
        cons_quals <- c(cons_quals, qF)
      } else if (bR != "-") {
        cons_chars <- c(cons_chars, bR)
        cons_quals <- c(cons_quals, qR)
      }
    } else {
      # Inside overlap region: Two-Regime Phred Quality Resolution
      if (bF != "-" && bR != "-") {
        if (bF == bR) {
          cons_chars <- c(cons_chars, bF)
          cons_quals <- c(cons_quals, max(qF, qR, na.rm = TRUE))
        } else if (max(qF, qR) >= 28 && abs(qF - qR) >= 8) {
          # Regime 1: Asymmetric post-homopolymer stutter rescue
          # One strand is high-confidence (Q >= 28) and exceeds the other by >= 8 Phred units
          if (qF >= qR) {
            cons_chars <- c(cons_chars, bF)
            cons_quals <- c(cons_quals, qF)
          } else {
            cons_chars <- c(cons_chars, bR)
            cons_quals <- c(cons_quals, qR)
          }
        } else {
          # Regime 2: Dual low-quality hairpin zone or high-quality co-amplified rDNA paralogs
          # Preserve 2-state IUPAC ambiguity for Felsenstein partial-likelihood marginalization
          cons_chars <- c(cons_chars, iupac_pair(bF, bR))
          cons_quals <- c(cons_quals, max(qF, qR, na.rm = TRUE))
        }
      } else if (bF != "-" && bR == "-") {
        # Indel inside overlap: follow the strand with higher local window quality
        if (!is.na(wF) && !is.na(wR) && wF >= wR) {
          cons_chars <- c(cons_chars, bF)
          cons_quals <- c(cons_quals, qF)
        }
      } else if (bR != "-" && bF == "-") {
        if (!is.na(wR) && !is.na(wF) && wR >= wF) {
          cons_chars <- c(cons_chars, bR)
          cons_quals <- c(cons_quals, qR)
        }
      }
    }
  }

  cons_str <- paste(cons_chars, collapse = "")
  cons_str <- sub("^N+", "", sub("N+$", "", cons_str))
  list(
    seq       = cons_str,
    mean_qual = mean(cons_quals, na.rm = TRUE),
    q30_pct   = 100 * mean(cons_quals >= 30, na.rm = TRUE)
  )
}

# 4. Process All 10 Accessions & Compare Against Galaxy Raw --------------------
qc_rows   <- list()
phred_seqs <- character(length(accessions))
names(phred_seqs) <- accessions

for (i in seq_along(accessions)) {
  acc <- accessions[i]
  fl  <- readLines(file.path(f_dir, paste0(acc, ".fastqsanger")), warn = FALSE)
  rl  <- readLines(file.path(r_dir, paste0(acc, ".fastqsanger")), warn = FALSE)

  seq_F  <- fl[2]
  qual_F <- utf8ToInt(fl[4]) - 33
  seq_R  <- rl[2]
  qual_R <- utf8ToInt(rl[4]) - 33
  seq_G  <- toupper(as.character(gal_dss[[acc]]))

  res_pw <- assemble_phred_weighted_consensus(seq_F, qual_F, seq_R, qual_R)
  seq_PW <- res_pw$seq
  phred_seqs[acc] <- seq_PW

  ambig_G  <- nchar(gsub("[ACGT]", "", seq_G))
  ambig_PW <- nchar(gsub("[ACGT]", "", seq_PW))

  gc_cnt_PW <- nchar(gsub("[^GC]", "", seq_PW))
  at_cnt_PW <- nchar(gsub("[^AT]", "", seq_PW))
  gc_pct_PW <- 100 * gc_cnt_PW / max(1, gc_cnt_PW + at_cnt_PW)

  qc_rows[[i]] <- tibble(
    ITS_Accession_ID        = acc,
    Forward_Len_bp          = nchar(seq_F),
    Forward_Mean_Phred      = round(mean(qual_F), 2),
    Forward_Bases_Q30       = sum(qual_F >= 30),
    Reverse_Len_bp          = nchar(seq_R),
    Reverse_Mean_Phred      = round(mean(qual_R), 2),
    Reverse_Bases_Q30       = sum(qual_R >= 30),
    Galaxy_Raw_Len_bp       = nchar(seq_G),
    Galaxy_Ambig_Count      = ambig_G,
    Galaxy_Ambig_Pct        = round(100 * ambig_G / nchar(seq_G), 2),
    Phred_Weighted_Len_bp   = nchar(seq_PW),
    Phred_Weighted_Mean_Q   = round(res_pw$mean_qual, 2),
    Phred_Weighted_Q30_Pct  = round(res_pw$q30_pct, 2),
    Phred_Ambig_Count       = ambig_PW,
    Phred_Ambig_Pct         = round(100 * ambig_PW / nchar(seq_PW), 2),
    Ambig_Reduction_Pct     = ifelse(ambig_G > 0, round(100 * (ambig_G - ambig_PW) / ambig_G, 1), 0.0),
    GC_Content_Pct          = round(gc_pct_PW, 2)
  )
}

qc_df <- bind_rows(qc_rows)
write.csv(qc_df, file.path(DIRS$tables, "its_consensus_qc_comparison.csv"), row.names = FALSE)

message(sprintf("[Consensus QC] Total IUPAC ambiguities across 10 Musa accessions reduced from %d (Galaxy Raw) -> %d (Phred-Weighted) (%.1f%% reduction)!",
                sum(qc_df$Galaxy_Ambig_Count), sum(qc_df$Phred_Ambig_Count),
                100 * (sum(qc_df$Galaxy_Ambig_Count) - sum(qc_df$Phred_Ambig_Count)) / sum(qc_df$Galaxy_Ambig_Count)))

# 5. Export N=11 Multi-FASTA Files (10 Musa + Ensete_ventricosum Outgroup) -----
phred_all11 <- c(phred_seqs, Ensete_ventricosum = outgrp_seq)
gal_all11   <- c(sapply(as.character(gal_dss), toupper), Ensete_ventricosum = outgrp_seq)

Biostrings::writeXStringSet(
  DNAStringSet(phred_all11),
  filepath = file.path(DIRS$processed, "its_consensus_phred_weighted.fasta")
)

Biostrings::writeXStringSet(
  DNAStringSet(gal_all11),
  filepath = file.path(DIRS$processed, "its_consensus_galaxy_raw.fasta")
)

saveRDS(
  list(
    qc_df       = qc_df,
    phred_all11 = phred_all11,
    gal_all11   = gal_all11
  ),
  file = file.path(DIRS$rds, "molecular_consensus_objects.rds")
)
message("[05] Saved its_consensus_phred_weighted.fasta, its_consensus_galaxy_raw.fasta, and QC comparison table.")


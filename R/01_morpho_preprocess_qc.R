# ==============================================================================
# R/01_morpho_preprocess_qc.R
# Step 1: Programmatic Anomaly Detection, T1/T2 Descriptor Deduplication,
# and Consensus Replicate Merging (16 rows -> 10 unique cultivars)
# ==============================================================================

source("R/00_config_aesthetics.R")
reset_project_seed()

message("=== [01] Starting Morphology Preprocessing & QC ===")

# 1. Load Raw Phenotypic CSV & Curated Metadata --------------------------------
raw_csv_path <- "data/raw/Synced Fruit Characteristics.csv"
meta_path    <- "data/processed/cultivar_metadata.csv"

raw_df <- read.csv(raw_csv_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
meta_df <- read.csv(meta_path, stringsAsFactors = FALSE)

# Map each of the 16 rows to its canonical Newick-safe Cultivar_Key
raw_df <- raw_df %>%
  mutate(
    Row_ID       = gsub(" ", "_", `Variety Key`),
    Base_Name    = sub("_2$", "", `Variety Key`),
    Cultivar_Key = gsub(" ", "_", Base_Name),
    Replicate    = ifelse(grepl("_2$", `Variety Key`), "Rep_2 (T2)", "Rep_1 (T1)")
  )

stopifnot(all(raw_df$Cultivar_Key %in% meta_df$Cultivar_Key[meta_df$In_Morpho]))

# 2. Programmatic Unit-Anomaly Detection ---------------------------------------
# Inspect the 6 core continuous morphometric columns across all 16 rows
core_cont_cols <- c(
  "Fruit Length (cm)",
  "Fruit Width (cm)",
  "Fruit Stalk Length (cm)",
  "Fruit Stalk Width (cm)",
  "Diameter of Flesh With Skin (cm)",
  "Diameter of Flesh Without Skin (cm)"
)

anomalies_list <- list()

for (col_name in core_cont_cols) {
  vals <- raw_df[[col_name]]
  med_val <- median(vals, na.rm = TRUE)
  iqr_val <- IQR(vals, na.rm = TRUE)
  upper_fence <- med_val + 3 * iqr_val

  for (i in seq_len(nrow(raw_df))) {
    v <- vals[i]
    if (!is.na(v) && v > upper_fence) {
      ckey <- raw_df$Cultivar_Key[i]
      rkey <- raw_df$Row_ID[i]
      # Find paired replicate value if available
      pair_vals <- vals[raw_df$Cultivar_Key == ckey & raw_df$Row_ID != rkey]
      pair_val  <- if (length(pair_vals) == 1) pair_vals[1] else NA_real_
      ratio_to_pair <- if (!is.na(pair_val) && pair_val > 0) v / pair_val else NA_real_

      # Identify physical unit discrepancy mechanism
      if (col_name %in% c("Fruit Width (cm)", "Diameter of Flesh With Skin (cm)", "Diameter of Flesh Without Skin (cm)") &&
          !is.na(ratio_to_pair) && ratio_to_pair > 2.5) {
        mechanism <- "Circumference recorded instead of diameter (value approx pi * diameter)"
        factor_div <- pi
        corrected_v <- round(v / pi, 3)
      } else if (col_name == "Fruit Stalk Width (cm)" && v > 8.0) {
        mechanism <- "Millimeters (mm) recorded instead of centimeters (cm) (10x scale factor)"
        factor_div <- 10.0
        corrected_v <- round(v / 10.0, 3)
      } else {
        mechanism <- "Extreme morphometric outlier (> median + 3*IQR)"
        factor_div <- 1.0
        corrected_v <- v
      }

      anomalies_list[[length(anomalies_list) + 1]] <- tibble(
        Row_ID            = rkey,
        Cultivar_Key      = ckey,
        Trait_Column      = col_name,
        Raw_Value_cm      = v,
        Paired_Rep_cm     = pair_val,
        Ratio_To_Pair     = round(ratio_to_pair, 3),
        Sample_Median_cm  = round(med_val, 3),
        Sample_IQR_cm     = round(iqr_val, 3),
        Upper_3IQR_Fence  = round(upper_fence, 3),
        Suspected_Cause   = mechanism,
        Correction_Factor = round(factor_div, 4),
        Corrected_Val_cm  = corrected_v
      )
    }
  }
}

anomalies_df <- bind_rows(anomalies_list)
write.csv(anomalies_df, file.path(DIRS$tables, "morpho_qc_anomalies.csv"), row.names = FALSE)
message(sprintf("[QC] Flagged %d unit/scale anomalies across %d rows -> results/tables/morpho_qc_anomalies.csv",
                nrow(anomalies_df), length(unique(anomalies_df$Row_ID))))

# 3. Create Raw vs. QC-Corrected 16-Row Tables ---------------------------------
apply_clean_columns <- function(df, apply_unit_qc = TRUE) {
  out <- df
  if (apply_unit_qc) {
    for (k in seq_len(nrow(anomalies_df))) {
      r_id <- anomalies_df$Row_ID[k]
      c_nm <- anomalies_df$Trait_Column[k]
      out[out$Row_ID == r_id, c_nm] <- anomalies_df$Corrected_Val_cm[k]
    }
  }

  # Derive skin thickness for Kalimbang where [T1: K7] is missing using 0.5*(Diam_With_Skin - Diam_No_Skin)
  skin_t1 <- out[["Skin Thickness (cm) [T1: K7]"]]
  skin_calc <- pmax(0.05, round(0.5 * (out[["Diameter of Flesh With Skin (cm)"]] - out[["Diameter of Flesh Without Skin (cm)"]]), 3))
  skin_unified <- ifelse(is.na(skin_t1), skin_calc, skin_t1)

  # Reconcile overlapping T1/T2 categorical descriptors where one block is missing
  # Note: Stalk Surface [T1: V20] & [T2: KT5] have 100% agreement on shared cultivars (except Kalimbang only in T2)
  stalk_surf_unified <- ifelse(!is.na(out[["Fruit Stalk Surface [T1: V20]"]]),
                               out[["Fruit Stalk Surface [T1: V20]"]],
                               out[["Fruit Stalk Surface [T2: KT5]"]])

  skin_young_unified <- ifelse(!is.na(out[["Skin Color When Young [T1: V21]"]]),
                               out[["Skin Color When Young [T1: V21]"]],
                               out[["Skin Color When Young [T2: KT6]"]])

  fruit_tip_unified  <- ifelse(!is.na(out[["Fruit Tip [T1: V19]"]]),
                               out[["Fruit Tip [T1: V19]"]],
                               out[["Fruit Tip [T2: KT1]"]])

  fruit_pos_unified  <- ifelse(!is.na(out[["Fruit Position [T1: V17]"]]),
                               out[["Fruit Position [T1: V17]"]],
                               out[["Fruit Position [T2: KT17]"]])

  fruit_shape_unified <- ifelse(!is.na(out[["Fruit Shape [T1: V18]"]]),
                                out[["Fruit Shape [T1: V18]"]],
                                out[["Fruit Shape [T2: KT20]"]])

  tibble(
    Row_ID              = out$Row_ID,
    Cultivar_Key        = out$Cultivar_Key,
    Replicate           = out$Replicate,
    # --- Core Complete Continuous Morphometrics (7 traits across all 10 cultivars) ---
    Fruit_Length_cm     = as.numeric(out[["Fruit Length (cm)"]]),
    Fruit_Width_cm      = as.numeric(out[["Fruit Width (cm)"]]),
    Stalk_Length_cm     = as.numeric(out[["Fruit Stalk Length (cm)"]]),
    Stalk_Width_cm      = as.numeric(out[["Fruit Stalk Width (cm)"]]),
    Diam_With_Skin_cm   = as.numeric(out[["Diameter of Flesh With Skin (cm)"]]),
    Diam_No_Skin_cm     = as.numeric(out[["Diameter of Flesh Without Skin (cm)"]]),
    Skin_Thickness_cm   = as.numeric(skin_unified),
    # --- Partial Continuous / Count Traits ---
    Fruit_Weight_g      = as.numeric(out[["Fruit Weight (g) [T2: KN1]"]]),
    Hands_Per_Bunch     = as.numeric(out[["Number of Hands per Bunch [T1: V27]"]]),
    Fruits_Per_Hand_Cat = as.numeric(out[["Number of Fruits per Hand [T1: V28]"]]),
    # --- Harmonized Categorical / Ordinal IPGRI Descriptors (0% missing across 10 cultivars) ---
    Fruit_Position_Cat  = as.integer(fruit_pos_unified),
    Fruit_Shape_Cat     = as.integer(fruit_shape_unified),
    Fruit_Tip_Cat       = as.integer(fruit_tip_unified),
    Stalk_Surface_Cat   = as.integer(stalk_surf_unified),
    Skin_Young_Cat      = as.integer(skin_young_unified),
    # --- Source-Specific Categorical Descriptors (handled by Gower's daisy with NA support) ---
    Bunch_Position_T1   = as.integer(out[["Bunch Position [T1: V15]"]]),
    Bunch_Shape_T1      = as.integer(out[["Bunch Shape [T1: V16]"]]),
    Relic_Tip_T2        = as.integer(out[["Relic Waste at Fruit Tip [T2: KT2]"]]),
    Skin_Ripe_T2        = as.integer(out[["Skin Color When Ripe [T2: KT7]"]]),
    Skin_Adhesion_T2    = as.integer(out[["Fruit Skin Adhesion [T2: KT9]"]]),
    Skin_Cracks_T2      = as.integer(out[["Cracks on Fruit Skin [T2: KT10]"]]),
    Pulp_Presence_T2    = as.integer(out[["Pulp in Fruit [T2: KT11]"]]),
    Pulp_Unripe_T2      = as.integer(out[["Pulp Color Before Ripe [T2: KT12]"]]),
    Flesh_Color_T2      = as.integer(out[["Fruit Flesh Color [T2: KT13]"]]),
    Flesh_Texture_T2    = as.integer(out[["Flesh Texture [T2: KT14]"]]),
    Dominant_Taste_T2   = as.integer(out[["Dominant Taste [T2: KT14]"]]),
    Seed_Surface_T2     = as.integer(out[["Seed Surface [T2: KT15]"]]),
    Seed_Shape_T2       = as.integer(out[["Seed Shape [T2: KT16]"]]),
    Fruits_Bunch_Cat_T2 = as.integer(out[["Number of Fruits in 1 Bunch [T2: KT18]"]]),
    Cross_Section_T2    = as.integer(out[["Cross-Sectional Shape [T2: KT21]"]])
  )
}

all16_qc  <- apply_clean_columns(raw_df, apply_unit_qc = TRUE)
all16_raw <- apply_clean_columns(raw_df, apply_unit_qc = FALSE)

write.csv(all16_qc,  file.path(DIRS$processed, "morpho_all16_qc.csv"),  row.names = FALSE)
write.csv(all16_raw, file.path(DIRS$processed, "morpho_all16_raw.csv"), row.names = FALSE)

# 4. Merge _2 Replicate Pairs into 10 Consensus Cultivars ----------------------
cont_vars <- c(
  "Fruit_Length_cm", "Fruit_Width_cm", "Stalk_Length_cm", "Stalk_Width_cm",
  "Diam_With_Skin_cm", "Diam_No_Skin_cm", "Skin_Thickness_cm",
  "Fruit_Weight_g", "Hands_Per_Bunch"
)

cat_vars <- setdiff(names(all16_qc), c("Row_ID", "Cultivar_Key", "Replicate", cont_vars))

first_non_na <- function(x) {
  vals <- x[!is.na(x)]
  if (length(vals) == 0) return(NA_integer_)
  vals[1]
}

build_consensus_10 <- function(df16) {
  # Preserve order in cultivar_metadata.csv
  morpho_meta <- meta_df %>% filter(In_Morpho)

  cons <- df16 %>%
    group_by(Cultivar_Key) %>%
    summarise(
      N_Replicates = n(),
      across(all_of(cont_vars), ~ ifelse(all(is.na(.x)), NA_real_, round(mean(.x, na.rm = TRUE), 4))),
      across(all_of(cat_vars),  ~ first_non_na(.x)),
      .groups = "drop"
    )

  morpho_meta %>%
    select(Cultivar_Key, Display_Name, Plot_Label_Markdown, Plot_Label_Plain,
           Genome_Group, In_ITS, ITS_Accession_ID) %>%
    inner_join(cons, by = "Cultivar_Key")
}

consensus_qc  <- build_consensus_10(all16_qc)
consensus_raw <- build_consensus_10(all16_raw)

stopifnot(nrow(consensus_qc) == 10, nrow(consensus_raw) == 10)
stopifnot(sum(is.na(consensus_qc[, c("Fruit_Length_cm", "Fruit_Width_cm", "Stalk_Length_cm",
                                     "Stalk_Width_cm", "Diam_With_Skin_cm", "Diam_No_Skin_cm",
                                     "Skin_Thickness_cm")])) == 0)

write.csv(consensus_qc,  file.path(DIRS$processed, "morpho_consensus_qc.csv"),  row.names = FALSE)
write.csv(consensus_raw, file.path(DIRS$processed, "morpho_consensus_raw.csv"), row.names = FALSE)

# 5. Export Structured Trait Dictionary ----------------------------------------
all_clean_vars <- c(cont_vars, cat_vars)
trait_dict <- tibble(
  Clean_Variable = all_clean_vars,
  Data_Type = case_when(
    Clean_Variable == "Hands_Per_Bunch" ~ "Discrete Count (Ratio)",
    Clean_Variable %in% cont_vars ~ "Continuous (Ratio)",
    grepl("_T1$", Clean_Variable) ~ "Categorical / Ordinal (T1 Block)",
    grepl("_T2$", Clean_Variable) ~ "Categorical / Ordinal (T2 Block)",
    TRUE ~ "Categorical / Ordinal (Harmonized T1+T2)"
  ),
  Missing_Count_In_10_Cultivars = as.integer(sapply(consensus_qc[, all_clean_vars], function(x) sum(is.na(x)))),
  Included_In_Track_A_Continuous = Clean_Variable %in% cont_vars[1:7],
  Included_In_Track_B_Gower      = TRUE
)

write.csv(trait_dict, file.path(DIRS$tables, "trait_dictionary.csv"), row.names = FALSE)
message("[01] Saved 10-cultivar consensus tables (morpho_consensus_qc.csv & morpho_consensus_raw.csv) and trait_dictionary.csv.")


# ==============================================================================
# R/00_config_aesthetics.R
# Deterministic Reproducibility, Directory Setup, and Locked Aesthetics Config
# Follows .agent/RULES/AGENT.md, .agent/SKILL/Aesthetics.md, Bioinformatics.md
# ==============================================================================

# Ensure user R library is available as fallback alongside renv for Bioconductor packages (Biostrings, DECIPHER)
user_r_lib <- "C:/Users/Asus/AppData/Local/R/win-library/4.5"
if (dir.exists(user_r_lib) && !(user_r_lib %in% .libPaths())) {
  .libPaths(c(.libPaths(), user_r_lib))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggtext)
  library(patchwork)
  library(cowplot)
  library(viridis)
})

# 1. Deterministic Reproducibility (AGENT.md Rule #4) --------------------------
SEED_FILE <- "seed.txt"
if (!file.exists(SEED_FILE)) {
  stop("seed.txt not found in project root. Deterministic seed is required.")
}
PROJECT_SEED <- as.integer(trimws(readLines(SEED_FILE, warn = FALSE)[1]))
stopifnot(!is.na(PROJECT_SEED))

suppressWarnings(RNGversion("3.6.0"))
set.seed(PROJECT_SEED)

reset_project_seed <- function() {
  suppressWarnings(RNGversion("3.6.0"))
  set.seed(PROJECT_SEED)
}

# 2. Output Directory Hierarchy ------------------------------------------------
DIRS <- list(
  processed = "data/processed",
  tables    = "results/tables",
  figures   = "results/figures",
  rds       = "results/rds",
  trees     = "results/trees"
)
for (d in DIRS) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# 3. Locked Project-Wide Color Mapping (Aesthetics.md Rules #1 & #2) -----------
# Okabe-Ito colorblind-safe palette locked across PCA, Morpho Dendrogram,
# ITS Phylogram, and Tanglegram.
GENOME_PALETTE <- c(
  "AA"              = "#E69F00", # Okabe-Ito Orange (Diploid acuminata: Mas)
  "AAA"             = "#0072B2", # Okabe-Ito Blue (Triploid acuminata: Ambon)
  "AAB"             = "#009E73", # Okabe-Ito Bluish Green (Hybrid triploid: Raja, Kapas)
  "ABB"             = "#CC79A7", # Okabe-Ito Reddish Purple (Balbisiana triploid: Kepok, Awak)
  "AA/AAB (Local)"  = "#56B4E9", # Okabe-Ito Sky Blue (South Kalimantan local cultivars)
  "Wild / Outgroup" = "#D55E00"  # Okabe-Ito Vermilion (Wild Musa outgroup: Hutan)
)

# Shape mapping to distinguish cultivars in both datasets vs morphology-only
DATASET_SHAPES <- c(
  "Matched (Morpho + ITS)" = 19, # Solid circle
  "Morphology Only"        = 17  # Solid triangle
)

# 4. Publication Typography & ggplot2 Theme (Aesthetics.md Rules #3, #4) -------
theme_thesis <- function(base_size = 10, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title          = ggtext::element_markdown(size = base_size + 1.5, face = "bold", hjust = 0, margin = margin(b = 4)),
      plot.subtitle       = ggtext::element_markdown(size = base_size - 0.5, color = "#333333", hjust = 0, margin = margin(b = 8)),
      plot.caption        = ggtext::element_markdown(size = base_size - 1.5, color = "#555555", hjust = 1, margin = margin(t = 6)),
      plot.tag            = element_text(size = 12, face = "bold", color = "black"),
      axis.title          = ggtext::element_markdown(size = base_size + 0.5, face = "bold", color = "black"),
      axis.text           = element_text(size = base_size - 1, color = "#111111"),
      legend.title        = ggtext::element_markdown(size = base_size - 0.5, face = "bold"),
      legend.text         = element_text(size = base_size - 1),
      legend.background   = element_rect(fill = "white", color = "#CCCCCC", linewidth = 0.3),
      panel.grid.minor    = element_blank(),
      panel.grid.major    = element_line(color = "#EFEFEF", linewidth = 0.3),
      panel.border        = element_rect(color = "#222222", linewidth = 0.6),
      strip.background    = element_rect(fill = "#F5F7FA", color = "#222222", linewidth = 0.5),
      strip.text          = element_text(size = base_size - 0.5, face = "bold")
    )
}

# 5. Dual-Format Publication Figure Exporter (Aesthetics.md Rule #7) -----------
save_thesis_plot <- function(plot_obj, basename_no_ext, width = 10, height = 7, dpi = 600) {
  pdf_path <- file.path(DIRS$figures, paste0(basename_no_ext, ".pdf"))
  png_path <- file.path(DIRS$figures, paste0(basename_no_ext, ".png"))
  svg_path <- file.path(DIRS$figures, paste0(basename_no_ext, ".svg"))

  tryCatch({
    ggsave(
      filename = pdf_path,
      plot     = plot_obj,
      width    = width,
      height   = height,
      units    = "in",
      device   = cairo_pdf
    )
  }, error = function(e) {
    ggsave(
      filename = pdf_path,
      plot     = plot_obj,
      width    = width,
      height   = height,
      units    = "in",
      device   = "pdf"
    )
  })
  ggsave(
    filename = png_path,
    plot     = plot_obj,
    width    = width,
    height   = height,
    units    = "in",
    dpi      = dpi,
    bg       = "white"
  )
  # Also attempt SVG if svglite is available, fallback silently
  tryCatch({
    ggsave(
      filename = svg_path,
      plot     = plot_obj,
      width    = width,
      height   = height,
      units    = "in"
    )
  }, error = function(e) NULL)

  message(sprintf("[Saved Figure] %s (.pdf & %d DPI .png)", basename_no_ext, dpi))
  invisible(list(pdf = pdf_path, png = png_path))
}


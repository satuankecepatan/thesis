# Project R Dependencies based on README.md and .agent rules
#
# 1. Morphology Tree Analysis (Fruit Traits)
library(cluster)       # Gower's distance (daisy), clustering algorithms
library(pvclust)       # Bootstrap cluster validation (AU / BP p-values)
library(FactoMineR)    # PCA and collinearity analysis
library(factoextra)    # PCA visualization and extraction

# 2. Molecular Phylogenetics (ITS Region)
library(ape)           # Phylogenetic analysis, tree I/O, cophenetic distances
library(phangorn)      # Robinson-Foulds distance (RF.dist), phylogenetic models
library(treeio)        # Tree data structures and metadata import
library(ggtree)        # Phylogenetic tree visualization

# 3. Tree vs. Tree Showdown & Statistics
library(dendextend)    # Tanglegram, entanglement, Baker's Gamma (cor_bakers_gamma)
library(vegan)         # Mantel test (mantel), Procrustes superimposition (protest)
library(Quartet)       # Quartet divergence and topological congruence

# 4. Data Science, Aesthetics, and Publication Graphics
library(tidyverse)     # Data manipulation (dplyr, tibble) and plotting (ggplot2)
library(ggrepel)       # Non-overlapping text labels
library(ggtext)        # Markdown typography in plot labels and titles
library(ggsci)         # Scientific journal color palettes
library(patchwork)     # Programmatic multi-panel plot assembly
library(cowplot)       # Plot publication styling
library(viridis)       # Colorblind-safe color palettes


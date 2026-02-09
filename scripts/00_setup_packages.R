# ======================================================
# Script: 00_setup_packages.R
# Project: Coral 16S Metabarcoding
# Author: Salome Chemutai
# Description: Install and load required R packages
# ======================================================

# List of required packages
packages <- c(
  "dada2",       # Denoising, chimera removal
  "phyloseq",    # Microbiome data handling
  "tidyverse",   # Data wrangling and plotting
  "ggplot2",     # Data visualization
  "vegan",       # Diversity metrics
  "Biostrings",  # DNA sequence manipulation
  "DECIPHER",    # Sequence alignment & taxonomic assignment
  "phangorn"     # Phylogenetic analysis
)

# Install packages that are not already installed
installed_packages <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load all packages
lapply(packages, library, character.only = TRUE)

# Optional: Check that all packages are loaded successfully
cat("All packages loaded successfully!\n")


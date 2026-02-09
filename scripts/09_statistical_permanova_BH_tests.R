### 09_statistical_permanova_BH_tests.R ###
# PERMANOVA tests for coral microbiome experiments/coral species

#### . Coral Bleaching: Treatment/Coral species effect ####

# Subset Samples
ps.bleaching <- subset_samples(
  ps.cleaned_N,
  Coral_species %in% c("Acropora_tenuis", "Acropora_verweyi") &
    Experiment == "Coral_bleaching"
)

# Remove taxa with zero abundance
ps.bleaching <- prune_taxa(taxa_sums(ps.bleaching) > 0, ps.bleaching)

# Agglomerate to Genus level
ps.genus <- tax_glom(ps.bleaching, taxrank = "Genus")

# Calculate Bray-Curtis distance
bray_dist <- phyloseq::distance(ps.genus, method = "bray")

# Extract metadata
metadata <- as(sample_data(ps.genus), "data.frame")

# Run PERMANOVA
set.seed(123)
permanova_result <- adonis2(
  bray_dist ~ Treatment, #(or Coral_species)
  data = metadata,
  permutations = 999
)

# Adjust p-values
raw_pvals <- permanova_result$`Pr(>F)`
adjusted_pvals <- p.adjust(raw_pvals, method = "BH")

# Summary table
permanova_summary <- data.frame(
  Factor = rownames(permanova_result),
  R2 = permanova_result$R2,
  Raw_P = raw_pvals,
  Adjusted_P = adjusted_pvals
)

print(permanova_summary)

#### Intertidal_2022: Treatment/Coral species effect ####

# Subset Samples
ps.intertidal <- subset_samples(
  ps.cleaned_N,
  Coral_species %in% c("Acropora_tenuis", "Stylophora_pistillata") &
    Experiment == "Intertidal_2022"
)

# Remove taxa with zero abundance
ps.intertidal <- prune_taxa(taxa_sums(ps.intertidal) > 0, ps.intertidal)

# Agglomerate to Genus level
ps.genus <- tax_glom(ps.intertidal, taxrank = "Genus")

# Calculate Bray-Curtis distance
bray_dist <- phyloseq::distance(ps.genus, method = "bray")

# Extract metadata
metadata <- as(sample_data(ps.genus), "data.frame")

# Run PERMANOVA
set.seed(123)
permanova_result <- adonis2(
  bray_dist ~ Treatment, #(or Coral_species)
  data = metadata,
  permutations = 999
)

# Adjust p-values
raw_pvals <- permanova_result$`Pr(>F)`
adjusted_pvals <- p.adjust(raw_pvals, method = "BH")

# Summary table
permanova_summary <- data.frame(
  Factor = rownames(permanova_result),
  R2 = permanova_result$R2,
  Raw_P = raw_pvals,
  Adjusted_P = adjusted_pvals
)

print(permanova_summary)

#### Nursery vs Natural Reef: Structure/Coral species effect ####

# Subset Samples
ps.nursery_natural <- subset_samples(
  ps.cleaned_N,
  Structure %in% c("Artificial_structure", "Natural_reef") &
    Experiment == "Intertidal_2022"
)

# Remove taxa with zero abundance
ps.nursery_natural <- prune_taxa(taxa_sums(ps.nursery_natural) > 0, ps.nursery_natural)

# Agglomerate to Genus level
ps.genus <- tax_glom(ps.nursery_natural, taxrank = "Genus")

# Calculate Bray-Curtis distance
bray_dist <- phyloseq::distance(ps.genus, method = "bray")

# Extract metadata
metadata <- as(sample_data(ps.genus), "data.frame")

# Run PERMANOVA
set.seed(123)
permanova_result <- adonis2(
  bray_dist ~ Structure, #(or Coral_species)
  data = metadata,
  permutations = 999
)

# Adjust p-values
raw_pvals <- permanova_result$`Pr(>F)`
adjusted_pvals <- p.adjust(raw_pvals, method = "BH")

# Summary table
permanova_summary <- data.frame(
  Factor = rownames(permanova_result),
  R2 = permanova_result$R2,
  Raw_P = raw_pvals,
  Adjusted_P = adjusted_pvals
)

print(permanova_summary)






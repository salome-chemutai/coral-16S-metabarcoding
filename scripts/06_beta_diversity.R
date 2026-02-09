# 06_beta_diversity.R
# Genus-level Bray PCoAs for Coral Bleaching, Intertidal 2022, and Nursery vs Natural Reef

# Load cleaned phyloseq object
ps <- readRDS("ps.cleaned_N.rds")

# ===============================================
### 1. Coral Bleaching ###
# ===============================================

# Subset samples
ps.bleaching <- subset_samples(
  ps.cleaned_N,
  Coral_species %in% c("Acropora_tenuis", "Acropora_verweyi") &
    Experiment == "Coral_bleaching"
)

# Remove zero-abundance taxa
ps.bleaching <- prune_taxa(taxa_sums(ps.bleaching) > 0, ps.bleaching)

# Agglomerate to Genus level
ps.genus.bleaching <- tax_glom(ps.bleaching, taxrank = "Genus")
ps.genus.bleaching <- prune_taxa(taxa_sums(ps.genus.bleaching) > 0, ps.genus.bleaching)

# Update metadata
meta <- as(sample_data(ps.genus.bleaching), "data.frame")
meta$Treatment <- ifelse(meta$Treatment == "Not_bleached", "Non_bleached", meta$Treatment)
sample_data(ps.genus.bleaching)$Treatment <- meta$Treatment

# Distance matrix & ordination
bray_dist_genus <- phyloseq::distance(ps.genus.bleaching, method = "bray")

# Perform PCoA
ord.pcoa.genus <- ordinate(ps.genus.bleaching, method = "PCoA", distance = bray_dist_genus)

# Plot Ordination
plot_ordination(ps.genus.bleaching, ord.pcoa.genus,
                color = "Treatment", shape = "Coral_species") +
  geom_point(size = 6, alpha = 0.9) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Genus-level PCoA - Coral Bleaching",
    color = "Treatment",
    shape = "Coral Species"
  ) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )


# ===============================================
### Intertidal 2022 ###
# ===============================================

# Subset samples
ps.intertidal <- subset_samples(
  ps.cleaned_N,
  Coral_species %in% c("Acropora_tenuis", "Stylophora_pistillata") &
    Experiment == "Intertidal_2022"
)

# Remove zero-abundance taxa
ps.intertidal <- prune_taxa(taxa_sums(ps.intertidal) > 0, ps.intertidal)

# Agglomerate to Genus level
ps.genus.intertidal <- tax_glom(ps.intertidal, taxrank = "Genus")
ps.genus.intertidal <- prune_taxa(taxa_sums(ps.genus.intertidal) > 0, ps.genus.intertidal)

# Update metadata
meta <- as(sample_data(ps.genus.intertidal), "data.frame")
sample_data(ps.genus.intertidal)$Treatment <- meta$Treatment

# Distance matrix & ordination
bray_dist_genus <- phyloseq::distance(ps.genus.intertidal, method = "bray")

# Perform PCoA
ord.pcoa.genus <- ordinate(ps.genus.intertidal, method = "PCoA", distance = bray_dist_genus)

# Plot Ordination
plot_ordination(ps.genus.intertidal, ord.pcoa.genus,
                color = "Treatment", shape = "Coral_species") +
  geom_point(size = 6, alpha = 0.9) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Genus-level PCoA - Intertidal 2022",
    color = "Treatment",
    shape = "Coral Species"
  ) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )


# ===============================================
### Nursery vs Natural Reef ###
# ===============================================

# Subset samples
ps.structure <- subset_samples(
  ps.cleaned_N,
  Structure %in% c("Artificial_structure", "Natural_reef") &
    Experiment == "Intertidal_2022"
)

# Remove zero-abundance taxa
ps.structure <- prune_taxa(taxa_sums(ps.structure) > 0, ps.structure)

# Agglomerate to Genus level
ps.genus.structure <- tax_glom(ps.structure, taxrank = "Genus")
ps.genus.structure <- prune_taxa(taxa_sums(ps.genus.structure) > 0, ps.genus.structure)

# Update metadata
meta <- as(sample_data(ps.genus.structure), "data.frame")
meta$Treatment <- ifelse(meta$Structure == "Artificial_structure", "Nursery", meta$Structure)
sample_data(ps.genus.structure)$Treatment <- meta$Treatment

# Distance matrix & ordination
dist.structure <- phyloseq::distance(ps.genus.structure, method = "bray")

# Perform PCoA 
ord.pcoa.genus <- ordinate(ps.genus.intertidal, method = "PCoA", distance = bray_dist_genus)

# Plot Ordination
plot_ordination(ps.genus.structure, ord.pcoa.genus,
                color = "Treatment", shape = "Coral_species") +
  geom_point(size = 6, alpha = 0.9) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Genus-level NMDS - Nursery vs Natural Reef",
    color = "Treatment",
    shape = "Coral Species"
  ) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )


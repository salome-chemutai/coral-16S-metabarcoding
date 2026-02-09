# 07_relative_abundance_phyla.R
# Phylum-level relative abundance plots for coral bleaching, intertidal experiment and nursery vs natural reef treatments/ coral species

# Load the cleaned phyloseq object
ps <- readRDS("ps.cleaned_N.rds")

### Coral Bleaching Experiment ###

# Subset to coral bleaching samples 
ps.bleaching <- subset_samples(
 ps.cleaned_N,
  Experiment == "Coral_bleaching" &
    Coral_species %in% c("Acropora_tenuis", "Acropora_verweyi") &
    Treatment %in% c("Not_bleached", "Bleached")
)

# Prune zero-abundance taxa
ps.bleaching <- prune_taxa(taxa_sums(ps.bleaching) > 0, ps.bleaching)

# Agglomerate to Phylum level
ps.phylum <- tax_glom(ps.bleaching, taxrank = "Phylum")

# Transform to relative abundance
ps.rel <- transform_sample_counts(ps.phylum, function(x) x / sum(x))

# Melt to long format
df <- psmelt(ps.rel)

# Top 10 phyla
top10 <- df %>%
  group_by(Phylum) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# Filter for top 10 phyla only
  df_top <- df %>% filter(Phylum %in% top10)

# Summarize and normalize within each treatment
  df_summed <- df_top %>%
  group_by(Treatment, Phylum) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Treatment) %>%
  mutate(Relative_Abundance = total_abundance / sum(total_abundance)) %>%
  ungroup()

# Rename "Not_bleached" to "Non_bleached"
  df_summed$Treatment <- recode(df_summed$Treatment, "Not_bleached" = "Non_bleached")

# Plot
 ggplot(df_summed, aes(x = Treatment, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = palette_colors) +
  labs(
    title ="Top 10 Phyla in  Coral Bleaching Experiment",
       x = "Structure", y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  )                          

### Intertidal 2022 Experiment ###
                                  
# Subset to intertidal 2022 samples 
  ps.intertidal <- subset_samples(
  ps.cleaned_N,
  Experiment == "Intertidal_2022" &
    Treatment %in% c("Natural_intertidal", "Intertidal_intertidal", "Intertidal_subtidal",
                     "Natural_subtidal", "Subtidal_subtidal", "Subtidal_intertidal")
)
# Prune zero-abundance taxa
ps.intertidal <- prune_taxa(taxa_sums(ps.intertidal) > 0, ps.intertidal)

# Agglomerate to Phylum level
ps.phylum <- tax_glom(ps.intertidal, taxrank = "Phylum")

# Transform to relative abundance
ps.rel <- transform_sample_counts(ps.phylum, function(x) x / sum(x))

# Melt to long format
df <- psmelt(ps.rel)

# Top 10 phyla                        
top10 <- df %>%
  group_by(Phylum) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)
                                  
# Filter for top 10 phyla only
  df_top <- df %>% filter(Phylum %in% top10)

# Summarize and normalize within each treatment
df_summed <- df_top %>%
  group_by(Treatment, Phylum) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Treatment) %>%
  mutate(Relative_Abundance = total_abundance / sum(total_abundance)) %>%
  ungroup()

# Plot
  ggplot(df_summed, aes(x = Treatment, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = palette_colors) +
  labs(
    title = "Top 10 Phyla in Intertidal 2022 Experiment",
    x = "Treatment",
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  )


### Nursery vs Natural Reef ###
                                  
# Subset samples
  ps.nursery_natural <- subset_samples(
  ps.cleaned_N,
  Structure %in% c("Artificial_structure", "Natural_reef") &
    Experiment == "Intertidal_2022"
)
# Prune zero-abundance taxa
ps.nursery_natural <- prune_taxa(taxa_sums(ps.nursery_natural) > 0, ps.nursery_natural)

# Agglomerate to Phylum level
ps.phylum <- tax_glom(ps.nursery_natural, taxrank = "Phylum")

# Transform to relative abundance
ps.rel <- transform_sample_counts(ps.phylum, function(x) x / sum(x))

# Melt to long format
df <- psmelt(ps.rel)

# Top ten    
top10 <- df %>%
  group_by(Phylum) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# Filter for top 10 phyla only
  df_top <- df %>% filter(Phylum %in% top10)

# Summarize and normalize within each treatment                                  
df_summed <- df_top %>%
  group_by(Structure, Phylum) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Structure) %>%
  mutate(Relative_Abundance = total_abundance / sum(total_abundance)) %>%
  ungroup()

# Rename "Artificial_structure" to "Nursery"
df_summed$Structure <- recode(df_summed$Structure, "Artificial_structure" = "Nursery")

# Plot
 ggplot(df_summed, aes(x = Treatment, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = palette_colors) +
  labs(
    title ="Top 10 Phyla in Nursery vs Natural Reef",
       x = "Structure", y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  )

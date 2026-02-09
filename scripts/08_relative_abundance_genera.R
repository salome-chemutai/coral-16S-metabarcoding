# 08_relative_abundance_genera.R
# Genus-level relative abundance plots for all experiments/coral species

# Load cleaned phyloseq object
ps <- readRDS("ps.cleaned_N.rds")

### Coral Bleaching Experiment ###

# Subset samples
ps.target <- subset_samples(
  ps.cleaned_N,
  Experiment == "Coral_bleaching" &
    Coral_species %in% c("Acropora_tenuis", "Acropora_verweyi")
)

# Prune zero-abundance taxa
ps.target <- prune_taxa(taxa_sums(ps.target) > 0, ps.target)

# Agglomerate to Genus level
ps.genus <- tax_glom(ps.target, taxrank = "Genus")

# Transform to relative abundance
ps.rel <- transform_sample_counts(ps.genus, function(x) x / sum(x))

# Melt to long format
df <- psmelt(ps.rel)

# Get top 10 genera overall
top10 <- df %>%
  group_by(Genus) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# Filter for top 10 only and rename Not_bleached to Healthy
df_top <- df %>%
  filter(Genus %in% top10)%>%
 mutate(Treatment = ifelse(Treatment == "Not_bleached", "Non_bleached", Treatment))

# Summarize and normalize
df_summed <- df_top %>%
  group_by(Treatment, Genus) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Treatment) %>%
  mutate(Relative_Abundance = total_abundance / sum(total_abundance)) %>%
  ungroup()

# Set factor levels for Treatment and Genus (most abundant at bottom)
# First, calculate overall abundance to sort genus globally
genus_order <- df_summed %>%
  group_by(Genus) %>%
  summarise(overall_abundance = sum(Relative_Abundance), .groups = "drop") %>%
  arrange(overall_abundance) %>%
  pull(Genus)

df_summed <- df_summed %>%
  mutate(
    Genus = factor(Genus, levels = genus_order),
    Treatment = factor(Treatment, levels = c("Non_bleached", "Bleached"))
  )

# Plot
ggplot(df_summed, aes(x = Treatment, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_fill_manual(values = palette_colors_50) +
  theme_classic() +
  theme_bw() +
  labs(
    #title = "Top 10 Genera in Bleached vs Non_bleached",
    x = "Treatment",
    y = "Relative Abundance",
    fill = "Genus"
  ) +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold")
  )

### Intertidal 2022 Experiment ###
                                  
# Subset samples
ps.target <- subset_samples(
   ps.cleaned_N,
  Experiment == "Intertidal_2022" &
    Coral_species %in% c("Acropora_tenuis", "Stylophora_pistillata")
)

# Prune zero-abundance taxa
ps.target <- prune_taxa(taxa_sums(ps.target) > 0, ps.target)

# Agglomerate to Genus level
ps.genus <- tax_glom(ps.target, taxrank = "Genus")

# Transform to relative abundance
ps.rel <- transform_sample_counts(ps.genus, function(x) x / sum(x))

# Melt to long format
df <- psmelt(ps.rel)

# Get top 10 genera
top10 <- df %>%
  group_by(Genus) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# Filter top 10
df_top <- df %>%
  filter(Genus %in% top10)

# Summarize and normalize
df_summed <- df_top %>%
  group_by(Treatment, Genus) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Treatment) %>%
  mutate(Relative_Abundance = total_abundance / sum(total_abundance)) %>%
  ungroup()

# Set factor levels
genus_order <- df_summed %>%
  group_by(Genus) %>%
  summarise(overall_abundance = sum(Relative_Abundance), .groups = "drop") %>%
  arrange(overall_abundance) %>%
  pull(Genus)

df_summed <- df_summed %>%
  mutate(
    Genus = factor(Genus, levels = genus_order),
    Treatment = factor(Treatment, levels = unique(Treatment))
  )

# Plot Intertidal 2022
ggplot(df_summed, aes(x = Treatment, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_fill_manual(values = palette_colors_50) +
  labs(title = "Top 10 Genera in Intertidal 2022 Experiment",
       x = "Treatment",
       y = "Relative Abundance",
       fill = "Genus") +
  theme_classic() +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  )

### Nursery vs Natural Reef ###

# Subset samples by Structure
ps.target <- subset_samples(
   ps.cleaned_N,
  Structure %in% c("Artificial_structure", "Natural_reef") &
    Experiment == "Intertidal_2022"
)

# Prune zero-abundance taxa
ps.target <- prune_taxa(taxa_sums(ps.target) > 0, ps.target)

# Agglomerate to Genus level
ps.genus <- tax_glom(ps.target, taxrank = "Genus")

# Transform to relative abundance
ps.rel <- transform_sample_counts(ps.genus, function(x) x / sum(x))

# Melt to long format
df <- psmelt(ps.rel)

# Get top 10 genera
top10 <- df %>%
  group_by(Genus) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# Filter top 10
df_top <- df %>%
  filter(Genus %in% top10)

# Summarize and normalize
df_summed <- df_top %>%
  group_by(Structure, Genus) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Structure) %>%
  mutate(Relative_Abundance = total_abundance / sum(total_abundance)) %>%
  ungroup()

# Set factor levels
genus_order <- df_summed %>%
  group_by(Genus) %>%
  summarise(overall_abundance = sum(Relative_Abundance), .groups = "drop") %>%
  arrange(overall_abundance) %>%
  pull(Genus)

df_summed <- df_summed %>%
  mutate(
    Genus = factor(Genus, levels = genus_order),
    Structure = factor(Structure, levels = c("Artificial_structure", "Natural_reef"))
  )

# Plot Nursery vs Natural Reef
ggplot(df_summed, aes(x = Structure, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_fill_manual(values = palette_colors_50) +
  labs(title = "Top 10 Genera in Nursery vs Natural Reef",
       x = "Structure",
       y = "Relative Abundance",
       fill = "Genus") +
  theme_classic() +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  )

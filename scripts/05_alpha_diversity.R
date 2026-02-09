# ============================================
# 05_alpha_diversity.R
# Purpose: Compute and plot alpha diversity metrics
# ============================================

### 1. Coral Bleaching ###

# Subset samples
ps.bleaching <- subset_samples(ps.cleaned_N,
                               Treatment %in% c("Not_bleached", "Bleached") &
                               Coral_species %in% c("Acropora_tenuis", "Acropora_verweyi"))
ps.bleaching <- prune_taxa(taxa_sums(ps.bleaching) > 0, ps.bleaching)

# Estimate alpha diversity
richness_df <- estimate_richness(ps.bleaching, measures=c("Observed","Shannon","Simpson"))
metadata <- as(sample_data(ps.bleaching), "data.frame")
richness_df <- cbind(richness_df, metadata)

# Rename treatment factor
richness_df$Treatment <- factor(richness_df$Treatment,
                                levels=c("Not_bleached","Bleached"),
                                labels=c("Non_bleached","Bleached"))

# Reshape and plot
richness_long <- richness_df %>%
  pivot_longer(cols=c("Observed","Shannon","Simpson"),
               names_to="Metric", values_to="Value")

ggplot(richness_long, aes(x=Treatment, y=Value, fill=Treatment)) +
  geom_boxplot(width=0.65, outlier.shape=NA, alpha=0.85) +
  stat_boxplot(geom="errorbar", width=0.3) +
  facet_wrap(~Metric, scales="free_y", nrow=1) +
  scale_fill_manual(values=c("Non_bleached"="forestgreen", "Bleached"="grey60")) +
  labs(title="Alpha Diversity Metrics in Coral Bleaching Experiment", x="Treatment", y="Index Value", fill="Treatment") +
  theme_minimal(base_size=16)

### 2. Intertidal 2022 Experiment ###

# Subset samples
ps.alpha <- subset_samples(
  ps.cleaned_N,
  Treatment %in% c("Natural_intertidal", "Intertidal_intertidal", "Intertidal_subtidal", "Natural_subtidal", "Subtidal_subtidal", "Subtidal_intertidal") &
    (Experiment == "Intertidal_2022")
)

# 2. Remove zero-abundance taxa
ps.alpha <- prune_taxa(taxa_sums(ps.alpha) > 0, ps.alpha)

# 3. Estimate alpha diversity
alpha_df <- estimate_richness(ps.alpha, measures = c("Observed", "Shannon", "Simpson"))

# 4. Add metadata back
metadata <- as(sample_data(ps.alpha), "data.frame")

alpha_df <- cbind(alpha_df, metadata)

# 5. Reshape to long format
alpha_long <- alpha_df %>%
  pivot_longer(cols = c("Observed", "Shannon", "Simpson"),
               names_to = "Metric", values_to = "Value")

# 6. Ensure Treatment is a factor with desired order
alpha_long$Treatment <- factor(alpha_long$Treatment, levels = c("Natural_intertidal", "Intertidal_intertidal", "Intertidal_subtidal", "Natural_subtidal", "Subtidal_subtidal", "Subtidal_intertidal"))

# 7. Plot faceted alpha diversity boxplots
ggplot(alpha_long, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9) +
  stat_boxplot(geom = "errorbar", width = 0.3)+
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = palette_colors) +
  labs(
    title = "Alpha Diversity in intertidal 2022 experiment",
    x = "Treatment",
    y = "Alpha Diversity Index",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), 
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

### Nursery vs Natural Reef ####

# Subset samples
ps.alpha <- subset_samples(
  ps.cleaned_N,
  Structure %in% c("Artificial_structure", "Natural_reef") &
    (Experiment == "Intertidal_2022")
)

# 2. Remove zero-abundance taxa
ps.alpha <- prune_taxa(taxa_sums(ps.alpha) > 0, ps.alpha)

# 3. Estimate alpha diversity
alpha_df <- estimate_richness(ps.alpha, measures = c("Observed", "Shannon", "Simpson"))

# 4. Add metadata back
metadata <- as(sample_data(ps.alpha), "data.frame")

alpha_df <- cbind(alpha_df, metadata)

# 5. Reshape to long format

alpha_long <- alpha_df %>%
  pivot_longer(cols = c("Observed", "Shannon", "Simpson"),
               names_to = "Metric", values_to = "Value")

# 6. Rename "Artificial_structure" to "Nursery"
alpha_long$Structure <- recode(alpha_long$Structure,
                               "Artificial_structure" = "Nursery")

# 7. Ensure Structure is a factor with desired order
alpha_long$Structure <- factor(alpha_long$Structure, levels = c("Natural_reef", "Nursery"))

# 8. Plot faceted alpha diversity boxplots
ggplot(alpha_long, aes(x = Structure, y = Value, fill = Structure)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Natural_reef" = "#4CAF50", "Nursery" = "#757575")) +
  labs(
    x = "Structure Type",
    y = "Alpha Diversity Index",
    fill = "Structure"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1), 
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )


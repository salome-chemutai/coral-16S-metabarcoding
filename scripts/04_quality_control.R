# ============================================
# 04_quality_control.R
# Purpose: Rarefaction curves and remove negative control OTUs
# ============================================

# Rarefaction curve
tab <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) tab <- t(tab)
rarecurve(tab, step=100, col=rainbow(nrow(tab)), lwd=2,
          ylab="Observed ASVs", xlab="Number of Reads", main="Rarefaction Curves", label=TRUE)

# Negative control removal
neg_controls <- subset_samples(ps, Sample_ID %in% c("Control_1", "Control_2"))
neg_controls <- prune_taxa(taxa_sums(neg_controls) > 0, neg_controls)
contaminant_taxa <- taxa_names(neg_controls)
ps.cleaned_N <- prune_taxa(!taxa_names(ps) %in% contaminant_taxa, ps)

# Save cleaned phyloseq object
saveRDS(ps.cleaned_N, file="ps.cleaned_N.rds")

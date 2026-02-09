# ============================================
# 02_taxonomy_assignment.R
# Purpose: Assign taxonomy using SILVA v138
# ============================================

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.2_toGenus_trainset.fa.gz",
                       tryRC = TRUE, multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.fa.gz")

# Save taxonomy
saveRDS(taxa, file="taxa.rds")

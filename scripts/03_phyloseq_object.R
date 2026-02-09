# ============================================
# 03_phyloseq_object.R
# Purpose: Create phyloseq object with metadata
# ============================================

# Load metadata
samdf <- as.data.frame(read_xlsx("16S_Samples_Metadata_752025.xlsx"))
rownames(samdf) <- sort(samdf$Sample_code)

# Make sure sequence table rownames match metadata
rownames(seqtab.nochim) <- rownames(samdf)

# Construct phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(samdf),
               tax_table(taxa))

# Save phyloseq object
saveRDS(ps, file="ps.rds")

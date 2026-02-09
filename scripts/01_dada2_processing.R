# ============================================
# 01_dada2_processing.R
# Purpose: Filter raw reads, denoise, merge, remove chimeras
# ============================================

setwd("~/16S_data/")  # Set working directory

# Forward and reverse read files
fnFs <- sort(list.files(pattern="_1_paired.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pattern="_2_paired.fastq.gz", full.names = TRUE))

# Sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)

# Filtered output paths
filt_path <- file.path(getwd(), "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filtering and trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(250,248),
                     maxN=0, maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Dereplication
derep_forward <- derepFastq(filtFs, verbose=TRUE)
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_forward) <- sample.names
names(derep_reverse) <- sample.names

# DADA2 inference
dadaFs <- dada(derep_forward, err=errF, multithread=TRUE)
dadaRs <- dada(derep_reverse, err=errR, multithread=TRUE)

# Merge paired reads
merged_reads <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(merged_reads)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save workspace
save.image("16S_data.RData")

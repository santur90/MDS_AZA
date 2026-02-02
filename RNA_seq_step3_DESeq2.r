# ==============================================================================
# RNA-seq Differential Expression Analysis using DESeq2
# Purpose: Compare gene expression between conditions
# ==============================================================================

# Set working directory
directory <- "RNA_seq/"
setwd(directory)

# ------------------------------------------------------------------------------
# 1. Load and prepare sample metadata
# ------------------------------------------------------------------------------
sampleTable <- read.table("samples.txt", h = T, sep = ",")
colnames(sampleTable) <- c("sampleName", "Treatment", "Type", "batch", "condition", "gender")

# ------------------------------------------------------------------------------
# 2. Load and prepare count data
# ------------------------------------------------------------------------------
counts <- read.table("counts.csv", h = T, row.names = 1, sep = ",", check.names = F)
head(counts)

# Match count matrix columns to sample table
counts <- as.matrix(counts[, colnames(counts) %in% sampleTable$sampleName])

# Verify sample matching
if (!all(colnames(counts) %in% sampleTable$sampleName)) {
  stop("Sample names in count matrix don't match sample table!")
}

# ------------------------------------------------------------------------------
# 3. Define parameters for analysis
# ------------------------------------------------------------------------------
# Comparison groups (verify these match your condition column)
condition_treatment <- "Ctrl"  # Treatment group
condition_control <- "MDS"    # Control/reference group

# Colors for visualization
cols <- c("MDS" = "red", "Ctrl" = "forestgreen")

# Thresholds for significance
qcutoff <- 0.05                    # Adjusted p-value cutoff
logfccut <- 0.58                   # log2 fold change cutoff (2^0.58 ≈ 1.5-fold)

# ------------------------------------------------------------------------------
# 4. Create DESeq2 object and run analysis
# ------------------------------------------------------------------------------
ddsHTSeq <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sampleTable,
  design = ~ batch + condition  # Adjust for batch effect
)

# Run DESeq2 analysis
dds <- DESeq(ddsHTSeq)

# Extract results
res <- results(
  dds, 
  contrast = c("condition", condition_treatment, condition_control),
  cooksCutoff = FALSE,
  alpha = qcutoff
)

# ------------------------------------------------------------------------------
# 5. Summary and QC
# ------------------------------------------------------------------------------
# Print summary statistics
cat("\n=== DESeq2 Results Summary ===\n")
summary(res)

# Check how many significant genes
cat("\nSignificant genes at padj <", qcutoff, ":", sum(res$padj < qcutoff, na.rm = TRUE), "\n")

# ------------------------------------------------------------------------------
# 6. Filter significant genes
# ------------------------------------------------------------------------------
# Remove NA values and filter by thresholds
sig <- as.data.frame(na.omit(res))
sig <- sig[which(sig$padj < qcutoff & abs(sig$log2FoldChange) > logfccut), ]

# Add gene annotations
sig$gene <- rownames(sig)
sig <- sig[order(sig$padj), ]  # Sort by adjusted p-value

cat("\nSignificant genes after fold change filter:", nrow(sig), "\n")
cat("  - Upregulated:", sum(sig$log2FoldChange > 0), "\n")
cat("  - Downregulated:", sum(sig$log2FoldChange < 0), "\n")

# ------------------------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------------------------
# Save all results
write.csv(
  as.data.frame(res), 
  file = "DESeq2_all_results.csv",
  row.names = TRUE
)

# Save significant genes only
write.csv(
  sig, 
  file = "DESeq2_significant_genes.csv",
  row.names = FALSE
)

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(
  normalized_counts, 
  file = "DESeq2_normalized_counts.csv"
)

# ------------------------------------------------------------------------------
# 8. Quality control plots
# ------------------------------------------------------------------------------
pdf("DESeq2_QC_plots.pdf", width = 10, height = 8)

# MA plot
plotMA(res, ylim = c(-5, 5), main = "MA Plot")
abline(h = c(-logfccut, logfccut), col = "red", lty = 2)

# Dispersion plot
plotDispEsts(dds, main = "Dispersion Estimates")

# PCA plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("condition", "batch"))

# Histogram of p-values
hist(res$pvalue[res$baseMean > 1], 
     breaks = 50, 
     col = "grey50", 
     border = "white",
     main = "Distribution of p-values",
     xlab = "p-value")

dev.off()

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:\n")
cat("  - DESeq2_all_results.csv\n")
cat("  - DESeq2_significant_genes.csv\n")
cat("  - DESeq2_normalized_counts.csv\n")
cat("  - DESeq2_QC_plots.pdf\n")

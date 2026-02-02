# =============================================================================
# Differential Methylation Analysis on Training Dataset
# Purpose: Load methylation data and identify DMRs between case and control groups
# =============================================================================

# Load required libraries
library(methylKit)      # Methylation data analysis
library(methylSig)      # Statistical testing for methylation differences
library(gplots)         # Graphics utilities
library(made4)          # Multivariate analysis
library(pheatmap)       # Heatmap visualization
library(chipenrich)     # ChIP-seq enrichment analysis
# =============================================================================
# Load Input Files
# =============================================================================
# Note: Input files were extracted from CpG reports files
# See bsseq::BSseq function documentation for file format details

message("Loading sample phenotype data...")
sample_table <- read.table("pData.txt", stringsAsFactors = FALSE, header = TRUE, 
                          sep = "\t", row.names = 1)

message("Loading methylation position data...")
position_data <- read.table("MethTile_pos.txt", header = TRUE, check.names = FALSE)

message("Loading coverage data...")
coverage_data <- read.table("MethTile_coverage.txt", header = TRUE, check.names = FALSE)

message("Loading methylation ratio data...")
methyl_data <- read.table("MethTile_methratio.txt", header = TRUE, check.names = FALSE)

message("Loading training sample metadata...")
training_pdata <- read.table("Train_samples.txt", header = TRUE, row.names = 1, 
                            check.names = FALSE)

# =============================================================================
# Construct BSseq Object from Training Data
# =============================================================================

message("\nConstructing BSseq object from training samples...")

# Filter coverage and methylation data to match training samples
training_samples <- as.character(rownames(training_pdata))
coverage_training <- coverage_data[, training_samples]
methyl_training <- methyl_data[, training_samples]

message("Training dataset contains ", length(training_samples), " samples")

# Create genomic ranges for methylation sites
genomic_ranges <- GenomicRanges::GRanges(
    seqnames = position_data$seqnames, 
    strand = rep("*", nrow(position_data)),
    ranges = IRanges::IRanges(start = position_data$starts + 1, width = 25)
)

# Construct BSseq object with training data
bsseq_training <- bsseq::BSseq(
    Cov = as.matrix(coverage_training),
    M = as.matrix(methyl_training),
    gr = genomic_ranges
)

# Attach phenotype data to BSseq object
S4Vectors::mcols(bsseq_training) <- training_pdata

message("BSseq object created with ", nrow(bsseq_training), " CpG sites")
# =============================================================================
# Identify Differentially Methylated Regions (DMRs)
# =============================================================================

# Configuration parameters
MIN_COVERAGE <- 3      # Require at least 3 samples per group with sufficient coverage
GROUP_COLUMN <- "Type"
CASE_LABEL <- "1"
CONTROL_LABEL <- "0"
N_CORES <- 5           # Number of cores for parallel processing

# Step 1: Filter loci by coverage per group
# Requirement: at least MIN_COVERAGE samples from each group
message("\nFiltering loci by group coverage...")
message("Requirement: at least ", MIN_COVERAGE, " samples per group")

bsseq_filtered <- methylSig::filter_loci_by_group_coverage(
    bs = bsseq_training,
    group_column = GROUP_COLUMN,
    c(CASE_LABEL = MIN_COVERAGE, CONTROL_LABEL = MIN_COVERAGE)
)
message("Retained ", nrow(bsseq_filtered), " loci after coverage filtering")

# Step 2: Test for differential methylation
# Compare case vs control groups, estimating dispersion from both groups
message("\nTesting for differential methylation...")
message("Comparison: case (", CASE_LABEL, ") vs control (", CONTROL_LABEL, ")")

dmr_results <- methylSig::diff_methylsig(
    bs = bsseq_filtered,
    group_column = GROUP_COLUMN,
    comparison_groups = c("case" = CASE_LABEL, "control" = CONTROL_LABEL),
    disp_groups = c("case" = TRUE, "control" = TRUE),
    local_window_size = 0,
    t_approx = FALSE,
    n_cores = N_CORES
)
message("Identified ", length(dmr_results), " significant DMRs")

# =============================================================================
# Format and Export DMR Results
# =============================================================================

message("\nFormatting DMR results...")

# Format genomic coordinates for DMRs
dmr_coordinates <- data.frame(
    position = paste0(
        GenomicRanges::seqnames(dmr_results), ":",
        start(dmr_results) - 1, ":",
        end(dmr_results)
    )
)

# Combine coordinates with statistical test results
dmr_results_table <- cbind(dmr_coordinates, as.data.frame(S4Vectors::mcols(dmr_results)))

# Create summary dataframe with key statistics
dmr_summary <- data.frame(
    position = dmr_results_table$position,
    meth_diff = dmr_results_table$meth_diff,
    pvalue = dmr_results_table$pvalue,
    qvalue = dmr_results_table$fdr,
    row.names = dmr_results_table$position
)

# Display top DMRs
cat("\nTop 10 differentially methylated regions:\n")
print(head(dmr_summary, 10))

# Export DMR results to file
output_file <- "Train_diff_methTile.values"
message("\nExporting results to: ", output_file)
write.table(
    dmr_summary, 
    file = output_file, 
    quote = FALSE, 
    sep = "\t", 
    row.names = TRUE, 
    col.names = TRUE
)
message("Results exported successfully!")

# Summary statistics
cat("\n=== DMR Analysis Summary ===\n")
cat("Total DMRs identified: ", nrow(dmr_summary), "\n")
cat("Hypermethylated (meth_diff > 0): ", sum(dmr_summary$meth_diff > 0), "\n")
cat("Hypomethylated (meth_diff < 0): ", sum(dmr_summary$meth_diff < 0), "\n")
cat("FDR < 0.05: ", sum(dmr_summary$qvalue < 0.05), "\n")
cat("============================\n")

# =============================================================================
# Methylation Data Analysis Pipeline
# Purpose: Identify differentially methylated regions (DMRs) between case and control
# =============================================================================

# Load required libraries
library(methylKit)      # Methylation data analysis
library(methylSig)      # Statistical testing for methylation differences
library(gplots)         # Graphics utilities
library(made4)          # Multivariate analysis
library(pheatmap)       # Heatmap visualization
library(chipenrich)     # ChIP-seq enrichment analysis
#' Construct BSseq Object from Methylation Data Files
#'
#' @param position_file Character. Path to file containing genomic positions
#' @param coverage_file Character. Path to file containing coverage counts
#' @param methyl_ratio_file Character. Path to file containing methylation ratios
#' @param phenotype_file Character. Path to file containing sample phenotype data
#' @param filter_samples Logical. Whether to filter samples based on phenotype file
#' @param region_width Integer. Width of genomic region for each methylation site (default: 25)
#'
#' @return BSseq object with methylation data and metadata
#'
load_methylation_data <- function(position_file, coverage_file, methyl_ratio_file, 
                                   phenotype_file, filter_samples = TRUE, region_width = 25) {
    # Input validation
    if (!file.exists(position_file)) {
        stop("Position file not found: ", position_file)
    }
    if (!file.exists(coverage_file)) {
        stop("Coverage file not found: ", coverage_file)
    }
    if (!file.exists(methyl_ratio_file)) {
        stop("Methylation ratio file not found: ", methyl_ratio_file)
    }
    if (!file.exists(phenotype_file)) {
        stop("Phenotype file not found: ", phenotype_file)
    }
    
    # Read input files
    message("Reading input files...")
    position_data <- read.table(position_file, header = TRUE, check.names = FALSE)
    coverage_data <- read.table(coverage_file, header = TRUE, check.names = FALSE)
    methyl_data <- read.table(methyl_ratio_file, header = TRUE, check.names = FALSE)
    phenotype_data <- read.table(phenotype_file, header = TRUE, row.names = 1, check.names = FALSE)
    
    # Filter samples to match phenotype data if requested
    if (filter_samples) {
        sample_names <- as.character(rownames(phenotype_data))
        coverage_data <- coverage_data[, sample_names]
        methyl_data <- methyl_data[, sample_names]
        message("Filtered to ", ncol(coverage_data), " samples")
    }

    # Create GRanges object for genomic coordinates
    genomic_ranges <- GenomicRanges::GRanges(
        seqnames = position_data$seqnames, 
        strand = rep("*", nrow(position_data)),
        ranges = IRanges::IRanges(start = position_data$starts + 1, width = region_width)
    )
    
    # Construct BSseq object
    message("Creating BSseq object...")
    bsseq_obj <- bsseq::BSseq(
        Cov = as.matrix(coverage_data),
        M = as.matrix(methyl_data),
        gr = genomic_ranges
    )
    
    # Attach phenotype data
    S4Vectors::mcols(bsseq_obj) <- phenotype_data
    
    # Print data summary
    message("\n=== Data Summary ===")
    message("Samples: ", ncol(bsseq_obj))
    message("CpG sites: ", nrow(bsseq_obj))
    cat("\nFirst few rows of phenotype data:\n")
    print(head(S4Vectors::mcols(bsseq_obj)))
    cat("\nFirst few genomic ranges:\n")
    print(head(GenomicRanges::granges(bsseq_obj)))
    
    return(bsseq_obj)
}
# =============================================================================
# Load Methylation Data
# =============================================================================

# Configuration parameters
PHENOTYPE_FILE <- "pData.txt"
POS_FILE <- "methTile_pos.txt"
COVERAGE_FILE <- "methTile_coverage.txt"
METHYL_RATIO_FILE <- "methTile_methratio.txt"
PHENO_DATA_FILE <- "methTile_pData.txt"

# Read sample phenotype information
sample_table <- read.table(PHENOTYPE_FILE, stringsAsFactors = FALSE, header = TRUE, sep = "\t")

# Load methylation data into BSseq object
methylation_data <- load_methylation_data(
    position_file = POS_FILE,
    coverage_file = COVERAGE_FILE,
    methyl_ratio_file = METHYL_RATIO_FILE,
    phenotype_file = PHENO_DATA_FILE,
    filter_samples = TRUE,
    region_width = 25
)
# =============================================================================
# Identify Differentially Methylated Regions (DMRs)
# =============================================================================

# Filter configuration
MIN_COVERAGE <- 3  # Require at least 3 samples per group with sufficient coverage
GROUP_COLUMN <- "Type"
CASE_LABEL <- "1"
CONTROL_LABEL <- "0"
N_CORES <- 5       # Number of cores for parallel processing

# Filter loci by coverage per group
# Requirement: at least MIN_COVERAGE samples from each group (case and control)
message("\nFiltering loci by group coverage...")
methylation_filtered <- methylSig::filter_loci_by_group_coverage(
    bs = methylation_data,
    group_column = GROUP_COLUMN,
    c(CASE_LABEL = MIN_COVERAGE, CONTROL_LABEL = MIN_COVERAGE)
)
message("Retained ", nrow(methylation_filtered), " loci after coverage filtering")

# Test for differential methylation
# Compare case vs control groups, estimating dispersion from both groups
message("\nTesting for differential methylation...")
dmr_results <- methylSig::diff_methylsig(
    bs = methylation_filtered,
    group_column = GROUP_COLUMN,
    comparison_groups = c("case" = CASE_LABEL, "control" = CONTROL_LABEL),
    disp_groups = c("case" = TRUE, "control" = TRUE),
    local_window_size = 0,
    t_approx = FALSE,
    n_cores = N_CORES
)
message("Identified ", length(dmr_results), " significant DMRs")

# Format DMR coordinates
dmr_coordinates <- data.frame(
    position = paste0(
        GenomicRanges::seqnames(dmr_results), ":",
        start(dmr_results) - 1, ":",
        end(dmr_results)
    )
)

# Combine coordinates with statistical results
dmr_summary <- cbind(dmr_coordinates, as.data.frame(S4Vectors::mcols(dmr_results)))

# Display top DMRs
cat("\nTop 10 differentially methylated regions:\n")
print(head(dmr_summary, 10))

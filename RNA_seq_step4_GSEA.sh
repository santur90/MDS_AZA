#!/bin/bash
set -euo pipefail  # Exit on error, undefined variables, pipe failures

# ==============================================================================
# GSEA Preranked Analysis for RNA-seq Data
# ==============================================================================
# Description: Performs Gene Set Enrichment Analysis on preranked gene lists
# Input: Tab-delimited file with gene IDs and statistical metrics
# Output: GSEA enrichment results
# ==============================================================================

# Configuration
INPUT_FILE="${1:-input.txt}"
OUTPUT_DIR="${2:-./GSEA_results}"
GMT_FILE="${3:-h.all.v2024.1.Hs.symbols.gmt}"  # Example: MSigDB hallmark gene sets
REPORT_LABEL="${4:-GSEA_analysis}"
TOP_PLOTS="${5:-20}"
NPERM="${6:-1000}"
MIN_GENESET_SIZE="${7:-15}"
MAX_GENESET_SIZE="${8:-500}"

# Derived filenames
RANK_FILE="${OUTPUT_DIR}/input.rnk"

# Validation
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file '$INPUT_FILE' not found!" >&2
    exit 1
fi

if [[ ! -f "$GMT_FILE" ]]; then
    echo "ERROR: GMT file '$GMT_FILE' not found!" >&2
    echo "Download from: https://www.gsea-msigdb.org/gsea/downloads.jsp" >&2
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=== Starting GSEA Preranked Analysis ==="
echo "Input file: $INPUT_FILE"
echo "GMT file: $GMT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Step 1: Create ranked gene list
echo "[1/2] Creating ranked gene list..."
sed -n '2,$p' "$INPUT_FILE" | \
    perl -lane '
        @x = split(/\./, $F[0]);  # Remove version from Ensembl IDs (e.g., ENSG00000000003.14)
        print join("\t", $x[0], $F[4]) if defined $F[4] && $F[4] ne "NA";
    ' > "$RANK_FILE"

# Validate rank file
if [[ ! -s "$RANK_FILE" ]]; then
    echo "ERROR: Rank file is empty! Check input file format." >&2
    echo "Expected format: Column 1 = Gene ID, Column 5 = Ranking metric (e.g., log2FC)" >&2
    exit 1
fi

GENE_COUNT=$(wc -l < "$RANK_FILE")
echo "   → Created rank file with $GENE_COUNT genes"

# Step 2: Run GSEA
echo "[2/2] Running GSEA Preranked..."
gsea-cli.sh GSEAPreranked \
    -nperm "$NPERM" \
    -rnk "$RANK_FILE" \
    -rpt_label "$REPORT_LABEL" \
    -plot_top_x "$TOP_PLOTS" \
    -out "$OUTPUT_DIR" \
    -collapse No_Collapse \
    -gmx "$GMT_FILE" \
    -set_min "$MIN_GENESET_SIZE" \
    -set_max "$MAX_GENESET_SIZE" \
    -zip_report false \
    -norm meandiv \
    -scoring_scheme weighted \
    -create_svgs true

echo ""
echo "=== GSEA Analysis Complete ==="
echo "Results saved to: $OUTPUT_DIR"
echo "Summary report: $OUTPUT_DIR/${REPORT_LABEL}.GseaPreranked.*/index.html"

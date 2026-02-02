#!/bin/bash

#==============================================================================
# Script: RNA-seq Step 2 - featureCounts
# Description: Count reads mapped to genomic features using featureCounts
#==============================================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# ============================================================================
# Configuration
# ============================================================================

ANNOTATION="annotation.gtf"
THREADS=5
STRANDEDNESS=0  # 0=unstranded, 1=stranded, 2=reverse-stranded
OUTPUT="all.id.txt"
MIN_QUALITY=x
MIN_OVERLAP=x

# ============================================================================
# Input Validation
# ============================================================================

echo "=========================================="
echo "featureCounts - Read Counting"
echo "=========================================="
echo "Start time: $(date)"
echo ""

# Check if annotation file exists
if [ ! -f "${ANNOTATION}" ]; then
    echo "ERROR: Annotation file '${ANNOTATION}' not found!"
    echo "Please provide a valid GTF annotation file."
    exit 1
fi

# Check if BAM files exist
BAM_COUNT=$(ls -1 *.bam 2>/dev/null | wc -l)
if [ ${BAM_COUNT} -eq 0 ]; then
    echo "ERROR: No BAM files found in current directory!"
    echo "Please ensure aligned BAM files are present."
    exit 1
fi

echo "Found ${BAM_COUNT} BAM file(s) to process:"
ls -1 *.bam
echo ""

# ============================================================================
# Run featureCounts
# ============================================================================

echo "Running featureCounts with parameters:"
echo "  Annotation: ${ANNOTATION}"
echo "  Threads: ${THREADS}"
echo "  Strandedness: ${STRANDEDNESS}"
echo "  Min mapping quality: ${MIN_QUALITY}"
echo "  Min overlap: ${MIN_OVERLAP}"
echo ""

featureCounts \
    -a ${ANNOTATION} \
    -o ${OUTPUT} \
    -p \
    -B \
    -C \
    -T ${THREADS} \
    -s ${STRANDEDNESS} \
    -t exon \
    -g gene_id \
    -Q ${MIN_QUALITY} \
    --minOverlap ${MIN_OVERLAP} \
    --fracOverlap 0 \
    --verbose \
    *.bam \
    2>&1 | tee featureCounts.log

# ============================================================================
# Check Results
# ============================================================================

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "SUCCESS: featureCounts completed!"
    echo "=========================================="
    echo "Output files:"
    echo "  Count matrix: ${OUTPUT}"
    echo "  Summary: ${OUTPUT}.summary"
    echo "  Log file: featureCounts.log"
    echo ""
    echo "End time: $(date)"
    
    # Display summary statistics
    if [ -f "${OUTPUT}.summary" ]; then
        echo ""
        echo "Assignment Summary:"
        echo "-------------------"
        cat ${OUTPUT}.summary
    fi
else
    echo ""
    echo "=========================================="
    echo "ERROR: featureCounts failed!"
    echo "=========================================="
    echo "Please check featureCounts.log for details."
    exit 1
fi

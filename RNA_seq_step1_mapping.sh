#!/bin/bash
set -euo pipefail

# ============================================
# RNA-seq Processing Pipeline
# ============================================

# Configuration parameters
SAMPLE="sample001"
THREADS=18
ADAPTER_CLIP_R1=10
ADAPTER_CLIP_R2=10
MIN_LENGTH=20
QUALITY_CUTOFF=20

# Directory setup
QC_DIR="QC"
TRIM_DIR="trimmed"
RRNA_DIR="rRNA_removed"
ALIGN_DIR="aligned"
LOGS_DIR="logs"

mkdir -p ${QC_DIR} ${TRIM_DIR} ${RRNA_DIR} ${ALIGN_DIR} ${LOGS_DIR}

# Log function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a ${LOGS_DIR}/pipeline.log
}

# Error handling function
error_exit() {
    echo "[ERROR] $1" | tee -a ${LOGS_DIR}/error.log
    exit 1
}

# ============================================
# Step 1: Adapter trimming and quality control
# ============================================
log "Starting adapter trimming for ${SAMPLE}..."
trim_galore \
    -q ${QUALITY_CUTOFF} \
    --phred33 \
    --stringency 6 \
    --length ${MIN_LENGTH} \
    -e 0.2 \
    --fastqc \
    --gzip \
    -o ${TRIM_DIR} \
    --paired \
    --clip_R1 ${ADAPTER_CLIP_R1} \
    --clip_R2 ${ADAPTER_CLIP_R2} \
    ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz \
    2>&1 | tee ${LOGS_DIR}/01_trim_galore.log

[ $? -eq 0 ] || error_exit "Adapter trimming failed"
log "Adapter trimming completed successfully"

# ============================================
# Step 2: Additional quality control (optional)
# ============================================
log "Running FastQC on trimmed reads..."
fastqc ${TRIM_DIR}/${SAMPLE}_R1_val_1.fq.gz \
    -o ${QC_DIR} \
    --extract \
    -t ${THREADS} \
    2>&1 | tee ${LOGS_DIR}/02_fastqc_R1.log

fastqc ${TRIM_DIR}/${SAMPLE}_R2_val_2.fq.gz \
    -o ${QC_DIR} \
    --extract \
    -t ${THREADS} \
    2>&1 | tee ${LOGS_DIR}/02_fastqc_R2.log

[ $? -eq 0 ] || error_exit "FastQC failed"
log "FastQC completed successfully"

# ============================================
# Step 3: rRNA removal
# ============================================
log "Removing rRNA contamination..."
bowtie2 \
    --very-sensitive-local \
    --no-unal \
    -I 1 \
    -X 1000 \
    -p ${THREADS} \
    -x /path/to/rRNA_index \
    -1 ${TRIM_DIR}/${SAMPLE}_R1_val_1.fq.gz \
    -2 ${TRIM_DIR}/${SAMPLE}_R2_val_2.fq.gz \
    --un-conc-gz ${RRNA_DIR}/${SAMPLE}_rRNAremoved.fq.gz \
    2> ${RRNA_DIR}/${SAMPLE}_rRNA_mapping_stats.txt | \
    samtools view -Sb -o ${RRNA_DIR}/${SAMPLE}_mapped_to_rRNA.bam -

[ $? -eq 0 ] || error_exit "rRNA removal failed"
log "rRNA removal completed successfully"

# Calculate rRNA contamination rate
TOTAL_READS=$(grep "reads; of these:" ${RRNA_DIR}/${SAMPLE}_rRNA_mapping_stats.txt | awk '{print $1}')
RRNA_READS=$(grep "aligned concordantly exactly 1 time" ${RRNA_DIR}/${SAMPLE}_rRNA_mapping_stats.txt | awk '{print $1}')
log "rRNA contamination: ${RRNA_READS}/${TOTAL_READS} reads"

# ============================================
# Step 4: Genome alignment with STAR
# ============================================
log "Aligning reads to genome with STAR..."
STAR \
    --runThreadN ${THREADS} \
    --genomeDir /path/to/star_index \
    --readFilesIn ${RRNA_DIR}/${SAMPLE}_rRNAremoved.fq.1.gz \
                  ${RRNA_DIR}/${SAMPLE}_rRNAremoved.fq.2.gz \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFileNamePrefix ${ALIGN_DIR}/${SAMPLE}_ \
    2>&1 | tee ${LOGS_DIR}/04_STAR_alignment.log

[ $? -eq 0 ] || error_exit "STAR alignment failed"
log "STAR alignment completed successfully"

# ============================================
# Step 5: Post-alignment processing
# ============================================
log "Indexing BAM file..."
samtools index ${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam

log "Generating alignment statistics..."
samtools flagstat ${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam \
    > ${ALIGN_DIR}/${SAMPLE}_alignment_stats.txt

samtools idxstats ${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam \
    > ${ALIGN_DIR}/${SAMPLE}_chromosome_stats.txt

[ $? -eq 0 ] || error_exit "Post-alignment processing failed"
log "Post-alignment processing completed successfully"

# ============================================
# Step 6: Gene quantification (optional)
# ============================================
log "Quantifying gene expression with featureCounts..."
featureCounts \
    -T ${THREADS} \
    -p \
    -t exon \
    -g gene_id \
    -a /path/to/annotation.gtf \
    -o ${ALIGN_DIR}/${SAMPLE}_gene_counts.txt \
    ${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam \
    2>&1 | tee ${LOGS_DIR}/06_featureCounts.log

[ $? -eq 0 ] || error_exit "Gene quantification failed"
log "Gene quantification completed successfully"

# ============================================
# Step 7: Generate summary report
# ============================================
log "Generating summary report..."
{
    echo "============================================"
    echo "RNA-seq Pipeline Summary Report"
    echo "Sample: ${SAMPLE}"
    echo "Date: $(date)"
    echo "============================================"
    echo ""
    echo "1. Trimming Statistics:"
    grep "Total reads processed:" ${LOGS_DIR}/01_trim_galore.log || echo "N/A"
    echo ""
    echo "2. rRNA Contamination:"
    cat ${RRNA_DIR}/${SAMPLE}_rRNA_mapping_stats.txt | head -n 5
    echo ""
    echo "3. Alignment Statistics:"
    cat ${ALIGN_DIR}/${SAMPLE}_alignment_stats.txt
    echo ""
    echo "4. Gene Quantification:"
    tail -n 5 ${LOGS_DIR}/06_featureCounts.log
    echo ""
    echo "============================================"
} > ${LOGS_DIR}/${SAMPLE}_summary_report.txt

log "Pipeline completed successfully for ${SAMPLE}!"
log "Summary report: ${LOGS_DIR}/${SAMPLE}_summary_report.txt"

# Clean up intermediate files (optional)
# rm -rf ${TRIM_DIR}/*.fq.gz
# rm -rf ${RRNA_DIR}/*.fq.gz

exit 0

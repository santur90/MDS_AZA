#!/bin/bash

# Settings
set -e # Exit on error, improve debugging
INPUT_FILE="input.fq"
TRIMMED_FILE="input_trimmed.fq"
GENOME_FOLDER="index" # Replace 'index' with genome index path
BUFFER_SIZE="30%"
MIN_LENGTH=x
ADAPTER="NNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

# Step 1: Adapter trimming
echo "Running cutadapt..."
cutadapt -e 0.2 -a $ADAPTER -m $MIN_LENGTH -O 6 $INPUT_FILE > $TRIMMED_FILE 2> cutadapt.log

# Step 2: Bismark mapping
echo "Running bismark mapping..."
bismark --genome $GENOME_FOLDER $TRIMMED_FILE

# Step 3: Methylation extractor
echo "Extracting methylation levels..."
bismark_methylation_extractor --gzip --bedGraph --buffer_size ${BUFFER_SIZE} \
  --cytosine_report --genome_folder $GENOME_FOLDER \
  ${TRIMMED_FILE%%.fq}_bismark_bt2.bam

# Step 4: Generate report and summary
echo "Generating reports..."
bismark2report
bismark2summary

# Done
echo "Script finished!"

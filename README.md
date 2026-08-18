# MDS_AZA

Analysis code supporting the MDS/AZA multi-omics study. The repository contains three linked analysis components:

1. ERRBS methylation profiling and differential methylation analysis
2. RNA-seq alignment, quantification, differential expression, and gene-set enrichment
3. Multi-omics machine-learning classification using methylation, expression, and clinical or mutation features

The scripts are organized by analysis stage and are intended to be run after the corresponding input tables, reference indexes, annotations, and software environments have been configured for the study.

## Repository contents

| Analysis | Main files | Primary outputs |
|---|---|---|
| ERRBS | `ERRBS_step1_bismark.sh`, `ERRBS_step2_methylSig.R` | Bismark reports, methylation tables, DMR results |
| RNA-seq | `RNA_seq_step1_mapping.sh`, `RNA_seq_step2_featureCounts.sh`, `RNA_seq_step3_DESeq2.r`, `RNA_seq_step4_GSEA.sh` | aligned BAMs, gene counts, DESeq2 results, GSEA reports |
| Machine learning | `Github_ML_step1.R`, `Github_ML_step2.R`, `Github_ML_step3.R` | train/test assignments, tuned classifier, test-set metrics, feature importance |

## Analysis flow

```text
ERRBS FASTQ
	-> adapter trimming
	-> Bismark alignment
	-> methylation extraction
	-> methylation tiles and phenotype metadata
	-> methylSig differential methylation

RNA-seq FASTQ
	-> adapter and quality trimming
	-> rRNA filtering
	-> STAR alignment
	-> featureCounts
	-> DESeq2
	-> preranked GSEA

Multi-omics matrices
	-> sample alignment
	-> stratified train/test split
	-> feature selection and classifier tuning
	-> independent test evaluation
```

## Requirements

The analyses require a Linux or HPC environment with the relevant command-line tools and R/Python packages. The main dependencies are:

- ERRBS: `cutadapt`, `Bismark`, `Bowtie2`, `samtools`, `BSseq`, `methylKit`, and `methylSig`
- RNA-seq: `Trim Galore`, `FastQC`, `Bowtie2`, `STAR`, `samtools`, `featureCounts`, `DESeq2`, and a GSEA installation
- Machine learning: Python 3, `numpy`, `pandas`, `scikit-learn`, `scipy`, `matplotlib`, and `mglearn`

Reference indexes and annotation files are not included. Use the same genome assembly for the Bismark/Bowtie2 indexes, STAR index, GTF annotation, and any downstream genomic resources.

## Input data

The supplementary tables define the study-specific inputs. At minimum, prepare:

### ERRBS

- paired-end FASTQ files
- a Bismark genome directory
- methylation position, coverage, and methylation-ratio tables
- phenotype metadata with the comparison column used by `methylSig`

### RNA-seq

- paired-end FASTQ files
- an rRNA Bowtie2 index
- a STAR genome index
- a gene annotation GTF
- sample metadata containing sample name, condition, and batch
- a GMT gene-set file for GSEA

### Machine learning

- phenotype or response-key table indexed by sample ID
- methylation matrix with features as rows and samples as columns
- expression matrix with features as rows and samples as columns
- optional mutation/clinical feature matrix

All matrices must use consistent sample identifiers. The ML scripts export train/test assignments; those assignments should be retained when comparing models or reproducing reported results.

## Running the analyses

Edit paths, sample identifiers, group labels, thresholds, and resource settings at the top of each script before execution. Example commands:

```bash
bash ERRBS_step1_bismark.sh
Rscript ERRBS_step2_methylSig.R

bash RNA_seq_step1_mapping.sh
bash RNA_seq_step2_featureCounts.sh
Rscript RNA_seq_step3_DESeq2.r
bash RNA_seq_step4_GSEA.sh <ranked_input> <gsea_output> <gene_sets.gmt>

python Github_ML_step1.R
python Github_ML_step2.R
python Github_ML_step3.R
```

The files named `Github_ML_step*.R` contain Python code from the original analysis and should be run with Python after renaming them to `.py` or invoking the appropriate Python interpreter. Their input paths and output paths are configured inside the scripts.

## Reproducibility notes

- Record genome assemblies, index versions, annotation releases, and software versions with each analysis run.
- Keep raw sequencing data and large result files outside Git; use the supplementary tables and repository scripts as the reproducible analysis record.
- Review read-level QC, mapping statistics, sample matching, batch effects, multiple-testing correction, and independent test-set performance before interpreting results.
- Machine-learning performance is predictive rather than causal and should be assessed with an external validation cohort whenever available.

## License and attribution

This repository contains project-specific analysis code. Add the study publication citation, data-access statement, and license here when the project metadata are finalized.

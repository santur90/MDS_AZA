# trim adaptor
trim_galore -q 20 --phred33 --stringency 6 --length x -e 0.2 --fastqc input_R1 input_R2 --gzip -o output_dir --paired --clip_R1 y --clip_R2 z
# fastQC
fastqc input_R1 -o QC/ --extract
fastqc input_R2 -o QC/ --extract
# remove rRNA
bowtie2 --very-sensitive-local --no-unal -I 1 -X 1000 -p 6 -x rRNA -1 val_1.fq -2 val_2.fq --un-conc-gz output_rRNAremoved.fq 2>output_Map2rRNAStat.xls | samtools view -S -b -o output_Map2rRNA.bam -
# mapping
STAR --runThreadN 18 --genomeDir star_index --readFilesIn input_rRNAremoved_R1.fq input_rRNAremoved_R2.fq --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix input_rRNAremoved --alignEndsType EndToEnd --readFilesCommand gunzip -c  --outSAMtype BAM SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded --outWigNorm RPM --quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within KeepPairs --alignEndsProtrude a ConcordantPair

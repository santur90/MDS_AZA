# cutadaptor
cutadapt -e 0.2 -a NNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m x -O 6 input.fq > input_trimmed.fq 2> cutadapt.log
# bismark mapping
bismark --genome index input_trimmed.fq
# bismark extractor
bismark_methylation_extractor --gzip --bedGraph --buffer_size 30% --cytosine_report --genome_folder x input_bismark_bt2.bam
# bismark report and summary
bismark2report
bismark2summary

featureCounts -a annotation.gtf -p -T 5 -s x -o all.id.txt -t exon -g gene_id *.bam >featureCounts.log

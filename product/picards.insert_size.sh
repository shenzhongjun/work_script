bam=100.sum_merged.rmdup.bam
/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/software/jre1.8.0_131/bin/java -jar /mnt/share01/tools/picards/picard-tools-1.111/picard.jar CollectInsertSizeMetrics \
I=$bam \
O=$bam.insert_size_metrics.txt \
H=$bam.insert_size_histogram.pdf \
M=0.5


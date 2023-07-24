#!/bin/bash
bam=$1
vcf=$2
out=`pwd`

python /mnt/share02/wangwp/Project/mutvision/mutationVisualize.py \
	--bam $bam \
	--mutlist $vcf \
	--html $out/mutation_visualize.html \
	--ref /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \
	--tlenfile $out/tlen.txt


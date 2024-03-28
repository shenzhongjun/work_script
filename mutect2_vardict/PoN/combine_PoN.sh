#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

cd /mnt/share05/clinical_project/projects/blood_tumor/test/zhouyj_test/GATK+/PoN/every_chip_100+_samples/step3.combine_PoN

mkdir -p ./tmp

/mnt/share02/zhouyj/software/bin/gatk \
    --java-options "-XX:ParallelGCThreads=4 -Xmx50G -Djava.io.tmpdir=./tmp" \
    CreateSomaticPanelOfNormals \
    -R /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \
    --germline-resource /mnt/share02/zhouyj/database/human_b37_reference/somatic-b37_af-only-gnomad.raw.sites.chr.vcf.gz \
    -V gendb:///mnt/share02/zhouyj/database/human_b37_reference/PoN/merge_bed_428samples/genomicsDB \
    -O merge_bed_428samples.PoN.vcf.gz

echo ==== _end_ at `date "+%F  %H:%M:%S"` ====

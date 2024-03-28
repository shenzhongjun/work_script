set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

cd /mnt/share05/clinical_project/projects/blood_tumor/test/zhouyj_test/GATK+/PoN/every_chip_100+_samples/step2.create_genomics_db

mkdir -p ./tmp

/mnt/share02/zhouyj/software/bin/gatk \
    --java-options "-XX:ParallelGCThreads=4 -Xmx50G -Djava.io.tmpdir=./tmp" \
    GenomicsDBImport \
    -R /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \
    -L final.bed \
    --genomicsdb-workspace-path /mnt/share02/zhouyj/database/human_b37_reference/PoN/merge_bed_428samples/genomicsDB \
    --sample-name-map sample.list

echo ==== _end_ at `date "+%F  %H:%M:%S"` ====

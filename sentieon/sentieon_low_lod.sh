#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

OUTDIR=$PWD/test
TUMOR=SJ202215NP69
NORMAL=SJ2022NCP69
TUMOR_FQ1=/mnt/share05/data/product_raw_data/E100042114_L01_DYDFZ-1575-503_1.fq.gz
TUMOR_FQ2=/mnt/share05/data/product_raw_data/E100042114_L01_DYDFZ-1575-503_2.fq.gz
NORMAL_FQ1=/mnt/share05/data/product_raw_data/E100042114_L01_DYDFZ-1575-505_1.fq.gz
NORMAL_FQ2=/mnt/share05/data/product_raw_data/E100042114_L01_DYDFZ-1575-505_2.fq.gz
ADAPTER=/mnt/share02/zhouyj/database/adapter/BGI_adapter.fa
FASTA=/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta
BED=/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NPC69.bed
DBSNP=/mnt/share01/tools/annovar/humandb/mutect2/dbsnp_138.b37.chr.vcf.gz
THREAD=64

mkdir -p $OUTDIR && cd $OUTDIR

export SENTIEON_LICENSE='10.90.1.100:8990'

# ******************************************
# Data QC
# ******************************************
mkdir -p $OUTDIR/1.QC/raw_data && cd $OUTDIR/1.QC/raw_data

ln -sf $TUMOR_FQ1 ./$TUMOR.raw.R1.fastq.gz
ln -sf $TUMOR_FQ2 ./$TUMOR.raw.R2.fastq.gz
ln -sf $NORMAL_FQ1 ./$NORMAL.raw.R1.fastq.gz
ln -sf $NORMAL_FQ2 ./$NORMAL.raw.R2.fastq.gz

mkdir -p $OUTDIR/1.QC/clean_data && cd $OUTDIR/1.QC/clean_data

/mnt/share01/tools/bin/fastp \
    -i $OUTDIR/1.QC/raw_data/$TUMOR.raw.R1.fastq.gz \
    -I $OUTDIR/1.QC/raw_data/$TUMOR.raw.R2.fastq.gz \
    -o $TUMOR.clean.R1.fastq.gz \
    -O $TUMOR.clean.R2.fastq.gz \
    --adapter_fasta $ADAPTER \
    --thread 16 \
    --json $TUMOR.fastp.json \
    --html $TUMOR.fastp.html

/mnt/share01/tools/bin/fastp \
    -i $OUTDIR/1.QC/raw_data/$NORMAL.raw.R1.fastq.gz \
    -I $OUTDIR/1.QC/raw_data/$NORMAL.raw.R2.fastq.gz \
    -o $NORMAL.clean.R1.fastq.gz \
    -O $NORMAL.clean.R2.fastq.gz \
    --adapter_fasta $ADAPTER \
    --thread 16 \
    --json $NORMAL.fastp.json \
    --html $NORMAL.fastp.html

# ******************************************************
# Mapping reads with BWA-MEM
# ******************************************************
mkdir -p $OUTDIR/2.Alignment/bwa && cd $OUTDIR/2.Alignment/bwa

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon bwa mem -R "@RG\\tID:rg_$TUMOR\\tSM:$TUMOR\\tPL:ILLUMINA" \
    -t $THREAD -K 10000000 $FASTA $OUTDIR/1.QC/clean_data/$TUMOR.clean.R1.fastq.gz $OUTDIR/1.QC/clean_data/$TUMOR.clean.R2.fastq.gz | \
    /mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon util sort -o tumor_sorted.bam -t $THREAD --sam2bam -i -

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon bwa mem -R "@RG\\tID:rg_$NORMAL\\tSM:$NORMAL\\tPL:ILLUMINA" \
    -t $THREAD -K 10000000 $FASTA $OUTDIR/1.QC/clean_data/$NORMAL.clean.R1.fastq.gz $OUTDIR/1.QC/clean_data/$NORMAL.clean.R2.fastq.gz | \
    /mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon util sort -o normal_sorted.bam -t $THREAD --sam2bam -i -

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -t $THREAD -i tumor_sorted.bam \
    --algo LocusCollector \
    --fun score_info \
    tumor_score.txt

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -t $THREAD -i tumor_sorted.bam \
    --algo Dedup \
    --score_info tumor_score.txt \
    --metrics tumor_dedup_metrics.txt \
    tumor_deduped.bam

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -t $THREAD -i normal_sorted.bam \
    --algo LocusCollector \
    --fun score_info \
    normal_score.txt

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -t $THREAD -i normal_sorted.bam \
    --algo Dedup \
    --score_info normal_score.txt \
    --metrics normal_dedup_metrics.txt \
    normal_deduped.bam

# ******************************************
# Somatic and Structural variant calling
# ******************************************
mkdir -p $OUTDIR/3.Variation/snv_indel && cd $OUTDIR/3.Variation/snv_indel

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -r $FASTA -t $THREAD \
    --interval $BED \
    --interval_padding 50 \
    -i $OUTDIR/2.Alignment/bwa/tumor_deduped.bam \
    -i $OUTDIR/2.Alignment/bwa/normal_deduped.bam \
    --algo TNscope \
    --tumor_sample $TUMOR --normal_sample $NORMAL \
    --dbsnp $DBSNP \
    --disable_detector sv --trim_soft_clip \
    --min_tumor_allele_frac 0.003 --filter_t_alt_frac 0.003 \
    --max_normal_alt_frac 0.001 --max_normal_alt_qsum 250 --max_normal_alt_cnt 10 \
    --assemble_mode 4 --sv_mask_ext 10 --max_fisher_pv_active 0.05 \
    --resample_depth 100000 \
    $TUMOR.somatic.vcf.gz

python /mnt/share02/zhouyj/script/sentieon/sentieon_vcf_filter.py \
  -i $TUMOR.somatic.vcf.gz \
  -o $TUMOR.somatic.filtered.vcf.gz

perl /mnt/share02/zhangxc/software/annovar/table_annovar.pl \
    $TUMOR.somatic.filtered.vcf.gz \
    /mnt/share01/tools/annovar/humandb/ -buildver hg19 \
    -out $TUMOR.somatic.filtered.anno \
    -remove -protocol refGene,cytoBand,clinvar_20200316,cosmic68,1000g2015aug_all,snp138,esp6500si_all,ljb26_all \
    -operation g,r,f,f,f,f,f,f -vcfinput

echo ==== end at `date "+%F  %H:%M:%S"` ====

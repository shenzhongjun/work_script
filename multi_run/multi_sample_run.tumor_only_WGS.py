#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
WGS文库从fasta原始数据开始，bwa比对，GATK Call变异
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "0.1.0"
__date__ = "2022-7-13 15:39:02"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import os
import subprocess


def get_args():
    parser = argparse.ArgumentParser(description='批量刷出多个样本的shell并投递运行')
    parser.add_argument('--sample_list', '-s', help='样本列表', required=True)
    parser.add_argument('--qsub', '-q', help='是否直接投递，默认否', action='store_true')
    parser.add_argument('--nt', '-t', help='qsub投递线程数', type=int)
    return parser.parse_args()


def make_shell(tumor, fq1, fq2, outdir, sh):
    fasta = "/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta"

    cmd = f"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

mkdir -p {outdir} && cd {outdir}

# ******************************************
# Data QC
# ******************************************
# mkdir -p {outdir}/1.QC/raw_data && cd {outdir}/1.QC/raw_data
# 
# ln -sf {fq1} ./{tumor}.raw.R1.fastq.gz
# ln -sf {fq2} ./{tumor}.raw.R2.fastq.gz
# 
# mkdir -p {outdir}/1.QC/clean_data && cd {outdir}/1.QC/clean_data
# 
# /mnt/share01/tools/bin/fastp \\
#     -i {outdir}/1.QC/raw_data/{tumor}.raw.R1.fastq.gz \\
#     -I {outdir}/1.QC/raw_data/{tumor}.raw.R2.fastq.gz \\
#     -o {tumor}.clean.R1.fastq.gz \\
#     -O {tumor}.clean.R2.fastq.gz \\
#     --adapter_fasta /mnt/share02/zhouyj/database/adapter/BGI_adapter.fa \\
#     --thread 16 --average_qual 20 --length_required 50 \\
#     --json {tumor}.fastp.json \\
#     --html {tumor}.fastp.html


# ******************************************************
# Mapping reads with BWA-MEM
# ******************************************************
# mkdir -p {outdir}/2.Alignment/bwa && cd {outdir}/2.Alignment/bwa && mkdir -p tmp
# 
# /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/software/bwa mem \\
#     -R '@RG\\tID:rg_{tumor}\\tSM:{tumor}\\tPL:ILLUMINA' \\
#     -t 24 -M -B 8 -O 12 -L 15 -U 20 -T 50 -d 75 -r 1.3 -E 4 -k 32 \\
#     -Y {fasta} \\
#     {outdir}/1.QC/clean_data/{tumor}.clean.R1.fastq.gz {outdir}/1.QC/clean_data/{tumor}.clean.R2.fastq.gz \\
#     | /mnt/share02/zhouyj/software/bin/sambamba view -S -f bam -t 16 -l 0 /dev/stdin \\
#     | /mnt/share02/zhouyj/software/bin/sambamba sort -t 24 -m 10G \\
#     --tmpdir=tmp -u -o {tumor}.sorted.bam /dev/stdin
#     
# /mnt/share02/zhouyj/software/bin/sambamba markdup \\
#     -r -t 12 --tmpdir tmp \\
#     --hash-table-size 540000 \\
#     --overflow-list-size 540000 \\
#     --sort-buffer-size 15000 \\
#     --io-buffer-size 1024 \\
#     {tumor}.sorted.bam \\
#     {tumor}.sorted.rmdup.bam \\
#     1>{tumor}.rmdup_stat.txt 2>&1

# ******************************************
# HaplotypeCaller and Annovar
# ******************************************
mkdir -p {outdir}/3.Variation/snv_indel && cd {outdir}/3.Variation/snv_indel && mkdir -p tmp
# 
# /mnt/share02/zhouyj/software/bin/gatk \\
#     --java-options "-XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=tmp" \\
#     HaplotypeCaller \\
#     -R {fasta} \\
#     --min-base-quality-score 30 \\
#     --max-reads-per-alignment-start 0 \\
#     -I {outdir}/2.Alignment/bwa/{tumor}.sorted.rmdup.bam \\
#     -O {tumor}.haplotype.vcf.gz
# 
# perl /mnt/share02/zhouyj/software/annovar/table_annovar.pl \\
#     {tumor}.haplotype.vcf.gz \\
#     /mnt/share02/zhouyj/software/annovar/humandb/ -buildver hg19 \\
#     -out {tumor}.haplotype \\
#     -vcfinput -remove \\
#     -protocol refGene -operation g

# ******************************************
# Mutect2 and Annovar（以25度组做对照，看40度和60度组中是否有新的变异）
# ******************************************
/mnt/share02/zhouyj/software/bin/gatk \\
    --java-options "-XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=tmp" \\
    Mutect2 \\
    -R {fasta} \\
    --germline-resource /mnt/share02/zhouyj/database/human_b37_reference/somatic-b37_af-only-gnomad.raw.sites.chr.vcf.gz \\
    --panel-of-normals /mnt/share02/zhouyj/database/human_b37_reference/somatic-b37_Mutect2-exome-panel.chr.vcf.gz \\
    --max-reads-per-alignment-start 0 \\
    -tumor {tumor} \\
    -normal W12878IDT1 \\
    -I {outdir}/2.Alignment/bwa/{tumor}.sorted.rmdup.bam \\
    -I /mnt/share05/clinical_project/projects/blood_tumor/test/zhouyj_test/CAP/WGS_library/W12878IDT1/2.Alignment/bwa/W12878IDT1.sorted.rmdup.bam \\
    --f1r2-tar-gz {tumor}.f1r2.tar.gz \\
    -O {tumor}.mutect2.raw.vcf.gz

/mnt/share02/zhouyj/software/bin/gatk \\
    --java-options "-XX:ParallelGCThreads=4 -Xmx10G -Djava.io.tmpdir=tmp" \\
    LearnReadOrientationModel \\
    -I {tumor}.f1r2.tar.gz \\
    -O {tumor}.read-orientation-model.tar.gz

/mnt/share02/zhouyj/software/bin/gatk \\
    --java-options "-XX:ParallelGCThreads=4 -Xmx10G -Djava.io.tmpdir=tmp" \\
    FilterMutectCalls \\
    -R {fasta} \\
    --orientation-bias-artifact-priors {tumor}.read-orientation-model.tar.gz \\
    --min-median-mapping-quality 20 \\
    -V {tumor}.mutect2.raw.vcf.gz \\
    -O {tumor}.mutect2.filter.vcf.gz

perl /mnt/share02/zhouyj/software/annovar/table_annovar.pl \\
    {tumor}.mutect2.filter.vcf.gz \\
    /mnt/share02/zhouyj/software/annovar/humandb/ -buildver hg19 \\
    -out {tumor}.mutect2.filter \\
    -vcfinput -remove \\
    -protocol refGene -operation g

echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def do_qsub(sh, nt):
    os.chdir(os.path.dirname(sh))
    if subprocess.call(f'qsub -l p={nt if nt > 1 else nt + 1} {sh}', shell=True) != 0:  # 运行成功返回0
        print(f'#### {sh} 投递失败！####\n')


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)

    with open(sample_list) as f:
        for line in f:
            tumor, fq1, fq2 = line.strip().split('\t')

            if not os.path.exists(f'{wd}/{tumor}'):
                os.mkdir(f'{wd}/{tumor}')
            make_shell(tumor, fq1, fq2, f'{wd}/{tumor}', f'{tumor}.wgs.sh')
            if argv.qsub:
                do_qsub(f'{wd}/{tumor}/{tumor}.wgs.sh', argv.nt)
            else:
                print(f'write {wd}/{tumor}/{tumor}.wgs.sh done.')

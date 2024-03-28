#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
从fasta原始数据开始运行sentieon，结果用annovar注释
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "0.1"
__date__ = "2022-7-8 13:55:57"
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


def make_shell(tumor, normal, platform, t_fq1, t_fq2, n_fq1, n_fq2, outdir, sh):
    if platform == 'T7':
        adapter = '/mnt/share02/zhouyj/database/adapter/BGI_adapter.fa'
        bed = "/mnt/share05/clinical_project/projects/blood_tumor/other/zhouyj_test/quality_assess/somatic_2022/rawdata/BGI/BGI.bed"
    elif platform == 'Illumina':
        adapter = '/mnt/share02/zhouyj/database/adapter/Illumina_adapter.fa'
        bed = '/mnt/share05/clinical_project/projects/blood_tumor/other/zhouyj_test/quality_assess/somatic_2022/rawdata/Illumina/Illumina.bed'
    else:
        raise Exception("平台设置错误，流程仅支持Illumina、T7平台")
    fasta = "/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta"
    dbsnp = "/mnt/share02/zhouyj/database/human_b37_reference/dbsnp_138.b37.chr.vcf.gz"

    cmd = f"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

mkdir -p {outdir} && cd {outdir}

export SENTIEON_LICENSE='10.90.1.100:8990'

# ******************************************
# Data QC
# ******************************************
mkdir -p {outdir}/1.QC/raw_data && cd {outdir}/1.QC/raw_data    # 后续改为python控制路径

ln -sf {t_fq1} ./{tumor}.raw.R1.fastq.gz
ln -sf {t_fq2} ./{tumor}.raw.R2.fastq.gz
ln -sf {n_fq1} ./{normal}.raw.R1.fastq.gz
ln -sf {n_fq2} ./{normal}.raw.R2.fastq.gz

mkdir -p {outdir}/1.QC/clean_data && cd {outdir}/1.QC/clean_data

/mnt/share01/tools/bin/fastp \\
    -i {outdir}/1.QC/raw_data/{tumor}.raw.R1.fastq.gz \\
    -I {outdir}/1.QC/raw_data/{tumor}.raw.R2.fastq.gz \\
    -o {tumor}.clean.R1.fastq.gz \\
    -O {tumor}.clean.R2.fastq.gz \\
    --adapter_fasta {adapter} \\
    --thread 16 \\
    --json {tumor}.fastp.json \\
    --html {tumor}.fastp.html

/mnt/share01/tools/bin/fastp \\
    -i {outdir}/1.QC/raw_data/{normal}.raw.R1.fastq.gz \\
    -I {outdir}/1.QC/raw_data/{normal}.raw.R2.fastq.gz \\
    -o {normal}.clean.R1.fastq.gz \\
    -O {normal}.clean.R2.fastq.gz \\
    --adapter_fasta {adapter} \\
    --thread 16 \\
    --json {normal}.fastp.json \\
    --html {normal}.fastp.html

# ******************************************************
# Mapping reads with BWA-MEM
# ******************************************************
mkdir -p {outdir}/2.Alignment/bwa && cd {outdir}/2.Alignment/bwa

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon bwa mem -R "@RG\\tID:rg_{tumor}\\tSM:{tumor}\\tPL:ILLUMINA" \\
    -t 32 -K 10000000 {fasta} {outdir}/1.QC/clean_data/{tumor}.clean.R1.fastq.gz {outdir}/1.QC/clean_data/{tumor}.clean.R2.fastq.gz | \\
    /mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon util sort -o tumor_sorted.bam -t 32 --sam2bam -i -

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon bwa mem -R "@RG\\tID:rg_{normal}\\tSM:{normal}\\tPL:ILLUMINA" \\
    -t 32 -K 10000000 {fasta} {outdir}/1.QC/clean_data/{normal}.clean.R1.fastq.gz {outdir}/1.QC/clean_data/{normal}.clean.R2.fastq.gz | \\
    /mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon util sort -o normal_sorted.bam -t 32 --sam2bam -i -

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -t 32 -i tumor_sorted.bam \\
    --algo LocusCollector \\
    --fun score_info \\
    tumor_score.txt
    
/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -t 32 -i tumor_sorted.bam \\
    --algo Dedup \\
    --score_info tumor_score.txt \\
    --metrics tumor_dedup_metrics.txt \\
    tumor_deduped.bam

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -t 32 -i normal_sorted.bam \\
    --algo LocusCollector \\
    --fun score_info \\
    normal_score.txt
    
/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -t 32 -i normal_sorted.bam \\
    --algo Dedup \\
    --score_info normal_score.txt \\
    --metrics normal_dedup_metrics.txt \\
    normal_deduped.bam

# ******************************************
# Somatic and Structural variant calling
# ******************************************
mkdir -p {outdir}/3.Variation/snv_indel && cd {outdir}/3.Variation/snv_indel

/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon driver -r {fasta} -t 32 \\
    --interval {bed} \\
    -i {outdir}/2.Alignment/bwa/tumor_deduped.bam \\
    -i {outdir}/2.Alignment/bwa/normal_deduped.bam \\
    --algo TNscope \\
    --tumor_sample {tumor} --normal_sample {normal} \\
    --dbsnp {dbsnp} \\
    --disable_detector sv --trim_soft_clip \\
    --min_tumor_allele_frac 0.005 --filter_t_alt_frac 0.008 \\
    --max_normal_alt_frac 0.001 --max_normal_alt_qsum 250 --max_normal_alt_cnt 10 \\
    --assemble_mode 4 --sv_mask_ext 10 --max_fisher_pv_active 0.05 \\
    --resample_depth 100000 --prune_factor 0 \\
    {tumor}.somatic.pre_filter.vcf.gz

/mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools filter -m + -s "insignificant" -e "(PV>0.25 && PV2>0.25) || (INFO/STR == 1 && PV>0.05)" \\
    {tumor}.somatic.pre_filter.vcf.gz \\
    | /mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools filter -m + -s "orientation_bias" -e "FMT/FOXOG[0] == 1" \\
    | /mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools filter -m + -s "strand_bias" -e "SOR > 3" \\
    | /mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools filter -m + -s "low_qual" -e "QUAL < 20" \\
    | /mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools filter -m + -s "short_tandem_repeat" -e "RPA[0]>=1" \\
    | /mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools filter -m + -s "noisy_region" -e "(ECNT>5 && ECNT<10 && (FMT/AF[0:0]<0.0125 || FMT/AD[0:1]<10)) || (ECNT>=10 && ECNT<20 && (FMT/AF[0:0]<0.015 || FMT/AD[0:1]<12)) || (ECNT>=20 && (FMT/AF[0:0]<0.02 || FMT/AD[0:1]<16))" \\
    | /mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" \\
    | /mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools sort \\
    | /mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04/bin/sentieon util vcfconvert - {tumor}.somatic.filter.vcf.gz

python /mnt/share01/tools/pipeline/clinical_pipeline/dna_pipeline/Releases/v0.9/script/sentieon_vcf_filter.py \\
    -i {tumor}.somatic.filter.vcf.gz \\
    -o {tumor}.somatic.filtered.vcf.gz

perl /mnt/share02/zhangxc/software/annovar/table_annovar.pl \\
    {tumor}.somatic.filtered.vcf.gz \\
    /mnt/share01/tools/annovar/humandb/ -buildver hg19 \\
    -out {tumor}.somatic.filtered.anno \\
    -remove -protocol refGene,cytoBand,clinvar_20200316,cosmic68,1000g2015aug_all,snp138,esp6500si_all,ljb26_all \\
    -operation g,r,f,f,f,f,f,f -vcfinput

echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def do_qsub(sh, nt):
    os.chdir(os.path.dirname(sh))
    if subprocess.call(f'qsub -V -cwd -l p={nt if nt > 1 else nt + 1} {sh}', shell=True) != 0:  # 运行成功返回0
        print(f'#### {sh} 投递失败！####\n')


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)

    with open(sample_list) as f:
        for line in f:
            tumor, normal, platform, t_fq1, t_fq2, n_fq1, n_fq2 = line.strip().split('\t')

            if not os.path.exists(f'{wd}/{tumor}'):
                os.mkdir(f'{wd}/{tumor}')
            make_shell(tumor, normal, platform, t_fq1, t_fq2, n_fq1, n_fq2, f'{wd}/{tumor}', f'{tumor}.sentieon.sh')
            if argv.qsub:
                do_qsub(f'{wd}/{tumor}/{tumor}.sentieon.sh', argv.nt)
            else:
                print(f'write {wd}/{tumor}/{tumor}.sentieon.sh done.')

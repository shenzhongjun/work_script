#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
血肿测试，从bam开始，包括wes、panel、tumor_only
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2021-11-30 15:45:21"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import subprocess
import argparse


def get_args():
    parser = argparse.ArgumentParser(description='批量刷出多个样本的shell并投递运行')
    parser.add_argument('--sample_list', '-s', help='样本列表', required=True)
    parser.add_argument('--qsub', '-q', help='是否直接投递，默认否', action='store_true')
    parser.add_argument('--nt', '-n', help='qsub投递线程数', type=int)
    return parser.parse_args()


def make_wes_shell(tumor, normal, t_bam, n_bam, outdir, sh):
    """
    Sentieon:This section describes TNscope parameters for tissue WES or panel (200-500x depth, AF > 1%)
    """
    cmd = f"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====
 
TUMOR_SM={tumor} #sample name
TUMOR_RGID="rg_$TUMOR_SM" #read group ID
NORMAL_SM={normal} #sample name
NORMAL_RGID="rg_$NORMAL_SM" #read group ID

TUMOR_BAM={t_bam}
NORMAL_BAM={n_bam}

FASTA="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta"
KNOWN_DBSNP="/mnt/share02/zhouyj/database/human_b37_reference/dbsnp_138.b37.chr.vcf.gz"
KNOWN_INDELS="/mnt/share02/zhouyj/database/human_b37_reference/1000G_phase1.indels.b37.chr.vcf.gz"
KNOWN_MILLS="/mnt/share02/zhouyj/database/human_b37_reference/Mills_and_1000G_gold_standard.indels.b37.chr.vcf.gz"
INTERVAL_FILE="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NT01T.bed"

SENTIEON_INSTALL_DIR=/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04
export SENTIEON_LICENSE=10.90.1.100:8990

NT=48
START_DIR="$PWD/TNscope_wes"
BCFTOOLS_BINARY="/mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools"

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
cd $WORKDIR

# ******************************************
# 5. Somatic and Structural variant calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \\
    ${{INTERVAL_FILE:+--interval_padding 50 --interval $INTERVAL_FILE}} \\
    -i $TUMOR_BAM -i $NORMAL_BAM \\
    --algo TNscope --disable_detector sv --trim_soft_clip \\
    --tumor_sample $TUMOR_SM --normal_sample $NORMAL_SM --dbsnp $KNOWN_DBSNP \\
    --assemble_mode 4 --sv_mask_ext 10 --max_fisher_pv_active 0.05 \\
    --min_tumor_allele_frac 0.01 --filter_t_alt_frac 0.01 \\
    --max_normal_alt_frac 0.005 --max_normal_alt_qsum 200 --max_normal_alt_cnt 5 \\
    {tumor}.output_tnscope.pre_filter.vcf.gz

# ******************************************
# 6. Variant filtration
# ******************************************
$BCFTOOLS_BINARY annotate ${{INTERVAL_FILE:+-R $INTERVAL_FILE}} -x "FILTER/triallelic_site" {tumor}.output_tnscope.pre_filter.vcf.gz | \\
    $BCFTOOLS_BINARY filter -m + -s "insignificant" -e "(PV>0.25 && PV2>0.25) || (INFO/STR == 1 && PV>0.05)" | \\
    $BCFTOOLS_BINARY filter -m + -s "orientation_bias" -e "FMT/FOXOG[0] == 1" | \\
    $BCFTOOLS_BINARY filter -m + -s "strand_bias" -e "SOR > 3" | \\
    $BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < 20" | \\
    $BCFTOOLS_BINARY filter -m + -s "short_tandem_repeat" -e "RPA[0]>=10" | \\
    $BCFTOOLS_BINARY filter -m + -s "noisy_region" -e "ECNT>5" | \\
    $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" | \\
    $BCFTOOLS_BINARY norm -f $FASTA -m +any | $BCFTOOLS_BINARY sort |
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - {tumor}.output_tnscope.filtered.vcf.gz

$BCFTOOLS_BINARY filter -i 'FILTER=="PASS"' {tumor}.output_tnscope.filtered.vcf.gz | \\
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - {tumor}.output_tnscope.filtered.PASS.vcf.gz

echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def make_panel_t_n_shell(tumor, normal, t_bam, n_bam, outdir, sh):
    """
    Sentieon:This section describes TNscope parameters for tissue WES or panel (200-500x depth, AF > 1%)
    """
    cmd = f"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

TUMOR_SM={tumor} #sample name
TUMOR_RGID="rg_$TUMOR_SM" #read group ID
NORMAL_SM={normal} #sample name
NORMAL_RGID="rg_$NORMAL_SM" #read group ID

TUMOR_BAM={t_bam}
NORMAL_BAM={n_bam}
FASTA="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta"
KNOWN_DBSNP="/home/zhouyj/development/human_b37_reference/dbsnp_138.b37.chr.vcf.gz"
INTERVAL_FILE="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NLP.bed"

SENTIEON_INSTALL_DIR=/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04
export SENTIEON_LICENSE=10.90.1.100:8990

NT=48
START_DIR="$PWD/TNscope_panel"
BCFTOOLS_BINARY="/mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools"

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
cd $WORKDIR

# ******************************************
# 5. Somatic and Structural variant calling
# ******************************************
# Consider adding `--disable_detector sv --trim_soft_clip` if not interested in SV calling
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \\
    ${{INTERVAL_FILE:+--interval_padding 50 --interval $INTERVAL_FILE}} \\
    -i $TUMOR_BAM -i $NORMAL_BAM \\
    --algo TNsco pe --disable_detector sv --trim_soft_clip \\
    --tumor_sample $TUMOR_SM --normal_sample $NORMAL_SM --dbsnp $KNOWN_DBSNP \\
    --assemble_mode 4 --sv_mask_ext 10 --max_fisher_pv_active 0.05 \\
    --min_tumor_allele_frac 0.01 --filter_t_alt_frac 0.01 \\
    --max_normal_alt_frac 0.005 --max_normal_alt_qsum 200 --max_normal_alt_cnt 5 \\
    {tumor}.output_tnscope.pre_filter.vcf.gz

# ******************************************
# 6. Variant filtration
# ******************************************
$BCFTOOLS_BINARY annotate ${{INTERVAL_FILE:+-R $INTERVAL_FILE}} -x "FILTER/triallelic_site" {tumor}.output_tnscope.pre_filter.vcf.gz | \\
    $BCFTOOLS_BINARY filter -m + -s "insignificant" -e "(PV>0.25 && PV2>0.25) || (INFO/STR == 1 && PV>0.05)" | \\
    $BCFTOOLS_BINARY filter -m + -s "orientation_bias" -e "FMT/FOXOG[0] == 1" | \\
    $BCFTOOLS_BINARY filter -m + -s "strand_bias" -e "SOR > 3" | \\
    $BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < 20" | \\
    $BCFTOOLS_BINARY filter -m + -s "short_tandem_repeat" -e "RPA[0]>=10" | \\
    $BCFTOOLS_BINARY filter -m + -s "noisy_region" -e "(TYPE="snp" && ECNT>5 && ECNT<10 && FMT/AF[0:0]<0.0125 && FMT/AD[0:1]<8) || (TYPE="snp" && ECNT>=10 && ECNT<20 && FMT/AF[0:0]<0.015 && FMT/AD[0:1]<8) || (TYPE="snp" && ECNT>=20 && FMT/AF[0:0]<0.02 && FMT/AD[0:1]<8)" || (TYPE!="snp" && ECNT>5) | \\
    $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -8" | \\
    $BCFTOOLS_BINARY norm -f $FASTA -m +any | $BCFTOOLS_BINARY sort |
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - {tumor}.output_tnscope.filtered.vcf.gz

$BCFTOOLS_BINARY filter -i 'FILTER=="PASS"' {tumor}.output_tnscope.filtered.vcf.gz | \\
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - {tumor}.output_tnscope.filtered.PASS.vcf.gz

echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def make_panel_t_only_shell(tumor, t_bam, outdir, sh):
    """
    Sentieon：This section describes TNscope parameters for ctDNA and other high depth cases (2000-5000x depth, AF > 0.3%, Tumor Only)
    """
    cmd = f"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

TUMOR_SM={tumor} #sample name
TUMOR_RGID="rg_$TUMOR_SM" #read group ID

TUMOR_BAM={t_bam}
FASTA="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta"
KNOWN_DBSNP="/home/zhouyj/development/human_b37_reference/dbsnp_138.b37.chr.vcf.gz"
INTERVAL_FILE="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NLP.bed"

SENTIEON_INSTALL_DIR=/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04
export SENTIEON_LICENSE=10.90.1.100:8990

NT=48
START_DIR="$PWD/TNscope_panel_tumor_only"
BCFTOOLS_BINARY="/mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools"

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
cd $WORKDIR

# ******************************************
# 5. Somatic and Structural variant calling
# ******************************************
# Consider adding `--disable_detector sv --trim_soft_clip` if not interested in SV calling
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \\
    ${{INTERVAL_FILE:+--interval_padding 50 --interval $INTERVAL_FILE}} \\
    -i $TUMOR_BAM \\
    --algo TNscope --disable_detector sv \\
    --tumor_sample $TUMOR_SM --dbsnp $KNOWN_DBSNP \\
    --assemble_mode 4 --disable_detector sv \\
    --min_tumor_allele_frac 0.005 --filter_t_alt_frac 0.007 \\
    --clip_by_minbq 1 --min_init_tumor_lod 3.0 --min_tumor_lod 3.0 --resample_depth 100000 \\
    {tumor}.output_tnscope.pre_filter.vcf.gz

# ******************************************
# 6. Variant filtration
# ******************************************
$BCFTOOLS_BINARY annotate ${{INTERVAL_FILE:+-R $INTERVAL_FILE}} {tumor}.output_tnscope.pre_filter.vcf.gz | \\
    $BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < 10" | \\
    $BCFTOOLS_BINARY filter -m + -s "short_tandem_repeat" -e "RPA[0]>=10" | \\
    $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" | \\
    $BCFTOOLS_BINARY sort | \\
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - {tumor}.output_tnscope.filtered.vcf.gz

$BCFTOOLS_BINARY filter -i 'FILTER=="PASS"' {tumor}.output_tnscope.filtered.vcf.gz | \\
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - {tumor}.output_tnscope.filtered.PASS.vcf.gz

echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def do_qsub(sh, nt):
    os.chdir(os.path.dirname(sh))
    if subprocess.call(f'qsub -V -cwd -l p={nt + 1} {sh}', shell=True) != 0:  # 运行成功返回0
        print(f'#### {sh} 投递失败！####\\n')


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)
    sample_dic = {}
    with open(sample_list) as f:
        f.readline()
        for line in f:
            wes, panel, normal, id = line.strip().split('\t')
            sample_dic[id] = {'wes': wes, 'panel': panel, 'normal': normal, 'id': id}

    for i in sample_dic:
        i = sample_dic[i]
        make_wes_shell(i['wes'], i['normal'],
                       f"{wd}/{i['id']}/bwa_alignment/{i['wes']}/{i['wes']}.final.bam",
                       f"{wd}/{i['id']}/bwa_alignment/{i['normal']}/{i['normal']}.final.bam",
                       f"{wd}/{i['id']}", 'run.wes.sh')
        make_panel_t_n_shell(i['panel'], i['normal'],
                             f"{wd}/{i['id']}/bwa_alignment/{i['panel']}/{i['panel']}.final.bam",
                             f"{wd}/{i['id']}/bwa_alignment/{i['normal']}/{i['normal']}.final.bam",
                             f"{wd}/{i['id']}", 'run.panel.sh')
        make_panel_t_only_shell(i['panel'],
                                f"{wd}/{i['id']}/bwa_alignment/{i['panel']}/{i['panel']}.final.bam",
                                f"{wd}/{i['id']}", 'run.panel.tumor_only.sh')
        for sh in ['run.wes.sh', 'run.panel.sh', 'run.panel.tumor_only.sh']:
            do_qsub(f"{wd}/{i['id']}/{sh}", argv.nt) if argv.qsub else print(f"write {wd}/{i['id']}/{sh} done.")

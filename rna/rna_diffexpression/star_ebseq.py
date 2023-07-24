#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
NRCt0自动刷出run.STAR_EBseq.sh脚本
使用方式：auto账号在项目目录下运行ebseq
"""

cmd = fr"""
control=Whole_Blood
tumor=`cat sample_list.txt |cut -f1|sed -n '2p'`
outdir=`pwd`
fq1=`ls $outdir/clean_data/$tumor/${{tumor}}_1.clean.fastq.gz`
fq2=`ls $outdir/clean_data/$tumor/${{tumor}}_2.clean.fastq.gz`

mkdir -p $outdir/STAR_EBseq/02.DEG
cd $outdir/STAR_EBseq/02.DEG

star_result=`ls $outdir/STAR_Fusion/$tumor/*ReadsPerGene.out.tab`
python /mnt/share02/lixx/05.database/05.GTex/3.test/STAR_EBseq/1.script/star_gtex_merge_EBseq.py \
    --star_ReadsPerGene $star_result \
    --gtex_path /mnt/share02/lixx/05.database/05.GTex/2.GTex_database \
    --tissue_type $control \
    --outdir $outdir/STAR_EBseq/02.DEG \
    --prefix $tumor
"""

with open('run.STAR_EBseq.sh', 'w') as w:
    w.write(cmd)

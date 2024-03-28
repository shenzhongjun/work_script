#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
根据订单号和样本号获得bam文件路径
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2024-3-26 09:47:42"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import subprocess

with open('sample_bam_list.txt') as f:
    for i in f:
        chip, sample, order, bam = i.strip().split('\t')
        if not os.path.exists('step1.call_snp_indel/shell'):
            subprocess.call('mkdir -p step1.call_snp_indel/shell', shell=True)
        with open(f'step1.call_snp_indel/shell/{sample}.call.sh', 'w') as w:
            cmd = f"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

sample={sample}
bam={bam}

cd /mnt/share05/clinical_project/projects/blood_tumor/test/zhouyj_test/GATK+/PoN/every_chip_100+_samples/step1.call_snp_indel

mkdir -p ./tmp

/mnt/share02/zhouyj/software/bin/gatk \\
    --java-options "-XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=./tmp" \\
    Mutect2 \\
    -R /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \\
    -I $bam \\
    -max-mnp-distance 0 \\
    -O $sample.vcf.gz

echo ==== _end_ at `date "+%F  %H:%M:%S"` ====
"""
            w.write(cmd)
        # python调用会把log生成到home目录？
        # subprocess.call(f'cd step1.call_snp_indel/shell && qsub {sample}.call.sh && cd -', shell=True)


#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
umi的bamdst深度查看
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
    parser.add_argument('--sample_list', help='样本列表', required=True)
    parser.add_argument('--qsub', help='是否直接投递，默认否', action='store_true')
    parser.add_argument('--nt', help='qsub投递线程数', type=int)
    return parser.parse_args()


def make_shell(outdir, sh):
    cmd = f"""#!/bin/bash
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

/mnt/share01/tools/bamdst/bamdst-master/bamdst \\
    -p /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NPC69.bed \\
    -o depthQC \\
    --flank 150 \\
    umi_consensus.bam

echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def do_qsub(sh, nt):
    os.chdir(os.path.dirname(sh))
    if subprocess.call(f'qsub -V -cwd -l p={nt} {sh}', shell=True) != 0:  # 运行成功返回0
        print(f'#### {sh} 投递失败！####\\n')


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)
    with open(sample_list) as f:
        for line in f:
            t = line.strip()
            sh = f'bamdst.{t}.sh'
            if not os.path.exists(f'{wd}/{t}/depthQC'):
                os.mkdir(f'{wd}/{t}/depthQC')
            make_shell(f'{wd}/{t}', sh)
            do_qsub(f'{wd}/{t}/{sh}', argv.nt) if argv.qsub else print(f'write {wd}/{t}/{sh} done.')


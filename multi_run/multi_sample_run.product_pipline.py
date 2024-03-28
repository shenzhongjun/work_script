#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个订单的生产流程目录并运行生产流程run.sh
sample_list.txt示例如下：
tumor_id	order_id	sample_type	sample_id	library_id
CAPT1	DDN23024022	T	CAPT1	DYDFZ-2938-502
CAPT1	DDN23024022	N	CAPT2DZ	DYDFZ-2938-506
CAPT3	DDN23024022	T	CAPT3	DYDFZ-2938-503
CAPT3	DDN23024022	N	CAPT4DZ	DYDFZ-2938-505
CAPT5	DDN23024022	T	CAPT5	DYDFZ-2938-504
CAPT5	DDN23024022	N	CAPT6DZ	DYDFZ-2938-507
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2023-9-13 14:34:34"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import subprocess
import argparse
import urllib.parse
import urllib.request
import json
import subprocess


def get_args():
    parser = argparse.ArgumentParser(description='批量刷出多个订单的目录并运行生产流程run.sh')
    parser.add_argument('--sample_list', '-s', help='样本列表', required=True)
    parser.add_argument('--project_id', '-p', help='产品类型，默认Ncet', default='Ncet')
    parser.add_argument('--chip_id', '-c', help='芯片类型，默认IDT', default='IDT')
    parser.add_argument('--owner', '-o', help='病人编号，不给的话根据订单号抓取，较为耗时')
    parser.add_argument('--qsub', '-q', help='刷出sjm脚本后直接投递', action='store_true')
    return parser.parse_args()


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)

    seq_index = '20230911'
    platform = 'T7'
    project_id = argv.project_id               # 'Ncet'
    chip_id = argv.chip_id                     # 'IDT'
    postfix = 'b8213'
    info = '新下机'
    own_id = argv.owner if argv.owner else None
    sample_dict = {}

    with open(sample_list) as f:
        f.readline()
        for line in f:
            tumor_id, order_id, sample_type, sample_id, library_id = line.strip().split('\t')
            fastqa, fastqb = subprocess.getoutput(f'ff {library_id} | grep -v html').strip().split('\n')
            if not own_id:
                request = urllib.request.Request(f"https://lims-api.chigene.cn/api/v1/tumor/orders/{order_id}",
                                                 headers={"token": "zhiyindongfang_zhongliu_2020"})
                response = urllib.request.urlopen(request)
                own_id = json.loads(response.read().decode('utf-8'))["data"]['family_member_code']
            # 研发项目tumor总是唯一的，以此为key
            if tumor_id not in sample_dict:
                sample_dict[tumor_id] = {}
            sample_dict[tumor_id][sample_type] = {
                'order_id': order_id,
                'sample_id': sample_id,
                'library_id': library_id,
                'fastqa': fastqa,
                'fastqb': fastqb,
                'own_id': own_id
            }

    for tumor_id in sample_dict:
        subprocess.call(f"mkdir -p {tumor_id}", shell=True)
        with open(f"{tumor_id}/sample_list.txt", 'w') as w:
            w.write('sample_id\tnormal_sample_id\tindex\tsample_type\tplatform\tfastqa\tfastqb\torder_id\town_id\t'
                    'library_id\tproject_id\tchip_id\tpostfix\tinfo\n')
            w.write(f"{sample_dict[tumor_id]['T']['sample_id']}\t{sample_dict[tumor_id]['N']['sample_id']}\t"
                    f"{seq_index}\tT\t{platform}\t"
                    f"{sample_dict[tumor_id]['T']['fastqa']}\t{sample_dict[tumor_id]['T']['fastqb']}\t"
                    f"{sample_dict[tumor_id]['T']['order_id']}\t{sample_dict[tumor_id]['T']['own_id']}\t"
                    f"{sample_dict[tumor_id]['T']['sample_id']}\t{project_id}\t{chip_id}\t{postfix}\t{info}\n")
            w.write(f"{sample_dict[tumor_id]['N']['sample_id']}\t{sample_dict[tumor_id]['N']['sample_id']}\t"
                    f"{seq_index}\tN\t{platform}\t"
                    f"{sample_dict[tumor_id]['N']['fastqa']}\t{sample_dict[tumor_id]['N']['fastqb']}\t"
                    f"{sample_dict[tumor_id]['T']['order_id']}\t{sample_dict[tumor_id]['N']['own_id']}\t"
                    f"{sample_dict[tumor_id]['N']['sample_id']}\t{project_id}\t{chip_id}\t{postfix}\t{info}\n")
        with open(f"{tumor_id}/run.sh", 'w') as w:
            w.write("""/mnt/share01/tools/miniconda/bin/python /mnt/share01/tools/pipeline/clinical_pipeline/dna_pipeline/Releases/v2.0/clinical_genomic_pipeline_v2.0.py \\
    --samplelist $PWD/sample_list.txt \\
    --project_dir $PWD \\
    --other
""")
        if argv.qsub:
            subprocess.call(f'cd {tumor_id} && sh run.sh && sjm *sjm && cd -', shell=True)

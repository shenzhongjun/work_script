#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
CAP补充实验，获取、整理、输出对应表中的样本信息
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2024-1-29 10:56:56"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description='CAP补充实验，获取、整理、输出对应表中的样本信息')
    parser.add_argument('--mapping', '-m', help='订单、肿瘤、对照对应信息', required=True)
    parser.add_argument('--raw', '-r', help='样本、订单、文库对应信息', required=True)
    parser.add_argument('--out', '-o', help='输出文件路径', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    argv = get_args()
    with open(argv.mapping) as mapping, open(argv.raw) as raw, open(argv.out, 'w') as out:
        out.write('tumor_id	order_id	sample_type	sample_id	library_id\n')
        mapping.readline()
        raw.readline()

        map_dic = {}
        for line in mapping:
            order, old_tumor, normal = line.strip().split('\t')
            map_dic[old_tumor] = {'order': order, 'normal': normal}

        raw_dic = {}
        capdz = ''
        capdz_lib = ''
        for line in raw:
            sample, order, library = line.strip().split('\t')
            if order == '':
                capdz = sample
                capdz_lib = library
            else:
                sample2 = sample.split('-')[0]
                # print(sample2)
                if sample2 in map_dic.keys():
                    raw_dic.setdefault(sample, {}).update(
                        {'tumor_id': sample,
                         'old_tumor': sample2,
                         'order_id': order,
                         'tumor_library': library}
                    )
                    if map_dic[sample2]['normal'] == '无对照':
                        raw_dic.setdefault(sample, {}).update(
                            {'normal_id': capdz, 'normal_library': capdz_lib}
                        )
                    else:
                        raw_dic.setdefault(sample, {}).update(
                            {'normal_id': map_dic[sample2]['normal']}
                        )
                else:
                    # print(sample, sample2)
                    for tumor in raw_dic:
                        if sample2 == raw_dic[tumor]['normal_id']:
                            raw_dic[tumor]['normal_library'] = library

        for sample in raw_dic:
            # print(sample, raw_dic[sample])
            out.write(f'{raw_dic[sample]["tumor_id"]}\t{raw_dic[sample]["order_id"]}\tT\t'
                      f'{raw_dic[sample]["tumor_id"]}\t{raw_dic[sample]["tumor_library"]}\n')
            out.write(f'{raw_dic[sample]["tumor_id"]}\t{raw_dic[sample]["order_id"]}\tN\t'
                      f'{raw_dic[sample]["normal_id"]}\t{raw_dic[sample]["normal_library"]}\n')
        # print(map_dic)
        # print(raw_dic)









#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
CAP补充实验，一个文件夹下所有投递的样本是否投递成功check

"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2024-1-30 14:15:01"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import sys

check_file = sys.argv[1]

with open(check_file) as f, open('failed_samples.txt', 'w') as w:
    samples = []
    for line in f:
        if line.startswith('A') and '/' not in line:
            samples.append(line.strip())
        elif line.startswith('A') and '/' in line:
            samples.remove(line.strip().split('/')[0])
        else:
            pass
    for i in samples:
        w.write(f'{i}\n')





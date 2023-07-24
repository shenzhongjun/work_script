#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
rMATs结果按样本统计可变剪切事件数
"""

import os

class Group(object):
    def __init__(self, name):
        self.name = name
        self.samples = []

    def write_results(self, wd, as_events):
        if not os.path.exists(f'{wd}/{self.name}'): os.mkdir(f'{wd}/{self.name}')
        for sample in self.samples:
            with open(f'{wd}/{self.name}/{sample.name}.可变剪切事件统计结果.txt', 'w') as w:
                w.write('AS事件类型\tAS事件总数\t显著差异事件数\n')
                for event in as_events:
                    for e in sample.as_events:
                        if e.name == event:
                            w.write(f'{event}\t{e.counts}\t{e.sig_counts}\n')


class Sample(object):
    def __init__(self, name, index):
        self.name = name
        self.index = index
        self.as_events = []


class ASEvent(object):
    def __init__(self, name, index):
        self.name = name
        self.index = index
        self.counts = 0
        self.sig_counts = 0

    def get_counts(self, line, group1):
        splits = line.split('\t')
        if self.name == 'MXE':
            fdr = splits[21]
            level = splits[22].split(',')[self.index] if group1 else splits[23].split(',')[self.index]
        else:
            fdr = splits[19]
            level = splits[20].split(',')[self.index] if group1 else splits[21].split(',')[self.index]
        if level != 'NA' and float(level) < 0.98:
            self.counts += 1
            if float(fdr) <= 0.05:
                self.sig_counts += 1


if __name__ == "__main__":
    wd = '/mnt/share_fsd/zhouyj/research/RNA/alternative_splicing_20220819/1.ebseq/ASprofile/rmats/stats'
    as_events = ['A5SS', 'A3SS', 'SE', 'RI', 'MXE']
    group1 = Group('group1')
    group2 = Group('group2')

    with open(f'{wd}/list1') as g1, open(f'{wd}/list2') as g2:
        index1 = index2 = 0
        for s1 in g1:
            sample1 = Sample(s1.strip(), index1)
            for e in as_events:
                event = ASEvent(e, index1)
                sample1.as_events.append(event)
            group1.samples.append(sample1)
            index1 += 1
        for s2 in g2:
            sample2 = Sample(s2.strip(), index2)
            for e in as_events:
                event = ASEvent(e, index2)
                sample2.as_events.append(event)
            group2.samples.append(sample2)
            index2 += 1

    for event in as_events:
        with open(f'{wd}/{event}.MATS.JC.txt') as f:
            f.readline()
            for line in f:
                for sample1 in group1.samples:
                    for e1 in sample1.as_events:
                        if e1.name == event:
                            e1.get_counts(line, True)
                for sample2 in group2.samples:
                    for e2 in sample2.as_events:
                        if e2.name == event:
                            e2.get_counts(line, False)

    group1.write_results(wd, as_events)
    group2.write_results(wd, as_events)

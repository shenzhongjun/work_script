#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/2/22 14:53
# @Author  : kailiu@acorndx.com
# @Descript:
# @File    : mutationVisulize.py
# @Software: PyCharm

import pysam
from argparse import ArgumentParser
from extractReads import getTargetFasta, getMutRelatedReads, statReadsTlen
from visHtml import htmlHeader, htmlCoordinate, htmlMutSite, htmlMutReads, htmlTail, posReference


def parserArgument():
    parser = ArgumentParser()
    parser.add_argument('--bam', action='store', dest='bam', required=True, help='bam file with index')
    parser.add_argument('--mutlist', action='store', dest='mutlist', required=True,
                        help='SNV/Indel Mutation with Chr/start/end/ref/allele')
    parser.add_argument('--html', action='store', dest='html', required=True, help='out html')
    parser.add_argument("--bpQ", action="store", dest="bpQ", default=30, type=int, help="Lowest base quality[30]")
    parser.add_argument("--mapQ", action="store", dest="mapQ", default=40, type=int, help="Lowest mapping quality[40]")
    parser.add_argument("--ref", action="store", dest="ref", required=True, help="reference fasta with index")
    parser.add_argument("--readslen", action="store", dest="readslen", default=150, type=int, help="length of reads[150]")
    parser.add_argument("--tlenfile", action="store", dest="tlenFile", required=True, help="file for recording fragment size")

    params = parser.parse_args()
    return params


# 获取变异列表
def getMutation(mutFile):
    # chr1\t10000\t10000\tG\tA (示例，非真实数据)
    Muts = []
    fileR = open(mutFile)
    for line in fileR:
        data = line.strip().split('\t')
        Muts.append(tuple(data))
    fileR.close()
    return Muts


def main():
    params = parserArgument()
    # length = int(params.readslen)
    length = params.readslen
    bam = params.bam
    ref = params.ref
    html = params.html
    bpQ = params.bpQ
    mapQ = params.mapQ
    tlenFile = params.tlenFile

    tlenFileW = open(tlenFile, 'w')
    title = ['chrom', 'start', 'end', 'ref', 'allele', 'tlenMean', 'tlenMax', 'tlenMin', 'tlenSTD', 'tlenCV', 'tlenMedian', 'alleleReads']
    tlenFileW.write('\t'.join(title) + '\n')
    Bam = pysam.AlignmentFile(bam, 'rb')  # bam 是 0-base 坐标，注意

    # -------------变异展示部分------------#
    Muts = getMutation(params.mutlist)
    htmlHeader(html, bpQ, mapQ)
    for mut in Muts:
        chrom, start, end, refbase, allele = mut[0:5]
        s = int(start) - length - 5
        e = int(start) + length + 5
        fasta = getTargetFasta(ref, '%s:%s-%s' % (chrom, s, e))
        posRef = posReference(s, e)
        htmlCoordinate(html, posRef, fasta)
        htmlMutSite(html, posRef, mut)
        mutReads, refReads, is_equals = getMutRelatedReads(Bam, mut, mapQ)
        passed = htmlMutReads(html, mutReads, refReads, bpQ, posRef, is_equals)
        content = statReadsTlen(mutReads)
        # if passed:
        tlenFileW.write('\t'.join(mut[0:5]) + '\t' + content + '\n')
        # print(mutReads)
    htmlTail(html)
    tlenFileW.close()
    # -------------变异展示部分结束------------#
    Bam.close()


if __name__ == '__main__':
    main()

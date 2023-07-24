#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/2/22 16:02
# @Author  : kailiu@acorndx.com
# @Descript:
# @File    : extractReads.py
# @Software: PyCharm

import pysam
import re
import subprocess
import numpy as np
from collections import Counter


# 获取变异附近的fasta序列
def getTargetFasta(ref, pos):
    cmd = ['samtools', 'faidx', ref, pos]
    status, fa = subprocess.getstatusoutput(' '.join(cmd))
    # status == 0,成功；status==1,错误退出。
    if not status:
        fasta = fa.strip().split('\n')
        fasta = ''.join(fasta[1:])
        # fasta = fasta.upper()
        return fasta
    else:
        return False


# 获得变异的类型
def getMutType(ref, allele):
    refL = len(ref)
    alleleL = len(allele)
    refBase = re.match(r'[ATGC]+', ref)
    alleleBase = re.match(r'[ATGC]+', allele)
    if refL == alleleL and refBase and alleleBase:
        return 'snv'  # 包括 AA>TT
    elif not alleleBase:
        return 'delete'
    else:
        return 'insert'  # 包括 TGATC>GAT 和 AA>TTT(OK)
    # if refL > alleleL:
    #     return 'delete'
    # elif refL < alleleL:
    #     return 'insert'
    # elif refL == alleleL and refBase and not alleleBase:
    #     return 'delete'
    # elif refL == alleleL and not refBase and alleleBase:
    #     return 'insert'
    # elif refL == alleleL and refBase and alleleBase:
    #     return 'snv'  # 包括 AA>TT
    # else:
    #     return False


# 检测两个区域有无交叠
def overLap(CoorA, CoorB):
    centerA = np.average(CoorA)
    centerB = np.average(CoorB)
    distance = abs(centerA - centerB)
    lenA = abs(CoorA[1] - CoorA[0])
    lenB = abs(CoorB[1] - CoorB[0])
    halfLen = np.average([lenA, lenB])
    if halfLen >= distance:  # 交叠
        return True
    else:
        return False


def move(ref_split, delete_seq, real_distances, min_index, nums, index, move_len):
    """
    :param ref_split: 从D开始到下一个D之前的所有reference碱基
    :param delete_seq: 当前这个D的所有缺失碱基
    :param real_distances: 每两个D之间的间隔列表
    :param min_index: 距离矩阵中最小的距离的下标
    :param nums: 保存cigar值的数值部分
    :param index: cigar字母中D的下标
    :param move_len: delete_seq在ref_split上移动的碱基个数
    :return: 在当前移动的计算过后，是否能够进行更大步长的移动
    """
    flag_1 = True  # 当前步长是否还能移动
    ref_1 = ref_split
    n = 0
    if len(delete_seq) >= move_len:
        # 如果缺失的碱基序列长度大于等于移动步长，则可以移动， 比如缺失2个，而步长为1，则可以移动
        # 而缺失2个，步长为3， 则不可以移动
        while flag_1:
            # 删除一个字符匹配
            ref_ = ref_1[n + move_len:]  # 移动步长为n个碱基
            if re.match(delete_seq, ref_):
                # 表示移动过后能够与序列开头匹配，则可以将D移动到这个位置
                n = n + move_len
            else:
                flag_1 = False
        if n > 0:
            # 这两个D之间的距离设为100， 防止再次计算这两个缺失
            real_distances[min_index] = 100
            # 前一个距离加n
            if min_index >= 1:
                real_distances[min_index - 1] = real_distances[min_index - 1] + n
            # 修改移动后的cigar数值
            if index[min_index] - 1 >= 0:
                nums[index[min_index] - 1] = nums[index[min_index] - 1] + n
            nums[index[min_index] + 1] = nums[index[min_index] + 1] - n
            return False
        else:
            # 表示可以进行下一个步长加一的循环
            return True
    else:
        return False


def change_cigar(reference: str, cigar: str, distance: int) -> str:
    """
    :param reference: 参考序列
    :param cigar: cigar值
    :param distance: 需要计算的两个缺失之间的最大距离
    :return: 新的cigar值
    """
    nums = re.findall(r"[\d]+", cigar)  # cigar值中的数字
    letters = re.findall(r"[A-Z=]+", cigar)  # cigar中的字母
    nums = [int(i) for i in nums]  # 转换为数字
    letter_counts = Counter(letters)  # 统计每个字母的个数
    # cigar中没有D， cigar不需要改变
    if "D" not in letter_counts.keys():
        new_cigar = cigar
        return new_cigar
    # 只有一个缺失块，cigar不需要改变
    if letter_counts["D"] == 1:
        new_cigar = cigar
        return new_cigar
    if letter_counts["D"] >= 2:
        print('There is more than 2 deletions.')
        # 该read中有大于2个缺失块
        index = []  # 储存D在letters中的下标
        # 遍历cigar值，找到等于D的下标
        for i, letter in enumerate(letters):
            if letter == "D":
                index.append(i)

        real_distances = []  # 储存两个D之间的距离
        for i, item in enumerate(index[:-1]):
            real_distance = sum(nums[: index[i + 1]]) - 1
            real_distances.append(real_distance)
        real_distances = np.array(real_distances)

        # 从reference上截取D到右临近的M的字符串
        flag = True  # 在距离列表里是否还有小于distance的值
        while flag:
            # 没有小于预设距离的值了
            if all([True if i > distance else False for i in real_distances]):
                flag = False
            # 找到当前距离列表里，距离最小的值的下标
            min_index = np.argmin(real_distances)
            i = index[min_index]  # 当前D的cigar字母或者cigar数字下标
            j = index[min_index + 1]  # 下一个D的cigar字母或者cigar数字下标

            # 2个D之间的计算
            start = sum(nums[0:i])  # 当前D的reference起始位置
            end = start + sum(nums[i: j])  # 下一个D的reference起始位置
            # 两个D之间的reference序列，包含当前的D序列
            ref_split = reference[start: end]
            # 缺失的序列
            delete_seq = reference[start: start + nums[i]]

            flag_loop = True
            n = 1  # 每次移动的步长
            while flag_loop:
                # 以当前步长移动时，是否能够更新cigar值，能的话就更新，不能的话，增加步长
                flag_loop = move(ref_split, delete_seq, real_distances, min_index, nums, index, n)
                n = n + 1
                # 这里最大步长为5， 也即考虑缺失中的重复块最大为5
                if n >= 5:
                    flag_loop = False
        # 新的cigar值
        new_cigar = [str(nums[i]) + letters[i] for i in range(len(nums))]
        new_cigar = "".join(new_cigar)
        # print("reference:", reference)
        # print("old cigar:", cigar)
        # print("new cigar:", new_cigar)
        return new_cigar

# 提取变异相关的reads


def getMutRelatedReads(Bam, mut, mapq):
    mutReads = []
    refReads = []  # 记录的是不支持变异的reads
    chrom, start, end, ref, allele = mut[0:5]
    # print(mut)
    mutType = getMutType(ref, allele)
    print(mut, mutType)

    # 提取变异相关reads
    if mutType == 'insert':
        extractStart = int(start) - 1  # mut 来自annovar注释的1-base坐标，所以需要start-1
        extractEnd = int(end) + 1
    else:
        extractStart = int(start) - 1  # 1-base > 0-base
        extractEnd = int(start) - 1 + len(ref)
    readsEntry = Bam.fetch(chrom, extractStart, extractEnd)

    readsCount = 0

    is_equals = []

    for reads in readsEntry:
        if reads.cigarstring == None:
            continue
        if len(reads.cigarstring) == 0:  # 没比对上的reads。
            continue
        if reads.is_duplicate:  # 去除duplication
            continue
        if reads.mapq < mapq:  # 去除低质量的reads
            continue
        is_equal = True
        readsCount += 1
        addReads = False
        mapPos = reads.pos + 1  # 0-base > 1-base
        sequence = reads.query_sequence
        refFasta = reads.get_reference_sequence().upper()
        reads.cigarstring = change_cigar(refFasta, reads.cigarstring, 7)
        alignment = re.findall(r'\d+\w', reads.cigarstring)
        # queryAStart = reads.query_alignment_start
        queryAStart = reads.pos + 1
        offLen = 0  # 记录碱基在reads上的坐标(正常list的下标)
        offSetRef = 0  # 坐标在ref上的偏移

        # 这里单独处理一下insert问题
        if mutType == 'insert':
            SEQ = {}
            index = queryAStart
            i = 0
            for pattern in alignment:
                #print(pattern,alignment,reads.qname,reads.cigarstring)
                if pattern[-1] == 'M' and pattern[:-1] != '0':
                    for i in range(int(pattern[:-1])):
                        #print(len(sequence),i+offLen)
                        index = queryAStart + i
                        SEQ[index] = sequence[i + offLen]
                    queryAStart = index + 1
                    offLen += i + 1
                elif pattern[-1] == 'D':
                    for i in range(int(pattern[:-1])):
                        index = queryAStart + i
                        SEQ[index] = ''
                    queryAStart = index + 1
                elif pattern[-1] == 'I':
                    for i in range(int(pattern[:-1])):
                        index = queryAStart - 1
                        SEQ[index] += sequence[i + offLen]
                    queryAStart = index + 1
                    offLen += i + 1
                else:
                    index = queryAStart
            # print(sorted(SEQ))
            content = ''
            # print('JIEWEI:',int(end),mapPos+len(SEQ)-1)
            # for i in range(max(int(start),mapPos), min(int(end),queryAStart)+1):
            addReads = False
            for i in range(int(start), int(end) + 1):
                if i not in SEQ:
                    print(i, 'not in ')
                    refReads.append(reads)
                    addReads = True
                    break
                else:
                    # print(i,end,SEQ[i])
                    content += SEQ[i]
            subSEQ = content
            #print(allele,subSEQ)
            if subSEQ[1:] == allele and not addReads:
                # print('marking..........')
                mutReads.append(reads)
                is_equals.append(is_equal)
            else:
                refReads.append(reads)
            continue

        # 接下来mapPos将会一直更新
        for pattern in alignment:
            if pattern[-1] == 'M':
                endPos = mapPos + int(pattern[:-1]) - 1

                if not overLap([mapPos, endPos], [int(start), int(end)]):
                    mapPos = endPos + 1
                    offLen += int(pattern[:-1])  # offLen: 这里指提取的sequence的index，
                    offSetRef += int(pattern[:-1])

                else:
                    if mutType == 'snv':
                        offLen += int(start) - mapPos
                        if sequence[offLen:offLen + len(allele)] == allele:  # 可能会漏掉一些变异在边上的reads
                            # print([mapPos + offLen, endPos], [int(start), int(end)])
                            if [int(start), int(start) + len(allele) - 1] != [int(start), int(end)]:
                                is_equal = False

                            mutReads.append(reads)
                            is_equals.append(is_equal)
                            addReads = True
                            break
                        else:
                            refReads.append(reads)
                            addReads = True
                            break
                    else:
                        mapPos = endPos + 1
                        offLen += int(pattern[:-1])  # 同时调整坐标
                        offSetRef += int(pattern[:-1])  # 同时调整坐标

            elif pattern[-1] == 'D':
                endPos = mapPos + int(pattern[:-1]) - 1

                if not overLap([mapPos, endPos], [int(start), int(end)]):
                    mapPos = endPos + 1
                    offSetRef += int(pattern[:-1])
                else:
                    if mutType == 'delete':
                        if refFasta[offSetRef:offSetRef + len(ref)] == ref and int(pattern[:-1]) >= len(ref):  # 此位置的缺失碱基和变异记录的一样
                            if [mapPos, endPos] != [int(start), int(end)]:
                                is_equal = False
                            # print(offLen,offSetRef,len(ref))
                            # print(refFasta)
                            # print(ref,refFasta[offSetRef:offSetRef+len(ref)],reads.query_name,reads.cigarstring)
                            mutReads.append(reads)
                            is_equals.append(is_equal)
                            addReads = True
                            break
                        else:
                            refReads.append(reads)
                            addReads = True
                            break
                    else:
                        mapPos = endPos + 1
                        offLen += 1

            elif pattern[-1] == 'I':
                # offLen += int(pattern[:-1])
                endPos = mapPos  # 插入是在两个碱基之间插入
                mapPos -= 1

                # print(sequence[offLen:offLen+ len(allele)])
                if not overLap([mapPos, endPos], [int(start), int(end)]):
                    mapPos = endPos  # 这里的处理和M、D不一样，因为插入并没有给比对的序列增加坐标
                    offLen += int(pattern[:-1])

                else:
                    if mutType == 'insert':
                        if sequence[offLen:offLen + len(allele)] == allele:
                            if [mapPos, endPos] != [int(start), int(end)]:
                                is_equal = False
                            mutReads.append(reads)
                            is_equals.append(is_equal)
                            addReads = True
                            break
                        else:
                            refReads.append(reads)
                            addReads = True
                            break
                    else:
                        mapPos = endPos + 1
                        offLen += int(pattern[:-1])

            elif pattern[-1] == 'S':
                offLen += int(pattern[:-1])

            else:
                print('have not record it', pattern)
        if not addReads:  # 最后捡漏
            # print('======')
            refReads.append(reads)
    show = 'supportAllele:{alleleReads}\tsupportRef:{refReads}\ttotalReads:{ReadsCount}\n'.format(
        alleleReads=len(set(mutReads)), refReads=len(set(refReads)), ReadsCount=readsCount)
    #print(show)

    return mutReads, refReads, is_equals


def statReadsTlen(mutReads):  # 统计支持变异的reads的插入片段长度,即fragment的长度
    if len(mutReads) == 0:
        return '0\t' * 7
    Tlen = []
    readsName = set()
    for reads in mutReads:
        if reads.query_name not in readsName:
            if reads.tlen == 0:  # 有些reads因为没有mate-reads所以推算不出tlen值。
                continue
            Tlen.append(abs(reads.tlen))
            readsName.update([reads.query_name])
        else:
            continue
    if len(Tlen) == 0:
        return '0\t' * 7
    AVE = round(np.average(Tlen), 2)
    std = np.std(Tlen)
    STD = round(std, 2)
    MED = np.median(Tlen)
    MAX = np.max(Tlen)
    MIN = np.min(Tlen)
    CV = round(STD / AVE, 2)
    content = '{ave}\t{max}\t{min}\t{std}\t{cv}\t{med}\t{readsCount}'.format(
        ave=AVE, max=MAX, min=MIN, std=STD, cv=CV, med=MED, readsCount=len(mutReads)
    )
    return content

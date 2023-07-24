#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
处理STAR FUSION未必对上的reads进行blastn的json格式结果。前处理步骤如下：
# 从bam中提取没有比对上的reads并转为fasta格式
sambamba view -F "unmapped" Aligned.out.bam|head -100 > unmapped.head.sam
samtools view -Sb unmapped.head.sam > unmapped.head.bam
bedtools bamtofastq -i unmapped.head.bam -fq unmapped.head.fastq
/home/zhouyj/bin/seqkit fq2fa unmapped.head.fastq > unmapped.head.fasta
# 使用unmapped.head.fasta在NCBI网站进行blstn序列比对，比对结果下载为3Z5BZCM4013-Alignment.json。
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023年4月19日15:45:07"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import string
import argparse
import json
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from wordcloud import WordCloud


def get_args():
    parser = argparse.ArgumentParser(description='处理STAR FUSION未必对上的reads进行blastn的json格式结果')
    parser.add_argument('--json', help='输入json文件', required=True)
    parser.add_argument('--out', help='结果输出路径', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    # ---提取关键信息写入txt---
    with open(args.json) as f:
        json_content = json.load(f)

    hit_list = []
    for query in json_content['BlastOutput2']:
        for hit in query['report']['results']['search']['hits']:
            hit_list.append(hit['description'][0])

    word_cloud_string = ''
    with open(f'{args.out}.xls', 'w') as w:
        w.write('blastid\taccession\ttitle\ttaxid\tsciname\n')
        for i in hit_list:
            w.write(f'{i["id"]}\t{i["accession"]}\t{i["title"]}\t{i["taxid"]}\t{i["sciname"]}\n')
            word_cloud_string += i["title"]

    # ---词云效果（From ChatGPT）---

    # 去除标点符号和数字
    text = word_cloud_string.translate(str.maketrans('', '', string.punctuation.replace(',', '')))     # + string.digits
    # 转换为小写字母
    # text = text.lower()
    # 分词并去除停用词
    stop_words = set(stopwords.words('english'))
    # stop_words.update({'program', 'version', 'reference', 'search', 'db', 'params', 'expect', 'sc_match', 'sc_mismatch', 'gap', 'filter', 'results', 'query', 'hit', 'num', 'description', 'id', 'accession', 'taxid', 'title', 'sciname', 'len', 'hsps', 'score', 'evalue', 'identity', 'align', 'seq', 'midline'})
    # 增加本项目停用词
    stop_words.update({'partial', 'sequence', 'chromosome', 'complete', 'genome', 'uncultured', 'gene', 'sequenceuncultured', 'bacterium', })
    tokens = word_tokenize(text)
    filtered_tokens = [token for token in tokens if token not in stop_words]
    text = ' '.join(filtered_tokens)
    # 创建 WordCloud 对象
    wc = WordCloud(width=800, height=400, max_words=100, background_color='white')
    # 生成词云图像
    wc.generate(text)
    # 保存词云图像到文件
    wc.to_file(f'{args.out}.wordcloud.png')

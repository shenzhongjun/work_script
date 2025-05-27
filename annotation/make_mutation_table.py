#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
VEP注释结果提取MANE转录本并添加VCF中必要信息，进一步生成肿瘤突变评级结果。
v0.2 开发版
2024-4-23：提高移码突变优先级
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023-4-11 15:49:46"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import io
import re
import argparse
import logging
import pandas as pd
from pandarallel import pandarallel
import time
import json

logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s %(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
pandarallel.initialize(progress_bar=False, nb_workers=16)       # 使用16个进程加速，与VEP保持一致


def get_args():
    db_path = '/mnt/share02/zhouyj/database/annotation'
    parser = argparse.ArgumentParser(
        description='VEP注释结果提取MANE转录本并添加VCF中必要信息，进一步生成肿瘤突变评级结果。')
    parser.add_argument('--anno', help='VEP生成的txt格式注释文件', required=True)
    parser.add_argument('--vcf', help='用于进行VEP注释的vcf文件', required=True)
    parser.add_argument('--out', help='结果输出路径', required=True)
    parser.add_argument('--mane', help='MANE转录本数据库路径',
                        default=f'{db_path}/custom/MANE/MANE.GRCh38.v1.1.summary.txt')
    parser.add_argument('--omim', help='OMIM数据库路径',
                        default=f'{db_path}/custom/OMIM/omim.txt')
    parser.add_argument('--hotspot', help='癌症热点突变数据库路径',
                        default=f'{db_path}/custom/hotspots_v2.clean.20230609.txt')
    parser.add_argument('--osgene', help='原癌基因抑癌基因数据库路径',
                        default=f'{ db_path}/custom/oncogenes_suppressorgenes.clean.20230609.txt')
    parser.add_argument('--riskgene', help='肿瘤遗传易感基因数据库路径',
                        default=f'{db_path}/custom/慧系版报告-solid_tumor_genetic_predisposed_gene.xlsx')
    parser.add_argument('--oncodb', help='OncoKB和CKB合并后的突变致癌性数据库路径',
                        default=f'{db_path}/custom/cbp_ckb_white.merge.20240315.txt')
    parser.add_argument('--json', help='配置结果输出列名及内容的json文件',
                        default=f'{os.path.dirname(os.path.realpath(__file__))}/make_mutation_table_test.json')
    return parser.parse_args()


def make_anno_loc(x):
    if x['id'].startswith('rs'):
        return x['id'].split(';')[0]
    elif 'TYPE=Complex' in x['info']:
        return f"{x['chr']}_{x['pos']}_{x['ref']}/{x['alt']}"
    elif len(x['ref']) < len(x['alt']):
        return f"{x['chr']}_{int(x['pos']) + 1}_-/{x['alt'][1:]}"
    elif len(x['ref']) > len(x['alt']):
        return f"{x['chr']}_{int(x['pos']) + 1}_{x['ref'][1:]}/-"
    else:
        return f"{x['chr']}_{x['pos']}_{x['ref']}/{x['alt']}"


def whether_somatic(x):
    if 'germline_risk' in x:
        return '突变来源不确定'
    elif 'germline' in x:
        return '胚系突变'
    else:
        return '体细胞突变'


def get_total_depth(x):
    if ',' in x.split(':')[1]:
        ref, alt = x.split(':')[1].split(',')
        return int(ref) + int(alt)
    else:
        return int(x.split(':')[1])


def get_exon_intron(x):
    if x['EXON'] != '-' and x['INTRON'] == '-':
        return x['EXON']
    elif x['EXON'] != '-' and x['INTRON'] != '-':
        if x['HGVSp'] != '-':
            return x['EXON']
        else:
            return x['INTRON']
    elif x['EXON'] == '-' and x['INTRON'] != '-':
        return x['INTRON']
    else:
        return '-'


def is_evidence_sufficient(review_status, mutation_type):
    """
    根据ClinVar审核状态和突变类型判断是否满足证据要求
    """
    # 定义审核状态到星级的映射
    status_to_stars = {
        'practice_guideline': 4,
        'reviewed_by_expert_panel': 3,
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,
        'criteria_provided,_conflicting_interpretations': 1,
        'criteria_provided,_single_submitter': 1,
        'no_assertion_criteria_provided': 0,
        'no_assertion_provided': 0,
        'no_interpretation_for_the_single_variant': 0,
        'no_classification_for_the_single_variant': 0
    }

    # 获取当前状态的星级
    stars = status_to_stars.get(review_status.lower(), 0)   # 尝试获取状态对应的星级，失败时视为0星

    # 根据突变类型判断
    if mutation_type.lower() == 'somatic':
        return stars >= 1  # 体细胞突变需要至少1星
    elif mutation_type.lower() == 'germline':
        # print(review_status, mutation_type, stars)
        return stars >= 2  # 胚系突变需要至少2星
    else:
        return False  # 未知突变类型返回False


def change_consequence(x):
    items = set(x['突变小类'].split(','))

    if not x['HGVSp']:
        for i in ['移码突变', '非移码缺失', '非移码插入']:  # 无p.的移码突变没有意义
            if i in items and '+' not in x['HGVSc'] and '-' not in x['HGVSc']:
                items.remove(i)
        if '延长突变' in items:  # 无p.的延长突变实质为同义突变
            items.remove('延长突变')
            items.add('同义突变')
        if '截断突变' in items:  # 无p.的截断突变实质可能为剪切位点（附近）突变，暂定为未知突变
            items.remove('截断突变')
            items.add('剪切位点附近突变')

    if '编码序列改变' in items or '预测蛋白序列改变' in items:
        items.remove('编码序列改变') if '编码序列改变' in items else items.remove('预测蛋白序列改变')
        if x['HGVSp']:
            if 'fs' in x['HGVSp']:
                items.add('移码突变')
            elif 'delins' in x['HGVSp'] or 'del' in x['HGVSp'] or 'ins' in x['HGVSp']:
                items.add('非移码突变')
            elif 'ext' in x['HGVSp']:
                items.add('延长突变')
            elif 'Ter' in x['HGVSp']:
                items.add('截断突变')
            elif 'Met1?' in x['HGVSp']:
                items.add('起始密码子突变')
            else:
                items.add('未知突变')
        else:
            items.add('剪切位点突变')

    if len(items) == 0:
        items.add('未知突变')

    if '起始密码子保留' in items or '终止密码子保留' in items:
        items.add('同义突变')

    if '非移码缺失' in items or '非移码插入' in items:
        items.add('非移码突变')

    if {'内含子突变', '5`UTR突变', '3`UTR突变'}.intersection(items):
        items.add('非编码区突变')

    if '剪切位点附近突变' in items and x['HGVSp']:
        items.remove('剪切位点附近突变')

    if '起始密码子突变' in items:
        if x['HGVSp'] not in ['p.Met1?', 'p.Val1?', 'p.Ile1?', 'p.Leu1?']:
            items.remove('起始密码子突变')
            if x['HGVSp'] and '?' in x['HGVSp']:
                x['HGVSp'] = x['HGVSp'].replace('?', 'del')
                x['氨基酸改变'] = x['氨基酸改变'].replace('?', 'del')
        elif '同义突变' in items:
            items.remove('起始密码子突变')
            x['HGVSp'] = x['HGVSp'].replace('?', '=')
            x['氨基酸改变'] = x['氨基酸改变'].replace('?', '=')

    if '同义突变' in items and ('移码突变' in items or '非移码突变' in items):
        items.remove('移码突变') if '移码突变' in items else items.remove('非移码突变')

    choose_order = ['同义突变', '截断突变', '延长突变', '移码突变', '起始密码子突变', '非移码突变', '错义突变',
                    '剪切位点突变', '剪切位点附近突变', '非编码区突变', '基因上游突变', '基因下游突变', '基因间区突变',
                    '未知突变']
    chosen_item = ''
    for item in choose_order:
        if item in items:
            chosen_item = item
            break

    x['突变小类'] = ','.join(list(items))  # 突变小类去重
    x['突变类型'] = chosen_item
    return x


def pathogenicity_rating(x, risk_gene_list):
    """
    OM = Oncogenic Moderate; OP = Oncogenic Supporting; OS = Oncogenic Strong; OVS = Oncogenic Very Strong;
    SBP = Somatic Benign Supporting; SBS = Somatic Benign Strong; SBVS = Somatic Benign Very Strong
    经与医学商讨，对遗传突变关注其致癌性而非致病性，所以也应用本SOP评级
    """
    # if x['突变来源'] != '胚系突变':
    evidence = []
    # 人群频率评分
    if float(x['MAX_AF']) >= 0.05:
        points = -8
        evidence.append('SBVS1')
    elif 0.01 <= float(x['MAX_AF']) < 0.05:
        points = -4
        evidence.append('SBS1')
    else:
        points = 1
        evidence.append('OP4')

    # 功能数据评分
    """
    OS2判断条件：OS1不适用时，
    OncoKB字段包含'Oncogenic'、'of function'，
    或ClinVar_ONC字段包含'Oncogenic'且ClinVar_ONCREVSTAT证据有效，
    或基因在遗传易感基因列表且ClinVar_CLNSIG字段包含'Pathogenic'、'Likely_pathogenic'且证据有效
    """
    if 'OS1' not in evidence and (
            ('Oncogenic' in x['OncoKB'] or 'of function' in x['OncoKB']) or
            ('Oncogenic' in x['ClinVar_ONC'] and is_evidence_sufficient(x['ClinVar_ONCREVSTAT'], 'somatic')) or
            (x['SYMBOL'] in risk_gene_list and ('Pathogenic' in x['ClinVar_CLNSIG'] or 'Likely_pathogenic' in x['ClinVar_CLNSIG']) and
              is_evidence_sufficient(x['ClinVar_CLNREVSTAT'], 'germline'))
    ):
        points += 4
        evidence.append('OS2')

    if ('Neutral' in x['OncoKB'] or 'no effect' in x['OncoKB']) or (
            ('Benign' in x['ClinVar_CLNSIG'] or 'Likely_benign' in x['ClinVar_CLNSIG']) and
            (is_evidence_sufficient(x['ClinVar_CLNREVSTAT'], 'germline')) and (x['SYMBOL'] in risk_gene_list)
    ):
        points += -4
        evidence.append('SBS2')

    # OM1蛋白功能域，去除重复区域、间区；2025年3月26日：两个数据库改为取交集
    ignore_keywords = ['重复区域', 'Repeat', 'repeat', 'Spacer', 'spacer']
    if ('OS1' not in evidence and 'OS3' not in evidence) and (
            (x['UniProt'] != '-' and all(a not in x['UniProt'] for a in ignore_keywords)) and
            (x['Interpro_domain'] != '-' and all(a not in x['Interpro_domain'] for a in ignore_keywords))):
        points += 2
        evidence.append('OM1')

    # 预测数据
    aaposion = -9 if x['氨基酸改变'] == '-' else re.match(r'p\..*([0-9]+).*', x['氨基酸改变']).group(1)
    mutation_ratio = float((int(x['Protein_len']) - int(aaposion)) / int(x['Protein_len']))
    if x['抑癌基因'] == '是' and (
            (x['突变类型'] == '剪切位点突变') or (
            (x['突变类型'] == '截断突变' or x['突变类型'] == '移码突变') and mutation_ratio >= 0.1)
    ):
        points += 8
        evidence.append('OVS1')
    if (
            x['抑癌基因'] == '是' and (
            (x['突变类型'] == '非移码突变') or (
            (x['突变类型'] == '截断突变' or x['突变类型'] == '移码突变') and mutation_ratio < 0.1))
    ) or (
            x['原癌基因'] == '是' and x['突变类型'] == '非移码突变'
    ):
        points += 2
        evidence.append('OM2')
    if x['突变类型'] == '同义突变' and '剪切位点' not in x['突变小类']:
        points += -1
        evidence.append('SBP2')

    # 癌症热点
    mut_count = 0 if x['Mutation_Count'] == '-' else int(x['Mutation_Count'])
    aa_count = 0 if x['AA_Change_Count'] == '-' else int(x['AA_Change_Count'])
    if mut_count >= 50 and aa_count >= 10:
        points += 4
        evidence.append('OS3')
    elif mut_count < 50 and aa_count >= 10:
        points += 2
        evidence.append('OM4')
    elif 0 < aa_count < 10:
        points += 1
        evidence.append('OP3')

    # 软件证据
    if x['CADD预测'] == x['FATHMM预测'] == 'D':
        points += 1
        evidence.append('OP1')
    elif x['CADD预测'] == x['FATHMM预测'] == 'B':
        points += -1
        evidence.append('SBP1')

    # 输出结论
    if points >= 10:
        x['pathogenicity'] = '致癌'
    elif 6 <= points < 10:
        x['pathogenicity'] = '可能致癌'
    elif 0 <= points < 6:
        x['pathogenicity'] = '意义不明'
    elif -6 <= points < 0:
        x['pathogenicity'] = '可能良性'
    else:
        x['pathogenicity'] = '良性'
    x['points'] = points
    x['evidence'] = ','.join(evidence)
    # else:
    #     x['pathogenicity'] = '-'
    #     x['points'] = '-'
    #     x['evidence'] = '-'
    return x


def trim_interpro(x):
    if x != '-':
        split_set = set(x.split(','))
        if '.' in split_set: split_set.remove('.')
        if len(split_set) >= 1:
            return ','.join(list(split_set))
        else:
            return '-'
    else:
        return '-'


def get_max_af(x):
    afs = ["AF", "EAS_AF", "AFR_AF", "AMR_AF", "EUR_AF", "SAS_AF", "gnomAD_exomes_controls_AF",
           "gnomAD_exomes_controls_EAS_AF", "gnomAD_exomes_controls_AFR_AF", "gnomAD_exomes_controls_AMR_AF",
           "gnomAD_exomes_controls_NFE_AF", "gnomAD_exomes_controls_SAS_AF",
           "WBBC_AF", "WBBC_North_AF", "WBBC_Central_AF", "WBBC_South_AF", "WBBC_Lingnan_AF"]
    values = []
    for i in afs:
        if x[i] == '-':
            values.append(0)
        elif 'gnomADe' in i and x['gnomADe_flag'] != '.':
            values.append(0)
        else:
            values.append(round(float(x[i]), 4))
    return max(values)


def get_alt_depth(x):
    if ',' in x.split(':')[1]:
        return int(x.split(':')[1].split(',')[1])
    else:
        return int(x.split(':')[2])


def get_soft_pred(x, soft):  # 数值越大越有害的软件
    threshold = []  # 判断有害无害的上下限，来源于文献https://www.sohu.com/a/534402602_121124574
    if soft == 'cadd':
        threshold = [22.7, 25.3]
    elif soft == 'vest4':
        threshold = [0.449, 0.764]
    elif soft == 'polyphen':
        threshold = [0.113, 0.978]
    elif soft == 'revel':
        threshold = [0.290, 0.664]
    elif soft == 'fathmm-xf':
        threshold = [0.5, 0.500001]

    if x == '-' or x == '.':
        return '-'
    elif float(x) >= threshold[1]:
        return 'D'  # 有害 f'{x},D'
    elif float(x) <= threshold[0]:
        return 'B'  # 良性 f'{x},B'
    else:
        return '-'


def get_soft2_pred(x, soft):  # 数值越小越有害的软件。2025-5-12注：fathmm已弃用，改成fathmm-xf
    threshold = []
    if soft == 'sift':
        threshold = [0, 0.08]
    elif soft == 'fathmm':
        threshold = [-4.14, 3.32]

    if x == '-' or x == '.':
        return '-'
    elif float(x) < threshold[0]:
        return 'D'  # 有害
    elif float(x) > threshold[1]:
        return 'B'  # 良性
    else:
        return '-'


def get_dbscSNV_pred(x):
    if x['ada_score'] == '-' or x['rf_score'] == '-':
        return '-'
    elif float(x['ada_score']) < 0.6 and float(x['rf_score']) < 0.6:
        return 'B'  # 良性 f"{x['ada_score']},{x['rf_score']},B"
    else:
        return 'D'  # 有害   f"{x['ada_score']},{x['rf_score']},D"


def get_mane_pred(x):
    if x['HGVSc'] in x['HGVSc_VEP'].split(','):
        mane_index = x['HGVSc_VEP'].split(',').index(x['HGVSc'])
        x['REVEL_score'] = x['REVEL_score'].split(',')[mane_index]
        x['Polyphen2_HDIV_score'] = x['Polyphen2_HDIV_score'].split(',')[mane_index]
        x['SIFT_score'] = x['SIFT_score'].split(',')[mane_index]
        # x['fathmm-XF_coding_score'] = x['fathmm-XF_coding_score'].split(',')[mane_index]
        x['VEST4_score'] = x['VEST4_score'].split(',')[mane_index]
    else:
        x['REVEL_score'] = '-'
        x['Polyphen2_HDIV_score'] = '-'
        x['SIFT_score'] = '-'
        # x['fathmm-XF_coding_score'] = '-'
        x['VEST4_score'] = '-'
    return x


if __name__ == "__main__":
    logging.info('开始运行：以VEP注释结果为基础生成Merge表')
    start_time = time.time()
    args = get_args()

    # =====================数据准备=====================
    logging.info('提取MANE转录本并与注释前vcf合并')
    with open(args.json) as f:
        json_content = json.load(f)
    # 读取注释前的vcf
    vcf_df = pd.read_table(
        args.vcf,
        comment='#',
        names=['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'tumor', 'normal'],
        dtype='str'
    )
    vcf_df['Uploaded_variation'] = vcf_df.apply(make_anno_loc, axis=1)
    # 读取VEP注释结果
    with open(args.anno) as f:
        anno = '\n'.join(line for line in f if not line.startswith('##'))   # 删除所有表头
        anno = re.sub('#Uploaded_variation', 'Uploaded_variation', anno)
    anno_df = pd.read_table(io.StringIO(anno), dtype='str')
    # 读取MANE转录本数据库
    mane_df = pd.read_table(args.mane, dtype='str')
    mane_df = mane_df[mane_df['MANE_status'] == 'MANE Select']
    mane_df['RefSeq_id'] = mane_df['RefSeq_nuc'].str.split('.', expand=True)[0]  # 截取转录本id，去掉版本号
    mane_df = mane_df[['RefSeq_id', 'Protein_len']]
    # 提取MANE转录本
    anno_df['RefSeq_id'] = anno_df['Feature'].str.split('.', expand=True)[0]
    anno_df = anno_df.merge(mane_df, on='RefSeq_id', how='left')
    # 针对脑瘤TERT启动子突变chr5_1295228_G/A和chr5_1295250_G/A注释到基因间区没有转录本和基因的情况单独设置
    # 注：VEP --distance 改成150也可以解决此问题，但可能带来其它假阳，所以直接在代码里改
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295228_G/A', 'SYMBOL'] = 'TERT'
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295228_G/A', 'Gene'] = '7015'
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295228_G/A', 'Feature'] = 'NM_198253'
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295228_G/A', 'HGVSc'] = 'NM_198253:c.-124G>A'
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295228_G/A', 'Protein_len'] = '1132'
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295250_G/A', 'SYMBOL'] = 'TERT'
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295228_G/A', 'Gene'] = '7015'
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295250_G/A', 'Feature'] = 'NM_198253'
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295250_G/A', 'HGVSc'] = 'NM_198253:c.-146G>A'
    anno_df.loc[anno_df['Uploaded_variation'] == 'chr5_1295250_G/A', 'Protein_len'] = '1132'
    anno_df = anno_df.dropna(subset=['Protein_len'])        # Protein_len为NA代表位于基因间区
    # 与注释前vcf合并以获取vcf信息
    merge_df = anno_df.merge(vcf_df, on='Uploaded_variation', how='left')

    # =====================数据整理=====================
    # -----------基本信息-----------
    logging.info('变异基本信息整理')
    merge_df['总深度_肿瘤'] = merge_df['tumor'].apply(get_total_depth)
    merge_df['突变深度_肿瘤'] = merge_df['tumor'].apply(get_alt_depth)
    merge_df['突变率_肿瘤'] = (merge_df['突变深度_肿瘤'] / merge_df['总深度_肿瘤']).fillna(0).round(4)
    merge_df['总深度_对照'] = merge_df['normal'].apply(get_total_depth)
    merge_df['突变深度_对照'] = merge_df['normal'].apply(get_alt_depth)
    merge_df['突变率_对照'] = (merge_df['突变深度_对照'] / merge_df['总深度_对照']).fillna(0).round(4)
    merge_df['突变可靠性'] = merge_df['filter'].apply(lambda x: '相对可靠' if 'PASS' in x else '不可靠')
    merge_df['突变来源'] = merge_df['filter'].apply(whether_somatic)
    merge_df['突变小类'] = merge_df['Consequence'].replace(json_content['cons_replace_dict'], regex=True)

    # -----------注释信息-----------
    logging.info('变异注释信息整理')
    merge_df['起始位置'] = merge_df['Location'].str.split(':', expand=True)[1].str.split('-', expand=True)[0]
    merge_df['终止位置'] = merge_df['Location'].apply(
        lambda x: x.split(':')[1].split('-')[1] if '-' in x.split(':')[1] else x.split(':')[1])
    merge_df['HGVSc'] = merge_df['HGVSc'].str.split(':', expand=True)[1]
    merge_df['HGVSp'] = merge_df['HGVSp'].str.split(':', expand=True)[1]
    merge_df['氨基酸改变'] = merge_df['HGVSp'].replace(json_content['amino_acids_dict'], regex=True)
    merge_df['STRAND'] = merge_df['STRAND'].replace(json_content['strand_replace_dict'])
    merge_df['EXON'] = merge_df['EXON'].apply(lambda x: 'exon' + x.replace('/', ',') if x != '-' else x)
    merge_df['INTRON'] = merge_df['INTRON'].apply(lambda x: 'IVS' + x.replace('/', ',') if x != '-' else x)
    merge_df['突变位置'] = merge_df.apply(get_exon_intron, axis=1)
    merge_df = merge_df.parallel_apply(change_consequence, axis=1)      # 结合注释信息，突变小类整合为突变类型

    # -----------公开数据库-----------
    logging.info('公开数据库注释信息整理')
    merge_df['rs号'] = merge_df['Existing_variation'].apply(lambda x: x.split(',')[0] if 'rs' in x else '-')
    # OMIM数据库，用geneid匹配
    logging.info('公开数据库注释信息整理 - 匹配OMIM数据库')
    omim_df = pd.read_table(args.omim, dtype='str', na_filter=False)
    omim_df = omim_df[['基因id', '疾病omim', '遗传方式', '疾病表型']].rename(columns={'基因id': 'Gene'})
    merge_df = merge_df.merge(omim_df, on='Gene', how='left')
    # Cancer Hotspot数据库，用gene+AAchange匹配
    logging.info('公开数据库注释信息整理 - 匹配Cancer Hotspot数据库')
    hotspot_df = pd.read_table(args.hotspot, dtype='str')
    hotspot_df = hotspot_df[['Hugo_Symbol', 'Mutation_Count', 'AA_Change', 'AA_Change_Count']]
    hotspot_df['hotspot_id'] = hotspot_df['Hugo_Symbol'].str.cat(hotspot_df['AA_Change'], sep=':')
    merge_df['hotspot_id'] = merge_df['SYMBOL'].str.cat(merge_df['氨基酸改变'], sep=':')
    merge_df = merge_df.merge(hotspot_df, on='hotspot_id', how='left')
    # oncodb合并的致癌性数据库，用gene+AAchange匹配CDS突变，用gene+c.匹配剪切位点突变
    logging.info('公开数据库注释信息整理 - 匹配合并的致癌性数据库')
    cbio_df = pd.read_table(args.oncodb)
    cbio_df = cbio_df.rename({'Gene': 'cbioGene', 'HGVSg': 'cbioHGVSg'}, axis=1)
    cds_df = cbio_df[~cbio_df['Mutation Type'].isin(['Splice_Site', 'Splice_Region'])].copy()
    splice_df = cbio_df[cbio_df['Mutation Type'].isin(['Splice_Site', 'Splice_Region'])].copy()
    cds_df['hotspot_id'] = cds_df['cbioGene'].str.cat(cds_df['Protein Change'], sep=':p.')  # 合并方法同Hotspot数据库
    splice_df['splice_id'] = splice_df['cbioGene'].str.cat(splice_df['cbioHGVSg'], sep=':')
    merge_df['splice_id'] = merge_df['SYMBOL'].str.cat(merge_df['HGVSg'], sep=':')
    merge_df = merge_df.merge(cds_df, on='hotspot_id', how='left')
    merge_df = merge_df.merge(splice_df, on='splice_id', how='left')
    merge_df['OncoKB'] = merge_df['OncoKB_x'].fillna(merge_df['OncoKB_y'])  # .fillna('-')
    merge_df = merge_df.drop(
        ['cbioHGVSg_x', 'Mut in Sample_x', 'cbioHGVSg_y', 'Mut in Sample_y'], axis=1
    ).drop_duplicates(keep='first')      # oncodb部分热点突变c.不同p.相同，合并后引入重复行，及时去重
    # 原癌抑癌基因数据库，用gene name匹配
    logging.info('公开数据库注释信息整理 - 匹配原癌抑癌基因数据库')
    osgene_df = pd.read_table(args.osgene)
    osgene_df['SYMBOL'] = osgene_df['基因名']
    merge_df = merge_df.merge(osgene_df, on='SYMBOL', how='left')
    # Interpro结构域去空值
    logging.info('公开数据库注释信息整理 - Interpro蛋白数据库格式整理')
    merge_df['Interpro_domain'] = merge_df['Interpro_domain'].apply(trim_interpro)
    # 获取人群频率最大值 todo:wbbc保留4位有效数字
    logging.info('公开数据库注释信息整理 - 获取最大人群频率')
    merge_df['MAX_AF'] = merge_df.parallel_apply(get_max_af, axis=1)

    # -----------有害性预测-----------
    logging.info('变异有害性预测信息整理')
    merge_df = merge_df.parallel_apply(get_mane_pred, axis=1)  # 获取MANE转录本的dbNSFP数据库注释结果
    merge_df = merge_df.copy()          # 上一步引发PerformanceWarning: DataFrame is highly fragmented警告，加.copy()解决
    merge_df['CADD预测'] = merge_df['CADD_phred_hg19'].apply(get_soft_pred, args=('cadd',))
    merge_df['FATHMM预测'] = merge_df['fathmm-XF_coding_score'].apply(get_soft_pred, args=('fathmm-xf',))
    merge_df['SIFT预测'] = merge_df['SIFT_score'].apply(get_soft2_pred, args=('sift',))
    merge_df['Polyphen2_HDIV预测'] = merge_df['Polyphen2_HDIV_score'].apply(get_soft_pred, args=('polyphen',))
    merge_df['VEST4预测'] = merge_df['VEST4_score'].apply(get_soft_pred, args=('vest4',))
    merge_df['REVEL预测'] = merge_df['REVEL_score'].apply(get_soft_pred, args=('revel',))
    merge_df['dbscSNV预测'] = merge_df.apply(get_dbscSNV_pred, axis=1)
    merge_df['MutationTaster_pred'] = merge_df['MutationTaster_pred'].str.split(',', expand=True)[0]

    # -----------生信备注列-----------
    logging.info('生信备注信息整理')
    merge_df['软件来源'] = merge_df['info'].str.extract('SOFT=(.*)')
    merge_df['超级白名单'] = merge_df['info'].apply(lambda x: '是' if 'SuperHotSpot' in x else '否')
    merge_df['重复区域'] = merge_df['info'].apply(lambda x: '是' if ';STR;' in x else '否')
    merge_df['VCF位置'] = merge_df.apply(lambda x: f'{x["chr"]}_{x["pos"]}_{x["ref"]}_{x["alt"]}', axis=1)

    # -----------肿瘤突变致病性评级-----------
    logging.info('进行肿瘤突变致病性评级')
    risk_gene_list = pd.read_excel(args.riskgene, sheet_name=0)['基因'].to_list()
    merge_df = merge_df.fillna('-')
    merge_df = merge_df.parallel_apply(pathogenicity_rating, axis=1, args=(risk_gene_list,))

    # =====================结果输出=====================
    # -----------英文列重命名并构建子Sheet-----------
    logging.info('英文列重命名并构建子Sheet')
    merge_df = merge_df.rename(columns=json_content['rename_dict'])
    merge_df = merge_df[json_content['final_col_list']]
    # key_df = merge_df.loc[merge_df['突变可靠性'] != '不可靠', json_content['key_col_list']]
    rating_df = pd.DataFrame(json_content['rating_dict'])
    # -----------写入txt和excel文件-----------
    logging.info('最终结果写入txt和excel文件')
    merge_df.to_csv(fr'{args.out}.txt', sep='\t', index=False)

    writer = pd.ExcelWriter(fr'{args.out}.xlsx')
    # key_df.to_excel(writer, sheet_name='主要信息', index=False)   # 2025年3月24日：这个sheet非必须，医学确认后去除
    merge_df.to_excel(writer, sheet_name='突变总表', index=False)
    rating_df.to_excel(writer, sheet_name='突变评级说明', index=False)
    writer.close()

    end_time = time.time()
    minutes, seconds = divmod(end_time - start_time, 60)
    logging.info(f'运行时间：{int(minutes)}分{int(seconds)}秒。')

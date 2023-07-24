#!/mnt/share01/tools/bin/python
# -*- coding: UTF-8 -*-

__author__ = "Fei Xue"
__email__ = "xuefei@zhuanhuayixue.org"
__version__ = "1a"
__date__ = "6/8/21 13:56 "
__copyright__ = "Copyright (C) 2021 All rights reserved"
__description__ = """
泛癌种：兼容panel双样本，单样本，RNA——panel， PDL1
"""
__update__ = "20210908:增加肠癌KRAS/NRAS/BRAF三基因同时阴性的处理方案" \
             "20210909:加入化疗药评分内容" \
             "20210915:兼容WES双样本过滤"


import os
import re
import json
import argparse
import pandas as pd
import numpy as np

DB_PATH = '/mnt/share01/tools/analysis_module/reports/annotation/solid_huixi/'
project_wes_list = ['NBCW', 'NCW1', 'NCW2', 'NCW3', 'NCW4', 'NCWs1', 'NCWs2', 'NCWs3', 'NCWs4', 'NPCW',
                    'NWCNS1', 'NWCNS2', 'NWCNS3', 'Ncec', 'Ncecm1', 'Ncecm2', 'Ncecm3', 'Ncet']
panel39_germline_genes = '/mnt/share02/yanyl/script/report_script/solid/pancancer/script/filter/panel39_germline_genes.xlsx'


def my_parser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--type1', help="", required=True)
    parser.add_argument('--type2', help="", default='type2')
    parser.add_argument('--panel_mutation_txt_civic_protein1', required=True)
    parser.add_argument('--panel_mutation_txt_civic_protein2', default=None)
    parser.add_argument('--fusion_file1', help="", required=True)
    parser.add_argument('--fusion_file2', help="", default=None)
    parser.add_argument('--rna_fusion_file', help="", default=None)
    parser.add_argument('--cnv_file1', help="", required=True)
    parser.add_argument('--cnv_file2', help="", default=None)
    parser.add_argument('--tumor_id1', help="肿瘤样本文库号1", required=True)
    parser.add_argument('--tumor_id2', help="肿瘤样本文库号2", default='tumor_id2')
    parser.add_argument('--tumor_wes_id1', help="wes双样本产品, wes组织样本文库号;不需要则不给出", default=None)
    parser.add_argument('--tumor_wes_id2', help="wes双样本产品，wes的ct样本文库号;不需要则不给出", default=None)
    parser.add_argument('--normal_id', help="对照样本文库号", required=True)
    parser.add_argument('--normal_bam_readcount', required=True)
    parser.add_argument('--msi_summary1', help="肿瘤样本文库号", required=True)
    parser.add_argument('--msi_summary2', help="肿瘤样本文库号", default=None)
    parser.add_argument('--PD_L1_json', help="PD_L1 json文件", required=True)
    parser.add_argument('--hla_file', help="hla分析结果文件", required=True)

    parser.add_argument('--cnv_1p19q_1', help="脑胶质瘤 NCCN 才用")
    parser.add_argument('--cnv_1p19q_2', default=None)
    parser.add_argument('--product', help="脑胶质瘤 NCCN 才用")

    parser.add_argument('--LimsInfo', help="LimsInfo.json", required=True)
    parser.add_argument('--project_id', help="project_id NCP5/NCP6", required=True)
    parser.add_argument('--outdir', help="", required=True)

    parser.add_argument('--castration', help="castration result", default=None)  # 前列腺癌去势治疗药物结果

    # 阈值
    parser.add_argument('--sv_read', help="重排支持reads数", type=int, default=6)
    parser.add_argument('--sv_rate', help="重排丰度", type=float, default=0.0005)  # 0.0005 四舍五入之后 >=0.1%

    # 知识库
    parser.add_argument('--lung_gene_list', default=f'{DB_PATH}/肺癌基因list.xlsx')
    parser.add_argument('--solid_tumor_targeting_site_drug',
                        default=f'{DB_PATH}/慧系版报告-solid_tumor_targeting_site_drug.xlsx')
    parser.add_argument('--solid_tumor_chemotherapy_annotation',
                        default=f'{DB_PATH}/慧系版报告-blood_tumor_chemotherapy_annotation.xlsx')
    parser.add_argument('--solid_tumor_ClinicalTrials', default=f'{DB_PATH}/慧系版报告-solid_tumor_ClinicalTrials.xlsx')
    parser.add_argument('--solid_tumor_gene_annotation', default=f'{DB_PATH}/慧系版报告-solid_tumor_gene_annotation.xlsx')
    parser.add_argument('--solid_tumor_chemotherapy_drug_list',
                        default=f'{DB_PATH}/慧系版报告-blood_tumor_chemotherapy_drug_list.xlsx')
    parser.add_argument('--solid_tumor_cancer_relation', default=f'{DB_PATH}/慧系版报告-疾病对照-solid_tumor.xlsx')
    parser.add_argument('--solid_tumor_NCCN_gene', default=f'{DB_PATH}/慧系版报告-solid_tumor_NCCN指南推荐基因汇总.xlsx')
    parser.add_argument('--solid_tumor_immunotherapy_gene', default=f'{DB_PATH}/慧系版报告-solid_tumor_免疫治疗疗效相关基因.xlsx')
    parser.add_argument('--solid_tumor_linchuangshiyan', default=f'{DB_PATH}/慧系版报告-solid_tumor_临床试验.xlsx')
    parser.add_argument('--solid_tumor_pathway', default=f'{DB_PATH}/brain_信号通路数据库.xlsx')
    parser.add_argument('--solid_tumor_pathway_393', default=f'{DB_PATH}/慧系版报告-信号通路数据库-393基因-20210917.xlsx')

    parser.add_argument('--tumor_genetic_predisposed_gene_list',
                        default='/mnt/share01/tools/analysis_module/reports/annotation/blood/zhixuan/germline_anno/tumor_genetic_predisposed_gene_list.txt')
    parser.add_argument('--solid_gene_list', default=f'{DB_PATH}/393-39各模块基因list.xlsx')
    parser.add_argument('--cancer', help="疾病大类，eg. 脑胶质/肺/肠/前列腺")
    args = parser.parse_args()
    return args


class Tools:

    @staticmethod
    def tmb_process(mutation_frame, sample_id, outdir, project_id):
        tmb_mutation_frame = mutation_frame[
            (mutation_frame["突变可靠性"] == "相对可靠") &
            (mutation_frame["结论"].isin(["体细胞较可靠突变,可能无害", "体细胞较可靠突变,可能有害"])) &
            (~mutation_frame["突变类型"].isin(["非编码区突变", "剪切位点附近的突变", "剪切位点突变", "同义突变"])) &
            (mutation_frame["突变率"] >= 0.05)]
        # tmb_mutation_frame.to_csv('tmb_mutation_frame.txt', sep='\t', index=False)
        tmb_num = tmb_mutation_frame.shape[0]

        if project_id in project_wes_list:
            tmb_value = tmb_num / 33
        else:
            tmb_value = tmb_num / 1.9
        open(f"{outdir}/TMB_stat_{sample_id}.xls", 'wt').write("sample_id\tTMB_somatic_count\tTMB\n%s\t%d\t%.2f" % (sample_id, tmb_num, tmb_value))
        return tmb_value

    @staticmethod
    def extrac_exon(x):
        need_list = []
        exon_list = re.findall('\(([a-z]*[A-Z]*[0-9]+)\)', x)
        for each in exon_list:
            if 'exon' in each:
                need_list.append(each)
        if len(need_list) == 0:
            return ','.join(set(exon_list))
        else:
            return ','.join(set(need_list))

    @staticmethod
    def rm_exon(x):
        x['核酸改变'] = re.sub('\([a-z]*[A-Z]*[0-9]+\)', '', x['核酸改变'])
        return x

    @staticmethod
    def source_of_evidence(x):
        """
        FDA NMPA NCCN 会议/共识 临床试验 个例报道'
        """
        x = x.fillna(value='')
        return_list = []
        if x['疾病'] in x['FDA']:
            return_list.append('FDA')
        if x['疾病'] in x['NMPA']:
            return_list.append('NMPA')
        if x['疾病'] in x['NCCN']:
            return_list.append('NCCN')
        if x['疾病'] in x['会议/共识']:
            return_list.append('会议共识')
        if x['疾病'] in x['临床试验']:
            return_list.append('临床试验')
        if x['疾病'] in x['个例报道']:
            return_list.append('个例报道')
        if len(return_list) == 0:
            return x['疾病'] + ':' + '-'
        else:
            return x['疾病'] + ':' + '/'.join(return_list)

    @staticmethod
    def cluster_columns(df, sep=','):
        return sep.join([str(x) for x in set(df.values)])

    @staticmethod
    def mut_status(x, gender, normal_id):
        depth = "突变深度(对照,{})".format(normal_id)
        freq = "突变率(对照,{})".format(normal_id)
        if x["染色体名称"] == "chrX":
            if gender == "男":
                if x[freq] > 0.75 and x[depth] >= 5:
                    return "Hem"
                else:
                    return "."
            else:
                if x[depth] >= 5:
                    if 0.1 < x[freq] < 0.85:
                        return "Het"
                    elif x[freq] >= 0.85:
                        return "Hom"
                    else:
                        return "."
                else:
                    return "."
        else:
            if x[depth] >= 5:
                if 0.1 < x[freq] < 0.85:
                    return "Het"
                elif x[freq] >= 0.85:
                    return "Hom"
                else:
                    return "."
            else:
                return "."

    @staticmethod
    def get_pathogenicity(x):
        set_info = list(set(str(x['Clinvar']).split(';')))
        if 'Pathogenic' in set_info:
            return '致病'
        elif 'Likely pathogenic' in set_info:
            return '可能致病'
        elif ['NA'] == set_info:
            return '-'
        else:
            return '无注释结果'

    @staticmethod
    def get_benign(x):
        set_info = list(set(str(x['Clinvar']).split(';')))
        if 'Benign' in set_info:
            flag = '良性'
        elif 'Likely benign' in set_info:
            flag = '可能良性'
        elif ['NA'] == set_info:
            flag = '-'
        else:
            flag = '无注释结果'
        return flag

    @staticmethod
    def convet_cosmic(x):
        if x['COSMIC数据库中有此突变的组织'] in ['', ' ', 'NA']:
            return "-"
        else:
            return "有"

    # 20210914更新兼容WES样本分析
    @staticmethod
    def get_conclusion(x, product, panel_id=None, wes_id=None):
        if product == 'wes':
            if pd.isnull(x[f"结论(待定,{panel_id})"]) or x[f"结论(待定,{panel_id})"] == 'NA' or (x[f"结论(待定,{panel_id})"] == '突变来源不确定'):
                return x["结论(待定,%s)" % wes_id]
            else:
                return x["结论(待定,%s)" % panel_id]
        else:
            return x["结论"]

    @staticmethod
    def get_vaf(x, product, panel_id=None, wes_id=None):
        if product == 'wes':
            if pd.isnull(x["突变率(待定,%s)" % wes_id]) or x["突变率(待定,%s)" % wes_id] == 'NA':
                return x["突变率(待定,%s)" % panel_id]
            else:
                if pd.isnull(x["突变率(待定,%s)" % panel_id]) or x["突变率(待定,%s)" % panel_id] == 'NA':
                    return x["突变率(待定,%s)" % wes_id]
                else:
                    if x["突变深度(待定,%s)" % panel_id] > x["突变深度(待定,%s)" % wes_id]:
                        return x["突变率(待定,%s)" % panel_id]
                    else:
                        return x["突变率(待定,%s)" % wes_id]
        else:
            return x["突变率(待定,%s)" % panel_id]

    @staticmethod
    def get_depth(x, product, panel_id=None, wes_id=None):
        if product == 'wes':
            if pd.isnull(x["突变深度(待定,%s)" % wes_id]) or x["突变深度(待定,%s)" % wes_id] == 'NA':
                return x["突变深度(待定,%s)" % panel_id]
            else:
                if pd.isnull(x["突变深度(待定,%s)" % panel_id]) or x["突变深度(待定,%s)" % panel_id] == 'NA':
                    return x["突变深度(待定,%s)" % wes_id]
                else:
                    if x["突变深度(待定,%s)" % panel_id] > x["突变深度(待定,%s)" % wes_id]:
                        return x["突变深度(待定,%s)" % panel_id]
                    else:
                        return x["突变深度(待定,%s)" % wes_id]
        else:
            return x["突变深度(待定,%s)" % panel_id]


class CancerFilter:

    def __init__(self, argv):
        self.argv = argv
        self.type1 = self.argv.type1
        self.type2 = self.argv.type2
        self.mutation_file1 = self.argv.panel_mutation_txt_civic_protein1
        self.mutation_file2 = self.argv.panel_mutation_txt_civic_protein2
        self.fusion_file1 = self.argv.fusion_file1
        self.fusion_file2 = self.argv.fusion_file2
        self.cnv_file1 = self.argv.cnv_file1
        self.cnv_file2 = self.argv.cnv_file2
        self.cnv_1p19q_1 = self.argv.cnv_1p19q_1
        self.cnv_1p19q_2 = self.argv.cnv_1p19q_2
        self.normal_bam_readcount = self.argv.normal_bam_readcount
        self.castration = self.argv.castration
        self.tumor_id1 = self.argv.tumor_id1
        self.tumor_id2 = self.argv.tumor_id2
        self.tumor_wes_id1 = self.argv.tumor_wes_id1
        self.tumor_wes_id2 = self.argv.tumor_wes_id2
        self.normal_id = self.argv.normal_id
        self.msi_summary1 = self.argv.msi_summary1
        self.msi_summary2 = self.argv.msi_summary2
        self.PD_L1_json = json.load(open(self.argv.PD_L1_json))
        self.gender = json.load(open(self.argv.LimsInfo))['data']['gender']
        self.hla_file = self.argv.hla_file
        self.product = self.argv.product
        self.project_id = self.argv.project_id
        self.outdir = os.path.realpath(self.argv.outdir) + '/'
        self.sv_read = self.argv.sv_read
        self.sv_rate = self.argv.sv_rate

        # rna 融合处理
        if os.path.isfile(self.argv.rna_fusion_file) and os.path.getsize(self.argv.rna_fusion_file) != 0:
            pd_rna_fusion_file = pd.read_csv(self.argv.rna_fusion_file, sep='\t')
        else:
            pd_rna_fusion_file = pd.DataFrame()

        if pd_rna_fusion_file.empty:
            self.rna_fusion_file = ''
        else:
            self.rna_fusion_file = self.argv.rna_fusion_file

        self.result_file_list = []

        if not self.mutation_file2:
            self.single_flag = True
        else:
            self.single_flag = False

        self.cancer = self.argv.cancer
        self.cancer_flag = self.argv.cancer

        self.solid_tumor_targeting_site_drug = self.argv.solid_tumor_targeting_site_drug
        self.lung_gene_list = self.argv.lung_gene_list
        self.solid_tumor_chemotherapy_annotation = self.argv.solid_tumor_chemotherapy_annotation
        self.solid_tumor_ClinicalTrials = self.argv.solid_tumor_ClinicalTrials
        self.solid_tumor_gene_annotation = self.argv.solid_tumor_gene_annotation
        self.solid_tumor_chemotherapy_drug_list = self.argv.solid_tumor_chemotherapy_drug_list
        self.solid_tumor_cancer_relation = self.argv.solid_tumor_cancer_relation
        self.solid_tumor_NCCN_gene = self.argv.solid_tumor_NCCN_gene
        self.solid_tumor_immunotherapy_gene = self.argv.solid_tumor_immunotherapy_gene
        self.solid_tumor_linchuangshiyan = self.argv.solid_tumor_linchuangshiyan
        self.solid_tumor_pathway = self.argv.solid_tumor_pathway
        self.solid_tumor_pathway_393 = self.argv.solid_tumor_pathway_393
        self.tumor_genetic_predisposed_gene_list = self.argv.tumor_genetic_predisposed_gene_list
        self.solid_gene_list = self.argv.solid_gene_list
        self.pd_solid_gene_list = pd.read_excel(self.solid_gene_list, header=None, sheet_name=None)

        # ['393基因list', '393-28融合基因list', '39基因list', '39-6融合基因list', '39-12拷贝数变异基因list',
        # '39-10胚系基因list', '前列腺101', '前列腺13化疗药基因', '尿路上皮98', '尿路上皮17化疗药基因', '肾81'
        # '肾32遗传性基因']

        if self.argv.cancer == '结直肠':
            self.mut_filter_gene = [x for x in sum(self.pd_solid_gene_list['39基因list'].fillna(value='').values.tolist(), []) if x != '']
            self.cnv_filter_gene = [x for x in sum(self.pd_solid_gene_list['39-12拷贝数变异基因list'].fillna(value='').values.tolist(), []) if x != '']
            self.fus_filter_gene = [x for x in sum(self.pd_solid_gene_list['39-6融合基因list'].fillna(value='').values.tolist(), []) if x != '']
            self.ger_filter_gene = [x for x in sum(self.pd_solid_gene_list['39-10胚系基因list'].fillna(value='').values.tolist(), []) if x != '']
    # else:
    # 	pass
    # [x for x in sum(self.pd_solid_gene_list['393基因list'].fillna(value='').values.tolist(), []) if x != '']
    # [x for x in sum(self.pd_solid_gene_list['393-28融合基因list'].fillna(value='').values.tolist(), []) if x != '']

    def filter_little_panel(self, pd_final):
        ''' add by rhd: 增加小panel基因列表过滤 '''
        little_panel = ['NCPP1', 'NCPP2', 'NCPP3', 'NCPK1', 'NCPK2', 'NCPK3', 'NCPU1', 'NCPU2', 'NCPU3']# NCPLt1 NCPLt2 NCPLt3
        sheet_dict = {'前列腺' : '前列腺101', '尿路上皮' : '尿路上皮98', '肾' : '肾81'}
        if self.project_id in little_panel:
            pd_panel_gene = pd.read_excel(self.solid_gene_list, sheet_name=sheet_dict[self.cancer], header=None)# 库表 癌症为小类
            filter_gene_list = [x[0] for x in pd_panel_gene.values.tolist()]
            pd_final = pd_final[pd_final['基因'].isin(filter_gene_list)]
        return pd_final

    def filter_mut(self, if_do_targeting_site_drug=True):

        def filter1(x):
            """ 错义突变 """
            _gene_name = x['基因']
            _type_mut = str(x['氨基酸改变']).replace('p.', '')
            if re.search('^[A-Z][0-9]+[A-Z]$', _type_mut):
                final_type_mut = _type_mut
            elif re.search('^[A-Z][0-9]+X,[0-9]+$', _type_mut):
                final_type_mut = _type_mut.split(',')[0]
            else:
                final_type_mut = _type_mut
            return final_type_mut

        def filter2(x):
            """ 去掉我们结果中fs前面的字母 """
            _gene_name = x['基因']
            _type_mut = str(x['氨基酸改变']).replace('p.', '')
            if _gene_name in type_fs + type_fs_balabala and 'fs' in _type_mut:
                if _gene_name in type_fs and 'fs*' in _type_mut:
                    final_type_mut = _type_mut.split('fs*')[0] + 'fs'
                else:
                    final_type_mut = _type_mut
            else:
                final_type_mut = _type_mut

            return final_type_mut

        def filter3(x):
            """ 特殊检测 需要判断位点范围 且流程突变结果中有相应的关键词"""
            my_gene = x['基因']
            my_type_mut = str(x['氨基酸改变']).replace('p.', '')
            my_chr = str(x['染色体名称'])
            my_start = int(x['起始位置'])
            my_end = int(x['终止位置'])
            final_type_mut = [my_type_mut]
            for db_gene, position in type_special.keys():
                # print(db_gene, position)
                dict_key_word = {}
                value_list = str(type_special[(db_gene, position)]).split(';')
                for each in value_list:
                    # if 'del' in each and 'ins' in each:
                    # 	dict_key_word[each] = 'del'
                    if 'del' in each:
                        dict_key_word[each] = 'del'
                    elif 'ins' in each:
                        dict_key_word[each] = 'ins'
                    else:
                        dict_key_word[each] = each
                db_chr = str(position.split(':')[0])
                db_start = int(position.split(':')[1].split('-')[0])
                db_end = int(position.split(':')[1].split('-')[1])
                if db_gene != my_gene:
                    continue
                # print(position, my_chr, db_chr, db_start, db_end, my_start, db_end)
                if (my_chr == db_chr) and (db_start <= my_start <= db_end) and (db_start <= my_end <= db_end):
                    for each_key_word in dict_key_word.keys():
                        """
                        each_key_word:  exon19del
                        each_key_value: del
                        """
                        each_key_value = dict_key_word[each_key_word]
                        if 'exon14skipping' == each_key_word:
                            final_type_mut.append('exon14skipping')
                        elif 'Promoter' == each_key_word and '启动子' in str(x['突变类型']):
                            final_type_mut.append('Promoter')
                        elif ('splice' in each_key_word) and (str(x['突变类型']) == '剪切位点突变'):
                            final_type_mut.append(each_key_word)
                        # elif 'Oncogenic Mutations' == each_key_word and (
                        # 		('Pathogenic' in str(x['Clinvar']).split(';') or 'Likely pathogenic' in str(x['Clinvar']).split(';')) or
                        # 		('CLASS=DM' in x['Hgmd'].split(';') or 'CLASS=DM?' in x['Hgmd'].split(';')) or
                        # 		(x['SIFT危害性预测'].split('(')[0] == 'D' and
                        # 		 x['Polyphen2_HDIV危害性预测'].split('(')[0] == 'D' and
                        # 		 x['Polyphen2_HVAR危害性预测'].split('(')[0] == 'D' and
                        # 		 x['PROVEAN危害性预测'].split('(')[0] == 'D' and
                        # 		 x['MutationTaster危害性预测'].split('(')[0] == 'D' and
                        # 	 	 x['M-CAP危害性预测'].split('(')[0] == 'D' and
                        # 	 	 x['M-REVEL危害性预测'].split('(')[0] == 'D') or
                        # 		('其他能预测为致病/可能致病的变异，即未被已知数据库如clinvar收录注释为良性/可能良性/意义未明的缺失、移码或者无义变异')
                        # ):
                        # 	pass
                        elif 'del' in my_type_mut and 'ins' in my_type_mut:
                            final_type_mut.append('del')
                        elif each_key_value in my_type_mut:
                            final_type_mut.append(each_key_word)
                    # else:
                    # 	pass
            return ';'.join(final_type_mut)

        return_tuple = []
        for_list = [[self.mutation_file1, self.tumor_id1, self.tumor_wes_id1], [self.mutation_file2, self.tumor_id2, self.tumor_wes_id2]]
        for (mutation_file, tumor_id, tumor_wes_id) in for_list:
            if not mutation_file:
                return_tuple.append(None)
                continue
            pd_panel_mutation = pd.read_csv(mutation_file, sep='\t', low_memory=False)
            pd_panel_mutation['外显子'] = pd_panel_mutation['核酸改变'].apply(lambda x: Tools.extrac_exon(x))
            pd_panel_mutation['结论'] = pd_panel_mutation.apply(lambda x: Tools.get_conclusion(x, self.product, panel_id=tumor_id, wes_id=tumor_wes_id), axis=1)
            # pd_panel_mutation.to_csv('突变改结论后.txt', sep='\t', index=False)
            pd_panel_mutation['突变率'] = pd_panel_mutation.apply(lambda x: Tools.get_vaf(x, self.product, panel_id=tumor_id, wes_id=tumor_wes_id), axis=1)
            pd_panel_mutation = pd_panel_mutation.apply(lambda x: Tools.rm_exon(x), axis=1)
            # pd_panel_mutation = pd_panel_mutation[
            # 	(pd_panel_mutation['结论'].isin(['体细胞较可靠突变,可能有害', '体细胞突变,可能有害'])) &
            # 	(pd_panel_mutation['突变可靠性'] == '相对可靠') &
            # 	(~pd_panel_mutation['突变类型'].isin(['同义突变', '非编码区突变', '剪切位点附近突变']))
            # 	]
            pd_panel_mutation_copy = pd_panel_mutation.copy()
            pd_panel_mutation = pd_panel_mutation[pd_panel_mutation['突变可靠性'] == '相对可靠']
            need_col_index = ['突变编号', '基因', '突变率', 'mRNA变体名称(NM)', '染色体名称', '起始位置', '终止位置', '核酸改变', '氨基酸改变', '突变类型',
                              '结论', 'COSMIC数据库中有此突变的组织', 'SIFT危害性预测', 'Polyphen2_HDIV危害性预测', 'REVEL危害性预测',
                              'MutationTaster危害性预测', '外显子', 'Clinvar', 'Hgmd']
            # 提取相应需要的列
            pd_panel_mutation = pd_panel_mutation[need_col_index]
            pd_panel_mutation['突变率'] = pd_panel_mutation['突变率'].map(lambda x: format(x, '.1%'))
            pd_panel_mutation = pd_panel_mutation.fillna(value='')
            pd_panel_mutation.to_csv('突变基础过滤后.txt', sep='\t', index=False)

            if not if_do_targeting_site_drug:
                pd_panel_mutation = self.filter_little_panel(pd_panel_mutation)
                return_tuple.append(pd_panel_mutation)
            else:
                """ 库文件 读取和准备"""
                solid_tumor_targeting_site_drug = pd.read_excel(self.solid_tumor_targeting_site_drug)
                solid_tumor_targeting_site_drug['证据来源'] = solid_tumor_targeting_site_drug.apply(
                    lambda x: Tools.source_of_evidence(x), axis=1)
                """ 下面循环是 总结库中不同修改"氨基酸改变"的基因，以便于后面判断哪些基因是需要进行哪种 修改"氨基酸改变"的方式"""
                type_fs = []
                type_fs_balabala = []
                type_special = {}
                for index, info in solid_tumor_targeting_site_drug.iterrows():
                    gene_name = info['基因']
                    type_mut = info['变异']

                    if info['变异类型'] == '缺失' and re.search('^[A-Z][0-9]+[A-Z]\*fs$', type_mut):
                        type_fs.append(gene_name)
                    elif info['变异类型'] == '缺失' and re.search('^[A-Z][0-9]+[A-Z]fs\*[0-9]+$', type_mut):
                        type_fs_balabala.append(gene_name)

                # if info['变异类型'] == '特殊检测':  # and info['疾病关键字'] == self.cancer_flag
                # 	position = info['位置']
                # 	if gene_name not in type_special:
                # 		type_special[gene_name] = {}
                # 	if position not in type_special[gene_name]:
                # 		type_special[gene_name][position] = type_mut

                pd_db_for_filter3 = solid_tumor_targeting_site_drug[solid_tumor_targeting_site_drug['变异类型'] == '特殊检测']
                pd_db_for_filter3 = pd_db_for_filter3.groupby(['基因', '位置'])['变异'].apply(lambda x: Tools.cluster_columns(x, sep=';'))
                type_special = pd_db_for_filter3.to_dict()
                """ 先进行基础的过滤，又因为突变结果文件比较大，所以再根据库文件里面的基因对突变文件进行过滤一下"""
                result = pd_panel_mutation.copy()
                result = result[result['基因'].isin(solid_tumor_targeting_site_drug['基因'])]

                result1 = result.copy()
                result2 = result.copy()
                result3 = result.copy()

                """ 使用三种 "氨基酸改变"修改方式进行修改结果为"氨基酸改变2"，这样方便和库文件直接进行 pd.merge """
                result1['氨基酸改变2'] = ''
                result1['氨基酸改变2'] = result1.apply(lambda x: filter1(x), axis=1)
                # result1.to_csv('result1.txt', sep='\t', index=False)
                result2['氨基酸改变2'] = ''
                result2['氨基酸改变2'] = result2.apply(lambda x: filter2(x), axis=1)
                # result2.to_csv('result2.txt', sep='\t', index=False)
                result3['氨基酸改变2'] = ''
                result3['氨基酸改变2'] = result3.apply(lambda x: filter3(x), axis=1)
                # result3.to_csv('result3-1.txt', sep='\t', index=False)
                result3 = result3.drop('氨基酸改变2', axis=1).join(
                    result3['氨基酸改变2'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename(
                        '氨基酸改变2'))
                # result3.to_csv('result3-2.txt', sep='\t', index=False)

                """ 三种文件进行合并，去除重复的行。注意这里是允许有重复的基因但是有不同的变异分类"""
                final_result = pd.concat([result1, result2, result3], join='inner')
                final_result = final_result.drop_duplicates()
                """ 准备好三种方法 得到合并的总突变文件之后，和库文件进行合并比较"""
                pd_mutation_final = pd.merge(solid_tumor_targeting_site_drug, final_result, left_on=['基因', '变异'],
                                             right_on=['基因', '氨基酸改变2'], how='inner')
                pd_mutation_final = pd_mutation_final.fillna(value='')
                pd_mutation_final = self.filter_little_panel(pd_mutation_final)
                pd_mutation_final.to_csv('突变过库后.txt', sep='\t', index=False)
                #### 如果xxx不足16个
                if pd_mutation_final.shape[0] < 16:
                    pd_panel_mutation_copy = pd_panel_mutation_copy[
                        (pd_panel_mutation_copy['突变可靠性'] == '相对可靠') &
                        (pd_panel_mutation_copy['结论'].isin(['体细胞较可靠突变,可能有害', '体细胞突变,可能有害'])) &
                        (~pd_panel_mutation_copy['突变类型'].isin(['同义突变', '非编码区突变', '剪切位点附近突变']))]
                    need_col_index = ['突变编号', '基因', '突变率', 'mRNA变体名称(NM)', '染色体名称', '起始位置', '终止位置', '核酸改变', '氨基酸改变', '突变类型',
                              '结论', 'COSMIC数据库中有此突变的组织', 'SIFT危害性预测', 'Polyphen2_HDIV危害性预测', 'REVEL危害性预测',
                              'MutationTaster危害性预测', '外显子', 'Clinvar', 'Hgmd']
                    pd_panel_mutation_copy = pd_panel_mutation_copy[need_col_index]
                    pd_panel_mutation['突变率'] = pd_panel_mutation['突变率'].map(lambda x: format(x, '.1%'))
                    pd_panel_mutation = pd_panel_mutation.fillna(value='')
                    pd_panel_mutation.to_csv('突变基础过滤后.txt', sep='\t', index=False)
                    #### 写到这里了 ######
            # exit(0)
                return_tuple.append(pd_mutation_final)
        return return_tuple

    def filter_cnv(self, if_do_targeting_site_drug=True):
        return_list = []
        for_list = [self.cnv_file1, self.cnv_file2]
        for cnv_file in for_list:
            if not cnv_file:
                return_list.append(None)
                continue
            pd_cnv_file = pd.read_excel(cnv_file)
            "chromosome	start	end	gene	log2	cn	depth	p_ttest	probes	weight"
            pd_cnv_file = pd_cnv_file[(pd_cnv_file['cn'] >= 4) & ~pd_cnv_file['gene'].isna()]  # 只要拷贝数扩增的
            if not pd_cnv_file.empty:
                pd_cnv_out = pd_cnv_file.drop('gene', axis=1).join(pd_cnv_file['gene'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('gene'))
                pd_cnv_out = pd_cnv_out[['gene', 'cn']]
                pd_cnv_out.columns = ['基因', '拷贝数']
                pd_cnv_out['变异（缺失/扩增）'] = '拷贝数扩增'
                pd_cnv_out = pd_cnv_out.fillna(value='')
            else:
                pd_cnv_out = pd.DataFrame(columns=['基因', '变异（缺失/扩增）', '拷贝数'])

            if if_do_targeting_site_drug:
                """
                solid_tumor_targeting_site_drug 文件读取与准备
                """
                db_solid_tumor_targeting_site_drug = pd.read_excel(self.solid_tumor_targeting_site_drug)
                db_solid_tumor_targeting_site_drug = db_solid_tumor_targeting_site_drug[
                    db_solid_tumor_targeting_site_drug['变异类型'] == '拷贝数变异']
                db_solid_tumor_targeting_site_drug['证据来源'] = db_solid_tumor_targeting_site_drug.apply(
                    lambda x: Tools.source_of_evidence(x), axis=1)

                """ 和 数据库 进行比较 """
                pd_cnv_out = pd.merge(db_solid_tumor_targeting_site_drug, pd_cnv_out, on='基因', how='inner')
                pd_cnv_out = pd_cnv_out.fillna(value='')
                pd_cnv_out = self.filter_little_panel(pd_cnv_out)
                return_list.append(pd_cnv_out)
            else:
                pd_cnv_out = self.filter_little_panel(pd_cnv_out)
                return_list.append(pd_cnv_out)

        return return_list

    def filter_fus(self, if_do_targeting_site_drug=True):
        return_list = []
        for_list = [[self.fusion_file1, self.type1], [self.fusion_file2, self.type2]]
        for fusion_file, sample_type in for_list:
            if not fusion_file:
                return_list.append(None)
                continue
            pd_fusion_file = pd.read_csv(fusion_file, sep='\t')
            pd_fusion_file = pd_fusion_file[['基因', 'fusion_gene_name', 'value']]

            if if_do_targeting_site_drug:
                header = ['基因', 'fusion_gene_name', 'value', '疾病关键字', '药物', '敏感性预测', '证据来源']
                """ （ solid_tumor_targeting_site_drug ） 文件读取和准备 """
                db_solid_tumor_targeting_site_drug = pd.read_excel(self.solid_tumor_targeting_site_drug)
                db_solid_tumor_targeting_site_drug = db_solid_tumor_targeting_site_drug[
                    db_solid_tumor_targeting_site_drug['变异类型'] == '融合']

                db_solid_tumor_targeting_site_drug['证据来源'] = db_solid_tumor_targeting_site_drug.apply(
                    lambda x: Tools.source_of_evidence(x), axis=1)

                pd_fusion_merge = pd.merge(db_solid_tumor_targeting_site_drug, pd_fusion_file, on='基因', how='inner')
                pd_fusion_merge = pd_fusion_merge[header]
                pd_fusion_merge = pd_fusion_merge.fillna(value='')
                pd_fusion_merge = pd_fusion_merge.sort_values('value').drop_duplicates(
                    subset=['基因', 'fusion_gene_name', '疾病关键字', '药物', '敏感性预测', '证据来源'], keep='last')
            else:
                pd_fusion_merge = pd_fusion_file[['基因', 'value', 'fusion_gene_name']]
                pd_fusion_merge = pd_fusion_merge.fillna(value='')
                pd_fusion_merge = pd_fusion_merge.sort_values('value').drop_duplicates(
                    subset=['基因', 'fusion_gene_name'], keep='last')
                pd_fusion_merge = pd_fusion_merge.drop_duplicates()

            if not pd_fusion_merge.empty:
                pd_fusion_merge['value'] = pd_fusion_merge['value'].map(lambda x: format(float(x) / 100, '.1%'))
            pd_fusion_merge = self.filter_little_panel(pd_fusion_merge)
            return_list.append(pd_fusion_merge)

        return return_list

    def filter_rna(self, if_do_targeting_site_drug=True):
        pd_rna_fus = pd.read_csv(self.rna_fusion_file, sep='\t')
        pd_rna_fus = pd_rna_fus.sort_values('FFPM').drop_duplicates(subset=['融合基因'], keep='last')
        pd_rna_fus = pd_rna_fus[['融合基因', 'FFPM']]
        pd_rna_fus.columns = ['fusion_gene_name', 'value']
        pd_rna_fus['fusion_gene_name'] = [str(x).replace('--', '/') for x in pd_rna_fus['fusion_gene_name']]
        pd_rna_fus['value'] = 'FFPM=' + pd_rna_fus['value'].map(str)
        if not pd_rna_fus.empty:
            pd_rna_fus['基因'] = pd_rna_fus['fusion_gene_name'].str.split('/', expand=True)[0]
            pd_rna_fus = pd_rna_fus[['基因', 'value', 'fusion_gene_name']]
        else:
            print('空')
            pd_rna_fus = pd.DataFrame(columns=['基因', 'value', 'fusion_gene_name'])
            self.rna_fusion_file = ''

        if if_do_targeting_site_drug:
            header = ['基因', 'fusion_gene_name', 'value', '疾病关键字', '药物', '敏感性预测', '证据来源']
            """ （ solid_tumor_targeting_site_drug ） 文件读取和准备 """
            db_solid_tumor_targeting_site_drug = pd.read_excel(self.solid_tumor_targeting_site_drug)
            db_solid_tumor_targeting_site_drug = db_solid_tumor_targeting_site_drug[
                db_solid_tumor_targeting_site_drug['变异类型'] == '融合']
            db_solid_tumor_targeting_site_drug['证据来源'] = db_solid_tumor_targeting_site_drug.apply(
                lambda x: Tools.source_of_evidence(x), axis=1)
            pd_rna_fus = pd.merge(db_solid_tumor_targeting_site_drug, pd_rna_fus, on='基因', how='inner')
            pd_rna_fus = pd_rna_fus[header]
        pd_rna_fus = self.filter_little_panel(pd_rna_fus)
        return pd_rna_fus

    def unit_1_1(self):
        """
        一、检测结果汇总 --- 靶向药物+其它重要基因相关检测结果
            靶向药物+其它重要基因
        """

        pd_mut_final1, pd_mut_final2 = self.filter_mut(if_do_targeting_site_drug=True)
        pd_cnv_final1, pd_cnv_final2 = self.filter_cnv(if_do_targeting_site_drug=True)
        pd_fus_final1, pd_fus_final2 = self.filter_fus(if_do_targeting_site_drug=True)

        pd_mut_final1.to_csv(f'{self.outdir}/all_drug_site1.txt', index=False, sep='\t')
        if not self.single_flag:
            pd_mut_final2.to_csv(f'{self.outdir}/all_drug_site2.txt', index=False, sep='\t')

        # mut
        header1 = ['基因', '外显子', '氨基酸改变', '突变率', '敏感性预测', '药物', '疾病关键字']
        header2 = ['基因', '外显子', '氨基酸改变', '突变率', '敏感性预测', '药物', '疾病关键字']
        if self.single_flag:
            pd_mut_final2 = pd.DataFrame(columns=header2)

        pd_mut_final1 = pd_mut_final1[header1].drop_duplicates()
        pd_mut_final2 = pd_mut_final2[header2].drop_duplicates()
        pd_mut_final1.columns = ['变异基因', '外显子', '氨基酸改变', self.type1, '敏感性预测', '药物', '疾病关键字']
        pd_mut_final2.columns = ['变异基因', '外显子', '氨基酸改变', self.type2, '敏感性预测', '药物', '疾病关键字']
        pd_mut = pd.merge(pd_mut_final1, pd_mut_final2, how='outer', on=['变异基因', '外显子', '氨基酸改变', '敏感性预测', '药物', '疾病关键字'])
        pd_mut = pd_mut[~(pd_mut['敏感性预测'] == '')].reset_index(drop=True)
        pd_mut['检测结果'] = [f"{pd_mut['外显子'][x]} {pd_mut['氨基酸改变'][x]}" for x in range(len(pd_mut))]
        pd_mut['类型'] = ['其它重要基因相关突变' if x in ['', ' '] else '靶向药物相关突变' for x in list(pd_mut['药物'])]
        pd_mut_target = pd_mut[pd_mut['类型'] == '靶向药物相关突变']
        pd_mut_other =  pd_mut[pd_mut['类型'] =='其它重要基因相关突变']
        pd_mut_other = pd_mut_other[pd_mut_other['疾病关键字'] == self.cancer_flag]
        pd_mut = pd.concat([pd_mut_target, pd_mut_other])
        pd_mut = pd_mut[['类型', '变异基因', '检测结果', self.type1, self.type2, '敏感性预测', '药物']]

        pd_mut = pd_mut.fillna(value='-')

        # 内分泌 add by rhd
        if '前列腺' in self.argv.cancer and self.project_id not in ['NCPM5', 'NCPM6', 'NCPM7', 'NCP393', 'NCP393c']:
            # print ('cancer type: 前列腺')
            '''输入文件为merge表+CNV（AR-V7扩增）'''
            '''AR-V7 突变结果外送，lims抓取；同时考虑AR-V7突变和拷贝数扩增'''
            pd_endo_final1, pd_endo_final2  = self.filter_mut(if_do_targeting_site_drug=True)
            pd_ARV7_final1, pd_ARV7_final2 = self.filter_cnv(if_do_targeting_site_drug=True)
            header1 = ['基因', '外显子', '氨基酸改变', '突变率', '内分泌药物敏感性', '内分泌药物']
            header2 = ['基因', '外显子', '氨基酸改变', '突变率', '内分泌药物敏感性', '内分泌药物']
            header3 = ['基因', '拷贝数', '内分泌药物敏感性', '内分泌药物']
            header4 = ['基因', '拷贝数', '内分泌药物敏感性', '内分泌药物']
            if self.single_flag:
                pd_endo_final2 = pd.DataFrame(columns=header2)
                pd_ARV7_final2 = pd.DataFrame(columns=header4)

            pd_endo_final1 = pd_endo_final1[header1].drop_duplicates()
            pd_endo_final2 = pd_endo_final2[header2].drop_duplicates()
            pd_endo_final1.columns = ['变异基因', '外显子', '氨基酸改变', self.type1, '内分泌药物敏感性', '内分泌药物']
            pd_endo_final2.columns = ['变异基因', '外显子', '氨基酸改变', self.type2, '内分泌药物敏感性', '内分泌药物']
            pd_endo_final1 = pd_endo_final1[~pd_endo_final1['内分泌药物'].isin([''])]
            pd_endo_final2 = pd_endo_final2[~pd_endo_final2['内分泌药物'].isin([''])]
            pd_endo = pd.merge(pd_endo_final1, pd_endo_final2, how='outer', on=['变异基因', '外显子', '氨基酸改变', '内分泌药物敏感性', '内分泌药物'])
            pd_endo['检测结果'] = [f"{pd_endo['外显子'][x]} {pd_endo['氨基酸改变'][x]}" for x in range(len(pd_endo))]
            pd_endo['类型'] = '内分泌治疗药物相关突变'
            pd_endo = pd_endo[['类型', '变异基因', '检测结果', self.type1, self.type2, '内分泌药物敏感性', '内分泌药物']]
            pd_endo = pd_endo.fillna(value='-')

            pd_ARV7_final1 = pd_ARV7_final1[header3].drop_duplicates()
            pd_ARV7_final2 = pd_ARV7_final2[header4].drop_duplicates()
            pd_ARV7_final1 = pd_ARV7_final1[pd_ARV7_final1['基因'] == 'AR']
            pd_ARV7_final1 = pd_ARV7_final1[pd_ARV7_final1['拷贝数'] >= 4]
            pd_ARV7_temp2 = pd_ARV7_final2[pd_ARV7_final2['基因'] == 'AR']
            pd_ARV7_final2 = pd_ARV7_temp2.copy()
            pd_ARV7_final2 = pd_ARV7_final2[pd_ARV7_final2['拷贝数'] >= 4]

            pd_ARV7_final1.columns = ['变异基因', self.type1, '内分泌药物敏感性', '内分泌药物']
            pd_ARV7_final2.columns = ['变异基因', self.type2, '内分泌药物敏感性', '内分泌药物']
            pd_ARV7 = pd.merge(pd_ARV7_final1, pd_ARV7_final2, how='outer', on=['变异基因', '内分泌药物敏感性', '内分泌药物'])
            pd_ARV7['检测结果'] = '拷贝数扩增'
            pd_ARV7['类型'] = '内分泌治疗药物相关突变'
            pd_ARV7 = pd_ARV7[['类型', '变异基因', '检测结果', self.type1, self.type2, '内分泌药物敏感性', '内分泌药物']]
            pd_ARV7[self.type1] = [f'cn={x}' if str(x) != 'nan' else '-' for x in list(pd_ARV7[self.type1])]
            pd_ARV7[self.type2] = [f'cn={x}' if str(x) != 'nan' else '-' for x in list(pd_ARV7[self.type2])]
            pd_ARV7 = pd_ARV7.fillna(value='-')

            pd_cas = pd.read_excel(self.castration)
            pd_HSD3B1 = pd.DataFrame(columns=['类型', '变异基因', '检测结果', '内分泌药物敏感性', '内分泌药物'])
            pd_HSD3B1['类型'] = ['内分泌治疗药物相关突变']
            pd_HSD3B1['变异基因'] = ['HSD3B1']
            pd_HSD3B1_type = '/'.join(pd_cas['检测结果'][0].split(' ')) if not pd_cas['检测结果'].empty else '-/-'
            pd_HSD3B1['检测结果'] = [f'1245 {pd_HSD3B1_type} 型']
            pd_HSD3B1[self.type1] = ['-']
            pd_HSD3B1[self.type2] = ['-']
            pd_HSD3B1['内分泌药物敏感性'] = [pd_cas['结果参考'][0]]
            pd_HSD3B1['内分泌药物'] = ['去势治疗疗效']

            pd_endo_final = pd.concat([pd_HSD3B1, pd_endo, pd_ARV7])
            pd_endo_final['敏感性预测'] = pd_endo_final['内分泌药物敏感性']
            pd_endo_final['药物']= pd_endo_final['内分泌药物']
            pd_endo_final.to_csv(f'{self.outdir}/内分泌治疗药物相关突变.txt', sep='\t', index=False)

        # cnv
        header1 = ['基因', '拷贝数', '敏感性预测', '药物', '疾病关键字']
        header2 = ['基因', '拷贝数', '敏感性预测', '药物', '疾病关键字']
        if self.single_flag:
            pd_cnv_final2 = pd.DataFrame(columns=header2)

        pd_cnv_final1 = pd_cnv_final1[header1].drop_duplicates()
        pd_cnv_final2 = pd_cnv_final2[header2].drop_duplicates()
        pd_cnv_final1.columns = ['变异基因', self.type1, '敏感性预测', '药物', '疾病关键字']
        pd_cnv_final2.columns = ['变异基因', self.type2, '敏感性预测', '药物', '疾病关键字']
        pd_cnv = pd.merge(pd_cnv_final1, pd_cnv_final2, on=['变异基因', '敏感性预测', '药物', '疾病关键字'], how='outer')
        pd_cnv = pd_cnv[~(pd_cnv['敏感性预测'] == '')].reset_index(drop=True)
        pd_cnv['检测结果'] = '拷贝数扩增'
        pd_cnv['类型'] = ['其它重要基因相关突变' if x in ['', ' '] else '靶向药物相关突变' for x in list(pd_cnv['药物'])]
        pd_cnv_target = pd_cnv[pd_cnv['类型'] == '靶向药物相关突变']
        pd_cnv_other = pd_cnv[pd_cnv['类型'] == '其它重要基因相关突变']
        pd_cnv_other = pd_cnv_other[pd_cnv_other['疾病关键字'] == self.cancer_flag]
        pd_cnv = pd.concat([pd_cnv_target, pd_cnv_other])
        pd_cnv = pd_cnv[['类型', '变异基因', '检测结果', self.type1, self.type2, '敏感性预测', '药物']]
        pd_cnv[self.type1] = [f'cn={x}' if str(x) != 'nan' else '-' for x in list(pd_cnv[self.type1])]
        pd_cnv[self.type2] = [f'cn={x}' if str(x) != 'nan' else '-' for x in list(pd_cnv[self.type2])]
        pd_cnv = pd_cnv.fillna(value='-')
        pd_cnv.replace('', np.NaN, inplace=True)
        pd_cnv = pd_cnv[~pd_cnv['药物'].isna()]

        # fus
        header1 = ['基因', 'fusion_gene_name', 'value', '敏感性预测', '药物', '疾病关键字']
        header2 = ['基因', 'fusion_gene_name', 'value', '敏感性预测', '药物', '疾病关键字']
        pd_fus_final1 = pd_fus_final1[header1].drop_duplicates()

        if os.path.isfile(self.rna_fusion_file):
            pd_rna_fusion = self.filter_rna(if_do_targeting_site_drug=True)
            pd_rna_fusion = pd_rna_fusion[header1].drop_duplicates()
            pd_fus_final1 = pd.concat([pd_fus_final1, pd_rna_fusion], axis=0)
            pd_fus_final1 = pd_fus_final1.groupby(['基因', 'fusion_gene_name', '敏感性预测', '药物', '疾病关键字'])['value'].apply(
                lambda x: Tools.cluster_columns(x, sep=';'))
            pd_fus_final1 = pd_fus_final1.reset_index()
            list_value_format = []
            for i in list(pd_fus_final1['value']):
                if ';' in i:
                    value_info = str(i).split(';')
                    value_rna = [x for x in value_info if 'FFPM' in x][0]
                    value_dna = [x for x in value_info if '%' in x][0]
                    list_value_format.append(f"{value_dna}(DNA);{value_rna}(RNA)")
                else:
                    list_value_format.append(i)
            pd_fus_final1['value'] = list_value_format
            pd_fus_final1 = pd_fus_final1[header1]

        if self.single_flag:
            pd_fus_final2 = pd.DataFrame(columns=header2)
        pd_fus_final2 = pd_fus_final2[header2].drop_duplicates()
        pd_fus_final1.columns = ['变异基因', '检测结果', self.type1, '敏感性预测', '药物', '疾病关键字']
        pd_fus_final2.columns = ['变异基因', '检测结果', self.type2, '敏感性预测', '药物', '疾病关键字']
        pd_fus = pd.merge(pd_fus_final1, pd_fus_final2, how='outer', on=['变异基因', '检测结果', '敏感性预测', '药物', '疾病关键字'])
        pd_fus = pd_fus[~(pd_fus['敏感性预测'] == '')].reset_index(drop=True)
        pd_fus['检测结果'] = [f'{x} 融合' for x in list(pd_fus['检测结果'])]
        pd_fus['类型'] = ['其它重要基因相关突变' if x in ['', ' '] else '靶向药物相关突变' for x in list(pd_fus['药物'])]
        pd_fus_target = pd_fus[pd_fus['类型'] == '靶向药物相关突变']
        pd_fus_other = pd_fus[pd_fus['类型'] == '其它重要基因相关突变']
        pd_fus_other = pd_fus_other[pd_fus_other['疾病关键字'] == self.cancer_flag]
        pd_fus = pd.concat([pd_fus_target, pd_fus_other])
        pd_fus = pd_fus[['类型', '变异基因', '检测结果', self.type1, self.type2, '敏感性预测', '药物']]
        pd_fus = pd_fus.fillna(value='-')

        if '前列腺' in self.argv.cancer and self.project_id not in ['NCPM5', 'NCPM6', 'NCPM7', 'NCP393', 'NCP393c' ]:
            out_unit_1 = pd.concat([pd_endo_final, pd_mut, pd_cnv, pd_fus], axis=0, ignore_index=True)
        else:
            out_unit_1 = pd.concat([pd_mut, pd_cnv, pd_fus], axis=0, ignore_index=True)

        if not out_unit_1.empty:
            out_unit_1 = out_unit_1.drop_duplicates()
            out_unit_1 = out_unit_1.groupby(['类型', '变异基因', '检测结果', self.type1, self.type2,
                                             '敏感性预测'])['药物'].apply(lambda x: Tools.cluster_columns(x, sep=','))
            out_unit_1 = out_unit_1.reset_index()

            if '前列腺' in self.argv.cancer and self.project_id not in ['NCPM5', 'NCPM6', 'NCPM7', 'NCP393', 'NCP393c']:
                out_unit_1['临床意义'] = [f"提示可能对{out_unit_1['药物'][x]}{out_unit_1['敏感性预测'][x]}" if (out_unit_1['类型'][x] == '靶向药物相关突变' or out_unit_1['类型'][x] == '内分泌治疗药物相关突变') else f"{out_unit_1['敏感性预测'][x]}" for x in range(len(out_unit_1))]
            else:
                out_unit_1['临床意义'] = [f"提示可能对{out_unit_1['药物'][x]}{out_unit_1['敏感性预测'][x]}"
                                      if out_unit_1['类型'][x] == '靶向药物相关突变' else f"{out_unit_1['敏感性预测'][x]}"
                                      for x in range(len(out_unit_1))]

            out_unit_1 = out_unit_1.drop(['敏感性预测', '药物'], axis=1)

            out_unit_1 = out_unit_1.groupby(['类型', '变异基因', '检测结果', self.type1, self.type2])['临床意义'].apply(
                lambda x: Tools.cluster_columns(x, sep=';'))
            out_unit_1 = out_unit_1.reset_index()

            sort_df_list = []
            if '前列腺' in self.argv.cancer and self.project_id not in ['NCPM5', 'NCPM6', 'NCPM7', 'NCP393', 'NCP393c' ]:
                tag_list = ['内分泌治疗药物相关突变', '靶向药物相关突变', '其它重要基因相关突变']
            else:
                tag_list = ['靶向药物相关突变', '其它重要基因相关突变']
            for tag in tag_list:
                df_tag = out_unit_1[out_unit_1['类型'] == tag]
                df_tag_mut = df_tag[df_tag['检测结果'].str.contains('exon|p.|IVS')]
                df_tag_cnv = df_tag[df_tag['检测结果'].str.contains('拷贝数扩增')]
                df_tag_fus = df_tag[df_tag['检测结果'].str.contains('融合')]
                df_tag_adt = df_tag[df_tag['检测结果'].str.contains('1245')] # add by rhd

                # 这里没有使用 self.type2 进行排序是因为 单样本分析时候 self.type2列基本为空列
                df_tag_mut = df_tag_mut.sort_values(by=['变异基因', self.type1], ascending=False)
                df_tag_cnv = df_tag_cnv.sort_values(by=['变异基因', self.type1], ascending=False)
                df_tag_fus = df_tag_fus.sort_values(by=['变异基因', self.type1], ascending=False)
                df_tag_adt = df_tag_adt.sort_values(by=['变异基因', self.type1], ascending=False)# rhd

                sort_df_list.extend([df_tag_mut, df_tag_cnv, df_tag_fus, df_tag_adt])

            pd_unit_list = pd.concat(sort_df_list, axis=0)
        else:
            pd_unit_list = pd.DataFrame(columns=['类型', '变异基因', '检测结果', self.type1, self.type2, '临床意义'])

        if self.single_flag:
            pd_unit_list = pd_unit_list[['类型', '变异基因', '检测结果', self.type1, '临床意义']]

        # 20210908肠癌新增(KRAS/NRAS/BRAF三基因阴性时在肠癌中需展示)
        special_gene_frame = pd_unit_list[((pd_unit_list['变异基因'] == 'KRAS') | (pd_unit_list['变异基因'] == 'NRAS') | (
                pd_unit_list['变异基因'] == 'BRAF'))]

        if self.cancer_flag == '肠' and special_gene_frame.empty:
            special_info = list()
            for gene in ['KRAS', 'BRAF', 'NRAS']:
                if self.single_flag:
                    special_info.append(['靶向药物相关突变', gene, '野生型', '-', '提示患者对于西妥昔单抗、帕尼单抗敏感'])
                else:
                    special_info.append(['靶向药物相关突变', gene, '野生型', '-', '-', '提示患者对于西妥昔单抗、帕尼单抗敏感'])
            pd_unit_list = pd_unit_list.append(pd.DataFrame(special_info, columns=list(pd_unit_list)), ignore_index=True)

        pd_unit_list.to_csv(f'{self.outdir}/检测结果汇总.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/检测结果汇总.txt')

    def unit_1_2(self):
        """ 1.2 免疫治疗相关检测结果 """
        """ 基础过滤 """

        # (1) tmb 获取
        pd_mut1 = pd.read_csv(self.mutation_file1, sep='\t', low_memory=False)
        pd_mut1['结论'] = pd_mut1.apply(lambda x: Tools.get_conclusion(x, self.product, panel_id=self.tumor_id1, wes_id=self.tumor_wes_id1), axis=1)
        pd_mut1['突变率'] = pd_mut1.apply(lambda x: Tools.get_vaf(x, self.product, panel_id=self.tumor_id1, wes_id=self.tumor_wes_id1), axis=1)
        tmb_value1 = Tools.tmb_process(pd_mut1, self.tumor_id1, self.outdir, self.project_id)
        tmb_value1 = round(float(tmb_value1), 2)
        tmb_level1 = 'High' if tmb_value1 >= 10 else 'Low'

        if self.single_flag:
            tmb_value2 = '-'
            tmb_level2 = '-'
            pd_unit_3_1_1 = pd.DataFrame([['肿瘤突变负荷(Tumor mutation Burden, TMB)', tmb_value1, tmb_level1]],
                                         columns=['检测指标', '检测结果(Muts/Mb)', '分级'])
        else:
            pd_mut2 = pd.read_csv(self.mutation_file2, sep='\t', low_memory=False)
            pd_mut2['结论'] = pd_mut2.apply(
                lambda x: Tools.get_conclusion(x, self.product, panel_id=self.tumor_id2, wes_id=self.tumor_wes_id2),
                axis=1)
            pd_mut2['突变率'] = pd_mut2.apply(
                lambda x: Tools.get_vaf(x, self.product, panel_id=self.tumor_id2, wes_id=self.tumor_wes_id2), axis=1)
            tmb_value2 = Tools.tmb_process(pd_mut2, self.tumor_id2, self.outdir, self.project_id)
            tmb_value2 = round(float(tmb_value2), 2)
            tmb_level2 = 'High' if tmb_value2 >= 10 else 'Low'

            pd_unit_3_1_1 = pd.DataFrame([[self.type1, tmb_value1, tmb_level1], [self.type2, tmb_value2, tmb_level2]],
                                         columns=['检测指标', '检测结果(Muts/Mb)', '分级'])
        pd_out_tmb = pd.DataFrame(
            [['肿瘤突变负荷(TMB，Muts/Mb)', '-', f'{tmb_value1}({tmb_level1})', f'{tmb_value2}({tmb_level2})']],
            columns=['相关性', '变异', self.type1, self.type2])
        pd_unit_3_1_1.to_csv(f'{self.outdir}/肿瘤突变负荷状态分析结果.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/肿瘤突变负荷状态分析结果.txt')

        # (2) MSI
        pd_msi1 = pd.read_csv(self.msi_summary1, sep='\t')
        msi1 = '-' if len(list(pd_msi1['MSI_Status'])) == 0 else list(pd_msi1['MSI_Status'])[0]
        dict_msi_str1 = {'MSS': '微卫星稳定(MSS)',
                         'MSI-H': '微卫星高度不稳定(MSI-H)',
                         'MSI-L': '微卫星低度不稳定(MSI-L)',
                         'MSI': '微卫星不稳定(MSI)'}
        msi_str1 = dict_msi_str1[msi1]

        if self.single_flag:
            msi2 = '-'
            msi_str2 = '-'
            pd_unit_3_1_2 = pd.DataFrame(['微卫星(MSS/MSI)状态分析', msi_str1], index=['检测指标', '检测结果']).T
        else:
            pd_msi2 = pd.read_csv(self.msi_summary2, sep='\t')
            msi2 = '-' if len(list(pd_msi2['MSI_Status'])) == 0 else list(pd_msi2['MSI_Status'])[0]
            msi_str2 = '微卫星稳定(MSS)' if msi2 == 'MSS' else '微卫星高度不稳定(MSI-H)' if msi2 == 'MSI-H' else '微卫星低度不稳定(MSI-L)'
            pd_unit_3_1_2 = pd.DataFrame([[self.type1, msi_str1], [self.type2, msi_str2]], columns=['检测指标', '检测结果'])

        pd_out_msi = pd.DataFrame([['微卫星状态分析(MSS/MSI)', '-', msi_str1, msi_str2]], columns=['相关性', '变异', self.type1, self.type2])

        pd_unit_3_1_2.to_csv(f'{self.outdir}/微卫星状态分析结果.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/微卫星状态分析结果.txt')

        # (3) PD-L1 : 外包检测，lims抓取
        # 阴性/阳性的判断规则为：如果TPS 0%,且CPS 0，则为阴性；否则为阳性，即TPS ≥1% ,或者CPS≥1
        # if os.path.isfile(self.argv.PD_L1_json):
        print(self.PD_L1_json['TPS']['info'])
        print(self.PD_L1_json['CPS']['info'])
        TPS = self.PD_L1_json['TPS']['value']
        CPS = self.PD_L1_json['CPS']['value']
        try:
            if (TPS == '<1' or float(TPS) == 0) and (CPS == '<1' or float(CPS) == 0):
                pd_l1_str = '阴性'
            else:
                pd_l1_str = f'阳性(TPS {TPS}%，CPS {CPS})'

            pd_unit_3_1_3 = pd.DataFrame([['Dako 22C3', f'{TPS}%', CPS, pd_l1_str.split('(')[0]]],
                                         columns=['抗体类型', '肿瘤细胞阳性比例分数(TPS)', '综合阳性分数(CPS)',
                                                  'PD-L1免疫组化检测结果'])
        except:
            pd_l1_str = '-'
            pd_unit_3_1_3 = pd.DataFrame(columns=['抗体类型', '肿瘤细胞阳性比例分数（TPS）', '综合阳性分数（CPS）'])
        # pd_l2_str
        pd_out_pd_l1 = pd.DataFrame([['PD-L1免疫组化检测结果', '-', pd_l1_str, pd_l1_str]],
                                    columns=['相关性', '变异', self.type1, self.type2])
        pd_unit_3_1_3.to_csv(f'{self.outdir}/免疫组化检测结果.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/免疫组化检测结果.txt')

        # (5) hla 分型
        pd_hla = pd.read_csv(self.hla_file, sep='\t')
        pd_hla = pd_hla['GeneType'].str.split('_', expand=True)
        pd_hla[0] = pd_hla[0].str[0:10]
        pd_hla[1] = pd_hla[1].str[0:10]
        row_num = pd_hla[~(pd_hla[0] == pd_hla[1])].shape[0]
        if row_num == 3:
            hla_str = '全杂合'
        elif 0 < row_num < 3:
            hla_str = '部分纯合'
        else:
            hla_str = '全纯合'
        pd_out_hla = pd.DataFrame([['HLA分型结果', '-', hla_str, hla_str]], columns=['相关性', '变异', self.type1, self.type2])

        # 突变相关性
        pd_mut1, pd_mut2 = self.filter_mut(if_do_targeting_site_drug=True)
        if self.single_flag:
            pd_mut2 = pd.DataFrame(columns=['基因', '外显子', '氨基酸改变', '突变率'])

        pd_mut1 = pd_mut1[['基因', '外显子', '氨基酸改变', '突变率']]
        pd_mut2 = pd_mut2[['基因', '外显子', '氨基酸改变', '突变率']]
        pd_mut1.columns = ['基因', '外显子', '氨基酸改变', self.type1]
        pd_mut2.columns = ['基因', '外显子', '氨基酸改变', self.type2]

        pd_mut1['变异'] = ''
        pd_mut2['变异'] = ''

        pd_mut1['变异2'] = ''
        pd_mut2['变异2'] = ''

        pd_mut1['变异'] = pd_mut1['基因'].str.cat([pd_mut1['外显子'], pd_mut1['氨基酸改变']], sep=' ')
        pd_mut1['变异2'] = pd_mut1['外显子'].str.cat([pd_mut1['氨基酸改变']], sep=' ')

        pd_mut2['变异'] = pd_mut2['基因'].str.cat([pd_mut2['外显子'], pd_mut2['氨基酸改变']], sep=' ')
        pd_mut2['变异2'] = pd_mut2['外显子'].str.cat([pd_mut2['氨基酸改变']], sep=' ')

        pd_mut1 = pd_mut1[['基因', '变异', '变异2', self.type1]]
        pd_mut2 = pd_mut2[['基因', '变异', '变异2', self.type2]]

        pd_mut1 = pd_mut1.drop_duplicates()
        pd_mut2 = pd_mut2.drop_duplicates()

        if pd_mut1.empty and pd_mut2.empty:
            pd_mut = pd.DataFrame(columns=['基因', '变异', '变异2', self.type1, self.type2])
        else:
            pd_mut = pd.merge(pd_mut1, pd_mut2, on=['基因', '变异', '变异2'], how='outer')
            pd_mut = pd_mut.fillna(value='-')

        pd_cnv1, pd_cnv2 = self.filter_cnv(if_do_targeting_site_drug=True)
        if self.single_flag:
            pd_cnv2 = pd.DataFrame(columns=['基因', '变异（缺失/扩增）', '拷贝数'])
        pd_cnv1 = pd_cnv1[['基因', '变异（缺失/扩增）', '拷贝数']]
        pd_cnv2 = pd_cnv2[['基因', '变异（缺失/扩增）', '拷贝数']]

        pd_cnv1.columns = ['基因', '变异（缺失/扩增）', self.type1]
        pd_cnv2.columns = ['基因', '变异（缺失/扩增）', self.type2]

        pd_cnv1['变异'] = [f'{x} 扩增' for x in pd_cnv1['基因']]
        pd_cnv2['变异'] = [f'{x} 扩增' for x in pd_cnv2['基因']]

        pd_cnv1['变异2'] = '扩增'
        pd_cnv2['变异2'] = '扩增'

        pd_cnv1 = pd_cnv1[['基因', '变异', '变异2', self.type1]]
        pd_cnv2 = pd_cnv2[['基因', '变异', '变异2', self.type2]]

        pd_cnv1 = pd_cnv1[pd_cnv1['基因'].isin(['MDM2', 'MDM4', 'FGF3', 'FGF4', 'FGF19', 'CCND1'])]
        pd_cnv2 = pd_cnv2[pd_cnv2['基因'].isin(['MDM2', 'MDM4', 'FGF3', 'FGF4', 'FGF19', 'CCND1'])]

        pd_cnv1 = pd_cnv1.drop_duplicates()
        pd_cnv2 = pd_cnv2.drop_duplicates()
        if pd_cnv1.empty and pd_cnv2.empty:
            pd_cnv = pd.DataFrame(columns=['基因', '变异', '变异2', self.type1, self.type2])
        else:
            pd_cnv = pd.merge(pd_cnv1, pd_cnv2, on=['基因', '变异', '变异2'], how='outer')
            pd_cnv = pd_cnv.fillna(value='-')

        pd_fusion1, pd_fusion2 = self.filter_fus(if_do_targeting_site_drug=True)
        # '基因', 'value', 'fusion_gene_name'
        if self.single_flag:
            pd_fusion2 = pd.DataFrame(columns=['基因', 'fusion_gene_name', 'value'])
        pd_fusion1 = pd_fusion1[['基因', 'fusion_gene_name', 'value']]
        pd_fusion2 = pd_fusion2[['基因', 'fusion_gene_name', 'value']]

        pd_fusion1.columns = ['基因', 'fusion_gene_name', self.type1]
        pd_fusion2.columns = ['基因', 'fusion_gene_name', self.type2]

        pd_fusion1['变异'] = [f"{list(pd_fusion1['基因'])[x]} {list(pd_fusion1['fusion_gene_name'])[x]}" for x in range(len(pd_fusion1))]
        pd_fusion2['变异'] = [f"{list(pd_fusion2['基因'])[x]} {list(pd_fusion2['fusion_gene_name'])[x]}" for x in range(len(pd_fusion2))]

        pd_fusion1['变异2'] = pd_fusion1['fusion_gene_name']
        pd_fusion2['变异2'] = pd_fusion2['fusion_gene_name']

        pd_fusion1 = pd_fusion1[['基因', '变异', '变异2', self.type1]]
        pd_fusion2 = pd_fusion2[['基因', '变异', '变异2', self.type2]]

        # 重排/融合：此模块经医学内部分析，只有ALK基因需要关注，在负相关sheet
        pd_fusion1 = pd_fusion1[pd_fusion1['基因'] == 'ALK']
        pd_fusion2 = pd_fusion2[pd_fusion2['基因'] == 'ALK']

        pd_fusion1 = pd_fusion1.drop_duplicates()
        pd_fusion2 = pd_fusion2.drop_duplicates()
        # ['基因', '变异', '变异2', 'tissue']
        if os.path.isfile(self.rna_fusion_file):
            header = list(pd_fusion1)
            pd_rna_fusion = self.filter_rna(if_do_targeting_site_drug=True)
            pd_rna_fusion['变异'] = pd_rna_fusion['基因'] + ' ' + pd_rna_fusion['fusion_gene_name']
            pd_rna_fusion['变异2'] = pd_rna_fusion['fusion_gene_name']
            pd_rna_fusion[self.type1] = pd_rna_fusion['value']
            pd_rna_fusion = pd_rna_fusion[header].drop_duplicates()
            pd_fusion1 = pd.concat([pd_fusion1, pd_rna_fusion], axis=0)
            pd_fusion1 = pd_fusion1[header]

        if pd_fusion1.empty and pd_fusion2.empty:
            pd_fusion = pd.DataFrame(columns=['基因', '变异', '变异2', self.type1, self.type2])
        else:
            pd_fusion = pd.merge(pd_fusion1, pd_fusion2, on=['基因', '变异', '变异2'], how='outer')
            pd_fusion = pd_fusion.fillna(value='-')

        "db_immunotherapy_gene : '基因' '检测结果' '临床意义' '相关癌种' '参考依据'"
        db_immunotherapy_gene = pd.read_excel(self.solid_tumor_immunotherapy_gene, sheet_name=['正相关', '负相关', '超进展'])
        db_positive_correlation = db_immunotherapy_gene['正相关'].ffill()
        db_negative_correlation = db_immunotherapy_gene['负相关'].ffill()
        db_super_progress = db_immunotherapy_gene['超进展'].ffill()

        db_positive_correlation['检测结果'] = db_positive_correlation['检测结果'].fillna(value='未检出')
        db_negative_correlation['检测结果'] = db_negative_correlation['检测结果'].fillna(value='未检出')
        db_super_progress['检测结果'] = db_super_progress['检测结果'].fillna(value='未检出')

        db_positive_correlation['相关性'] = '免疫疗效相关性基因突变:正相关'
        db_negative_correlation['相关性'] = '免疫疗效相关性基因突变:负相关'
        db_super_progress['相关性'] = '免疫疗效相关性基因突变:超进展'

        db_merge_immunotherapy_gene = pd.concat([db_positive_correlation, db_negative_correlation, db_super_progress], axis=0)
        pd_mut_merge = pd.merge(db_merge_immunotherapy_gene, pd_mut, on='基因', how='inner')
        pd_cnv_merge = pd.merge(db_merge_immunotherapy_gene, pd_cnv, on='基因', how='inner')
        pd_fus_merge = pd.merge(db_merge_immunotherapy_gene, pd_fusion, on='基因', how='inner')
        pd_fus_merge = pd_fus_merge.fillna(value='-')

        pd_immunotherapy = pd.concat([pd_mut_merge, pd_cnv_merge, pd_fus_merge], axis=0)
        pd_positive_correlation = pd_immunotherapy[pd_immunotherapy['相关性'] == '免疫疗效相关性基因突变:正相关']
        pd_negative_correlation = pd_immunotherapy[pd_immunotherapy['相关性'] == '免疫疗效相关性基因突变:负相关']
        pd_super_progress = pd_immunotherapy[pd_immunotherapy['相关性'] == '免疫疗效相关性基因突变:超进展']

        if pd_positive_correlation.empty:
            pd_positive_correlation = pd.DataFrame(columns=['基因', '检测结果', '临床意义', '相关癌种', '参考依据', '相关性', '变异', '变异2'])
        if pd_negative_correlation.empty:
            pd_negative_correlation = pd.DataFrame(columns=['基因', '检测结果', '临床意义', '相关癌种', '参考依据', '相关性', '变异', '变异2'])
        if pd_super_progress.empty:
            pd_super_progress       = pd.DataFrame(columns=['基因', '检测结果', '临床意义', '相关癌种', '参考依据', '相关性', '变异', '变异2'])

        pd_positive_correlation = pd_positive_correlation.drop_duplicates()
        pd_negative_correlation = pd_negative_correlation.drop_duplicates()
        pd_super_progress       = pd_super_progress.drop_duplicates()

        out_positive_correlation = pd.merge(db_positive_correlation, pd_positive_correlation[['基因', '变异', '变异2']],
                                            on='基因', how='left')
        out_negative_correlation = pd.merge(db_negative_correlation, pd_negative_correlation[['基因', '变异', '变异2']],
                                            on='基因', how='left')
        out_super_progress       = pd.merge(db_super_progress, pd_super_progress[['基因', '变异', '变异2']], on='基因', how='left')

        out_positive_correlation = out_positive_correlation[['基因', '变异2', '临床意义', '相关癌种', '参考依据', '相关性']]
        out_negative_correlation = out_negative_correlation[['基因', '变异2', '临床意义', '相关癌种', '参考依据', '相关性']]
        out_super_progress = out_super_progress[['基因', '变异2', '临床意义', '相关癌种', '参考依据', '相关性']]

        out_positive_correlation.columns = ['基因', '检测结果', '临床意义', '相关癌种', '参考依据', '相关性']
        out_negative_correlation.columns = ['基因', '检测结果', '临床意义', '相关癌种', '参考依据', '相关性']
        out_super_progress.columns = ['基因', '检测结果', '临床意义', '相关癌种', '参考依据', '相关性']

        if not out_positive_correlation.empty:
            out_positive_correlation = out_positive_correlation.groupby(
                ['基因', '临床意义', '相关癌种', '参考依据', '相关性'])['检测结果'].apply(Tools.cluster_columns, sep=';')
            out_positive_correlation = out_positive_correlation.reset_index()

        if not out_negative_correlation.empty:
            out_negative_correlation = out_negative_correlation.groupby(
                ['基因', '临床意义', '相关癌种', '参考依据', '相关性'])['检测结果'].apply(Tools.cluster_columns, sep=';')
            out_negative_correlation = out_negative_correlation.reset_index()

        if not out_super_progress.empty:
            out_super_progress = out_super_progress.groupby(
                ['基因', '临床意义', '相关癌种', '参考依据', '相关性'])['检测结果'].apply(Tools.cluster_columns, sep=';')
            out_super_progress = out_super_progress.reset_index()

        out_positive_correlation = out_positive_correlation[['基因', '检测结果', '临床意义', '相关癌种', '参考依据', '相关性']]
        out_negative_correlation = out_negative_correlation[['基因', '检测结果', '临床意义', '相关癌种', '参考依据', '相关性']]
        out_super_progress = out_super_progress[['基因', '检测结果', '临床意义', '相关癌种', '参考依据', '相关性']]

        out_positive_correlation['检测结果'] = ['未检出' if str(x) == 'nan' else str(x) for x in
                                            list(out_positive_correlation['检测结果'])]
        out_negative_correlation['检测结果'] = ['未检出' if str(x) == 'nan' else str(x) for x in
                                            list(out_negative_correlation['检测结果'])]
        out_super_progress['检测结果'] = ['未检出' if str(x) == 'nan' else str(x) for x in
                                      list(out_super_progress['检测结果'])]

        out_positive_correlation.to_csv(f'{self.outdir}/免疫治疗疗效正相关基因检测结果.txt', sep='\t', index=False)
        out_negative_correlation.to_csv(f'{self.outdir}/免疫治疗疗效负相关基因检测结果.txt', sep='\t', index=False)
        out_super_progress.to_csv(f'{self.outdir}/免疫治疗疗效超进展相关基因检测结果.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/免疫治疗疗效正相关基因检测结果.txt')
        self.result_file_list.append(f'{self.outdir}/免疫治疗疗效负相关基因检测结果.txt')
        self.result_file_list.append(f'{self.outdir}/免疫治疗疗效超进展相关基因检测结果.txt')

        boolean_positive_correlation = not pd_positive_correlation.empty
        boolean_negative_correlation = not pd_negative_correlation.empty
        boolean_super_progress       = not pd_super_progress.empty

        if not boolean_positive_correlation:
            pd_positive_correlation = pd.DataFrame(
                {'相关性': ['免疫疗效相关性基因突变:正相关'], '变异': ['未检出'], self.type1: ['-'], self.type2: ['-']})
        else:
            pd_positive_correlation = pd.DataFrame(
                {'相关性': ['免疫疗效相关性基因突变:正相关'],
                 '变异': [';'.join([str(x) for x in list(pd_positive_correlation['变异'])])],
                 self.type1: [';'.join([str(x) for x in list(pd_positive_correlation[self.type1])])],
                 self.type2: [';'.join([str(x) for x in list(pd_positive_correlation[self.type2])])]
                 })
        if not boolean_negative_correlation:
            pd_negative_correlation = pd.DataFrame(
                {'相关性': ['免疫疗效相关性基因突变:负相关'], '变异': ['未检出'], self.type1: ['-'], self.type2: ['-']})
        else:
            pd_negative_correlation = pd.DataFrame(
                {'相关性': ['免疫疗效相关性基因突变:负相关'],
                 '变异': [';'.join([str(x) for x in list(pd_negative_correlation['变异'])])],
                 self.type1: [';'.join([str(x) for x in list(pd_negative_correlation[self.type1])])],
                 self.type2: [';'.join([str(x) for x in list(pd_negative_correlation[self.type2])])]
                 })
        if not boolean_super_progress:
            pd_super_progress       = pd.DataFrame(
                {'相关性': ['免疫疗效相关性基因突变:超进展'],
                 '变异': ['未检出'], self.type1: ['-'], self.type2: ['-']})
        else:
            pd_super_progress       = pd.DataFrame(
                {'相关性': ['免疫疗效相关性基因突变:超进展'],
                 '变异': [';'.join([str(x) for x in list(pd_super_progress['变异'])])],
                 self.type1: [';'.join([str(x) for x in list(pd_super_progress[self.type1])])],
                 self.type2: [';'.join([str(x) for x in list(pd_super_progress[self.type2])])]
                 })

        pd_immunotherapy_final = pd.concat([pd_positive_correlation, pd_negative_correlation, pd_super_progress], axis=0)
        pd_out_immunotherapy = pd_immunotherapy_final[['相关性', '变异', self.type1, self.type2]]
        pd_out_immunotherapy.to_csv(f'{self.outdir}/免疫治疗疗效所有相关基因检测结果.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/免疫治疗疗效所有相关基因检测结果.txt')

        # 结论 conclusion
        list_conclusion_a = []
        if self.single_flag:
            for_list = [[tmb_value1, msi1]]
        else:
            for_list = [[tmb_value1, msi1], [tmb_value2, msi2]]

        if TPS == '' and CPS == '':
            for tmb_value, msi in for_list:
                if (not boolean_negative_correlation) and \
                        (
                                (float(tmb_value) >= 10) or
                                (msi == 'MSI-H') or
                                boolean_positive_correlation
                        ):
                    _conclusion_a = '较好'
                elif (float(tmb_value) < 10) and \
                        (msi in ['MSI-L', 'MSS']) and \
                        (not boolean_positive_correlation and boolean_negative_correlation):
                    _conclusion_a = '较差'
                else:
                    _conclusion_a = '一般'
                list_conclusion_a.append(_conclusion_a)
        else:
            for tmb_value, msi in for_list:
                if (not boolean_negative_correlation) and \
                        ( (float(tmb_value) >= 10) or (msi == 'MSI-H') or (pd_l1_str == '阳性') or boolean_positive_correlation ):
                    _conclusion_a = '较好'
                elif (float(tmb_value) < 10) and \
                        (msi in ['MSI-L', 'MSS']) and \
                        (pd_l1_str == '阳性') and \
                        (not boolean_positive_correlation and boolean_negative_correlation):
                    _conclusion_a = '较差'
                else:
                    _conclusion_a = '一般'

                list_conclusion_a.append(_conclusion_a)

        if '较好' in list_conclusion_a:
            conclusion_a = '较好'
        elif '一般' in list_conclusion_a:
            conclusion_a = '一般'
        elif '较差' in list_conclusion_a:
            conclusion_a = '较差'

        if boolean_super_progress:
            conclusion_b = '警惕可能出现的免疫治疗超进展。'
        else:
            conclusion_b = ''

        conclusion = f'提示患者从免疫治疗获益效果{conclusion_a}，若使用相关药物，请加强随访。{conclusion_b}'

        pd_out_unit_1_2 = pd.concat([pd_out_tmb, pd_out_msi, pd_out_pd_l1, pd_out_immunotherapy, pd_out_hla], axis=0)
        pd_out_unit_1_2['结论'] = conclusion
        if self.single_flag:
            pd_out_unit_1_2 = pd_out_unit_1_2[['相关性', '变异', self.type1, '结论']]
        pd_out_unit_1_2.to_csv(f'{self.outdir}/免疫治疗相关检测结果.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/免疫治疗相关检测结果.txt')

    def unit_1_3(self):
        """ 遗传性突变 雅丽指导双样本只用 组织的 """
        pd_tumor_genetic_predisposed_gene_list = pd.read_csv(self.tumor_genetic_predisposed_gene_list, sep='\t')
        pd_mut = pd.read_csv(self.mutation_file1, sep='\t', low_memory=False, keep_default_na=False)
        # pd_mut = pd_mut.fillna(value='') 不可以进行替换NA  因为 pathogenicity 这块需要提取NA的信息
        pd_mut = pd_mut.apply(lambda x: Tools.rm_exon(x), axis=1)

        pd_mut['结论'] = pd_mut.apply(
            lambda x: Tools.get_conclusion(x, self.product, panel_id=self.tumor_id1, wes_id=self.tumor_wes_id1), axis=1)

        # (1) 纯合、杂合、半合子 判断
        pd_mut['突变状态'] = ''
        pd_mut["突变状态"] = pd_mut.apply(lambda x: Tools.mut_status(x, self.gender, self.normal_id), axis=1)

        # (2) 致病性、可能致病性、NA 判断  去掉其他无意义注释

        pd_mut['pathogenicity'] = ''
        pd_mut["pathogenicity"] = pd_mut.apply(lambda x: Tools.get_pathogenicity(x), axis=1)

        # (3) 根据 (2) 构造 '突变意义描述' 内容
        pd_mut['突变意义描述'] = ''
        pd_mut['突变意义描述'] = '该变异在ClinVar数据库中收录为“' + pd_mut["pathogenicity"] + '”变异，遗传该突变可能增加' + pd_mut['疾病表型名称'] + '的发病风险'
        pd_mut = pd_mut[['结论', '突变级别', '突变类型', 'pathogenicity', '基因', '核酸改变', '氨基酸改变', '突变状态', '突变意义描述']]
        # 检不出结果可以解开下面两行这样测试
        # pd_mut['突变级别'] = '1'
        # pd_mut['pathogenicity'] = '致病'
        # pd_mut.to_csv('遗传性肿瘤基因检测结果temp.txt', sep='\t', index=False)

        # (4) 进行'遗传性突变'的过滤
        pd_mut = pd_mut[
            (pd_mut['结论'] == '遗传性突变') &
            pd_mut['突变级别'].isin(['1', '2']) &
            ~pd_mut['突变类型'].isin(['非编码区突变', '同义突变']) &
            (pd_mut['pathogenicity'].isin(['致病', '可能致病', '-']))
            ]

        # (5) 还需要进行基因过滤，这里使用merge进行过滤
        # 20211105 rhd 增加肾癌81panel 32个遗传性肾癌基因过滤
        pd_merge = pd.merge(pd_mut, pd_tumor_genetic_predisposed_gene_list, on='基因', how='inner')
        panel39_code = ['NCPC1', 'NCPC2', 'NCPC3']
        if self.project_id in panel39_code:
            panel39_gene_frame = pd.read_excel(panel39_germline_genes, sheet_name='panel39')
            pd_merge = pd.merge(pd_merge, panel39_gene_frame)
        # (6) 提取此模块需要的列
        pd_merge = pd_merge[['基因', '核酸改变', '氨基酸改变', '突变状态', '突变意义描述']]

        def filter_germline_gene(pd_merge):
            little_panel = ['NCPK1', 'NCPK2', 'NCPK3']
            sheet_dict = {'肾': '肾32遗传性基因'}
            if self.project_id in little_panel:
                pd_panel_gene = pd.read_excel(self.solid_gene_list, sheet_name=sheet_dict[self.cancer], header=None)# 库表 癌症为小类
                filter_gene_list = [x[0] for x in pd_panel_gene.values.tolist()]
                pd_merge = pd_merge[pd_merge['基因'].isin(filter_gene_list)]
            return pd_merge
        pd_merge = filter_germline_gene(pd_merge)

        # (7) 进行列的合并
        pd_merge['遗传性肿瘤基因检测结果'] = ''
        pd_merge['遗传性肿瘤基因检测结果'] = pd_merge['基因'].str.cat([pd_merge['核酸改变'],
                                                          pd_merge['氨基酸改变'].fillna(value='').replace('NA', '-'), pd_merge['突变状态']], sep=' ')
        # pd_merge['遗传性肿瘤基因检测结果'] = pd_merge['基因'].fillna(value='-') + ' ' + pd_merge['氨基酸改变'].fillna(value='-')
        pd_merge = pd_merge[['遗传性肿瘤基因检测结果', '突变意义描述']]
        pd_merge.to_csv(f'{self.outdir}/遗传性肿瘤基因检测结果.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/遗传性肿瘤基因检测结果.txt')

    def unit_nccn_glioma(self):
        """ NCCN """
        """ 脑胶质瘤没有双样本 所以这里过滤函数返回两个 只用了第一个返回的内容"""
        pd_mutation_no_filter = pd.read_csv(self.mutation_file1, sep='\t', low_memory=False)
        pd_mutation_no_filter['突变率'] = pd_mutation_no_filter.apply(
            lambda x: Tools.get_vaf(x, self.product, panel_id=self.tumor_id1, wes_id=self.tumor_wes_id1), axis=1)
        pd_mutation_no_filter['突变深度'] = pd_mutation_no_filter.apply(
            lambda x: Tools.get_depth(x, self.product, panel_id=self.tumor_id1, wes_id=self.tumor_wes_id1), axis=1)
        pd_mutation_no_filter = pd_mutation_no_filter.fillna(value='')
        pd_mutation_no_filter['突变率'] = pd_mutation_no_filter['突变率'].map(lambda x: format(x, '.1%'))
        pd_mutation1, pd_mutation2 = self.filter_mut(if_do_targeting_site_drug=False)
        pd_mutation1 = pd_mutation1.fillna(value='')
        pd_cnv1, pd_cnv2 = self.filter_cnv(if_do_targeting_site_drug=False)
        pd_fusion1, pd_fusion2 = self.filter_fus(if_do_targeting_site_drug=False)

        pd_mutation_no_filter['检测结果'] = ''
        pd_mutation1['检测结果'] = ''
        pd_cnv1['检测结果'] = ''
        pd_fusion1['检测结果'] = ''

        def mutation_convert(x):
            # 如果p.为空（即剪切突变），则其位置输出c._突变丰度
            if x['氨基酸改变'] == '':
                return x['核酸改变'] + ' ' + str(x['突变率'])
            else:
                return x['氨基酸改变'] + ' ' + str(x['突变率'])

        pd_mutation_no_filter['检测结果'] = pd_mutation_no_filter.apply(lambda x:  mutation_convert(x), axis=1)
        pd_mutation1['检测结果'] = pd_mutation1.apply(lambda x: mutation_convert(x), axis=1)
        pd_cnv1['检测结果'] = 'cn=' + pd_cnv1['拷贝数'].map(str)
        pd_fusion1['检测结果'] = pd_fusion1['fusion_gene_name'] + ' ' + pd_fusion1['value']

        pd_nccn_glioma = pd.read_excel(self.solid_tumor_NCCN_gene, sheet_name=self.cancer_flag)
        pd_nccn_glioma = pd_nccn_glioma.fillna(method='ffill')
        pd_nccn_glioma.set_index(["基因"], inplace=True)
        pd_nccn_glioma['检测结果'] = '未检出'
        dict_nccn_glioma = pd_nccn_glioma.to_dict()
        dict_nccn_result = dict_nccn_glioma['检测结果']
        dict_nccn_type = dict_nccn_glioma['变异类型']

        for each_gene in dict_nccn_result.keys():

            if each_gene == '1p19q':
                # (1) 1p19q
                result_1p19q = open(self.cnv_1p19q_1).read().strip().split(':')[1]
                dict_nccn_result[each_gene] = result_1p19q
            elif each_gene == 'MGMT':
                # (2) MGMT 等待外送结果
                dict_nccn_result[each_gene] = '-'
            elif each_gene == 'EGFR':
                # (3) EGFR
                # 此为关注该基因的扩增，需生信分析该基因的扩增情况，检测结果栏呈现“检出”或者“未检出”
                pd_cnv_EGFR = pd_cnv1[pd_cnv1['基因'] == 'EGFR']
                dict_nccn_result[each_gene] = '未检出' if pd_cnv_EGFR.empty else '检出'
            # dict_nccn_result[each_gene] = '未检出' if pd_cnv_EGFR.empty else ';'.join(list(pd_cnv_EGFR['检测结果']))

            elif each_gene == 'H3F3A':
                # (4) H3F3A * merge 有就行 只看这个
                # 此基因只匹配关注p.Lys28Met（p.L28M），其他变异不纳入，检测结果栏呈现“检出”或者“未检出”
                pd_mutation_H3F3A = pd_mutation_no_filter[(pd_mutation_no_filter['基因'] == 'H3F3A') &
                                                          (pd_mutation_no_filter['氨基酸改变'] == 'p.L28M')]
                dict_nccn_result[each_gene] = '未检出' if pd_mutation_H3F3A.empty else '检出'

            elif each_gene == 'TERT':
                # (5) TERT * merge 有就行
                # 此基因发生非编码区突变时，只报出c.-146G>A位置（注意：旧模板此处写成了p.-146G>A是错误的）；
                # 如检出编码区突变，则只报出“结论中”为“体细胞较可靠突变,可能有害、体细胞突变,可能有害”的有害突变。
                pd_mutation_no_filter_TERT = pd_mutation_no_filter[
                    (pd_mutation_no_filter['突变类型'] == '非编码区突变') &
                    (pd_mutation_no_filter['基因'] == 'TERT')
                    ]

                if self.product == 'panel':
                    mut_depth = 20
                else:
                    mut_depth = 8

                pd_mutation_no_filter_TERT = pd_mutation_no_filter_TERT[
                    (pd_mutation_no_filter_TERT["突变深度"] >= mut_depth) &
                    (pd_mutation_no_filter_TERT["总深度(对照,%s)" % self.normal_id] > 9) &
                    (pd_mutation_no_filter_TERT["突变深度(对照,%s)" % self.normal_id] < 4) &
                    (pd_mutation_no_filter_TERT["突变率(对照,%s)" % self.normal_id] < 0.03)
                    ]

                pd_mutation_TERT = pd_mutation1[(pd_mutation1['基因'] == 'TERT')]
                pd_mutation_merge_TERT = pd.concat([pd_mutation_no_filter_TERT, pd_mutation_TERT], axis=0)
                dict_nccn_result[each_gene] = '未检出' if pd_mutation_merge_TERT.empty else ';'.join(list(pd_mutation_merge_TERT['检测结果']))

            elif each_gene in ['RELA', 'NTRK1', 'NTRK2', 'NTRK3', 'NAB2', 'C19MC', 'EWSR1']:
                # (6) 'RELA', 'NTRK1', 'NTRK2', 'NTRK3', 'NAB2', 'C19MC', 'EWSR1'
                # 基因变异类型写的融合：只关注融合，其他变异类型不纳入，检测结果栏呈现 包含该基因的融合对或者“未检出”
                pd_fusion_some = pd_fusion1[pd_fusion1['基因'] == each_gene]
                dict_nccn_result[each_gene] = '未检出' if pd_fusion_some.empty else ';'.join(list(pd_fusion_some['检测结果']))

            elif each_gene == 'BRAF':
                # (7) BRAF：检测类型为基因突变/基因融合：既关注点突变也同时关注融合，检测结果依据检出情况填写
                pd_mutation_BRAF = pd_mutation1[(pd_mutation1['基因'] == each_gene)]
                pd_fusion_BRAF = pd_fusion1[pd_fusion1['基因'] == each_gene]
                pd_merge_BRAF = pd.concat([pd_mutation_BRAF, pd_fusion_BRAF], axis=0)
                dict_nccn_result[each_gene] = '未检出' if pd_merge_BRAF.empty else ';'.join(list(pd_merge_BRAF['检测结果']))

            else:
                # (8) 除上述特殊检测之外的其他基因变异类型为基因突变的，该表格呈现突变有明确意义的位点，
                # 但指南等并未全部明确突变位点，故此栏呈现当前基因致病性结论为“体细胞较可靠突变,可能有害、体细胞突变,可能有害位点”

                list_type = dict_nccn_type[each_gene].split('/')

                if '点突变' in list_type or '基因突变' in list_type:
                    pd_mut_other = pd_mutation1[pd_mutation1['基因'] == each_gene]
                else:
                    pd_mut_other = pd.DataFrame(columns=list(pd_mutation1))

                if '基因融合' in list_type:
                    pd_fus_other = pd_fusion1[pd_fusion1['基因'] == each_gene]
                else:
                    pd_fus_other = pd.DataFrame(columns=list(pd_fusion1))

                if '基因扩增' in list_type:
                    pd_cnv_other = pd_cnv1[pd_cnv1['基因'] == each_gene]
                else:
                    pd_cnv_other = pd.DataFrame(columns=list(pd_cnv1))

                pd_merge_other = pd.concat([pd_mut_other, pd_fus_other, pd_cnv_other], axis=0)
                dict_nccn_result[each_gene] = '未检出' if pd_merge_other.empty else ';'.join(list(pd_merge_other['检测结果']))

        dict_nccn_glioma['检测结果'] = dict_nccn_result
        pd_nccn = pd.DataFrame(dict_nccn_glioma)
        pd_nccn.to_csv(f'{self.outdir}/NCCN.txt', sep='\t')
        pd_nccn = pd.read_csv(f'{self.outdir}/NCCN.txt', sep='\t')
        pd_nccn.columns = ['基因', '变异类型', '检测结果', '临床意义']
        pd_nccn.to_csv(f'{self.outdir}/NCCN.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/NCCN.txt')

    def unit_nccn_lung(self):
        """
        NCCN指南推荐基因检测 靶向药物+其它重要基因
        """

        pd_mut1, pd_mut2 = self.filter_mut(if_do_targeting_site_drug=False)
        pd_cnv1, pd_cnv2 = self.filter_cnv(if_do_targeting_site_drug=False)
        pd_fus1, pd_fus2 = self.filter_fus(if_do_targeting_site_drug=False)
        pd_solid_tumor_NCCN_gene = pd.read_excel(self.solid_tumor_NCCN_gene, self.cancer_flag)

        # 20210908新增NCCN数据库合并单元格读取问题解决
        if self.cancer_flag != '肺':
            pd_solid_tumor_NCCN_gene = pd_solid_tumor_NCCN_gene.fillna(method='ffill')

        if self.single_flag:
            header = ['基因', '突变率', 'mRNA变体名称(NM)', '染色体名称', '起始位置', '终止位置', '核酸改变', '氨基酸改变',
                      '突变类型', '结论', 'COSMIC数据库中有此突变的组织', 'SIFT危害性预测', 'Polyphen2_HDIV危害性预测', 'REVEL危害性预测',
                      'MutationTaster危害性预测', '外显子']
            pd_mut2 = pd.DataFrame(columns=header)
            pd_cnv2 = pd.DataFrame(columns=list(pd_cnv1))
            pd_fus2 = pd.DataFrame(columns=list(pd_fus1))

        for_list = [[pd_mut1, pd_cnv1, pd_fus1, self.tumor_id1, self.type1],
                    [pd_mut2, pd_cnv2, pd_fus2, self.tumor_id2, self.type2]]

        list_mut = []
        list_cnv = []
        list_fus = []
        for pd_mut_tmp, pd_cnv_tmp, pd_fus_tmp, tumor_id, sample_type in for_list:
            pd_mut = pd_mut_tmp.copy()
            pd_cnv = pd_cnv_tmp.copy()
            pd_fus = pd_fus_tmp.copy()

            header = ['基因', '氨基酸改变', '突变率']
            pd_mut = pd_mut.fillna(value='')
            pd_mut = pd_mut[header]
            pd_mut.columns = ['基因', '氨基酸改变', sample_type]
            pd_mut['氨基酸改变'] = [x if str(x) not in ['', ' '] else '(-)' for x in pd_mut['氨基酸改变']]
            pd_mut = pd_mut.drop_duplicates()
            list_mut.append(pd_mut)

            header = ['基因', '拷贝数']
            pd_cnv = pd_cnv[header]
            pd_cnv.columns = ['基因', sample_type]
            pd_cnv[sample_type] = [f'cn={x}' for x in list(pd_cnv[sample_type])]
            pd_cnv = pd_cnv.drop_duplicates()
            list_cnv.append(pd_cnv)

            header = ['基因', 'fusion_gene_name', 'value']
            pd_fus = pd_fus[header]
            pd_fus.columns = ['基因', 'fusion_gene_name', sample_type]
            pd_fus[sample_type] = pd_fus[sample_type].str.replace('/', '-')
            pd_fus = pd_fus.drop_duplicates()
            list_fus.append(pd_fus)

        pd_mut_merge = pd.merge(list_mut[0], list_mut[1], on=['基因', '氨基酸改变'], how='outer')
        list_mut_info = []
        for index, line in pd_mut_merge.iterrows():
            _list = []
            if str(line[self.type1]) not in ['', ' ', 'nan']:
                _list.append(f'{line[self.type1]}({self.type1})')
            if str(line[self.type2]) not in ['', ' ', 'nan']:
                _list.append(f'{line[self.type2]}({self.type2})')
            info = line['基因'] + ' ' + line['氨基酸改变'] + ' ' + '/'.join(_list)
            list_mut_info.append(info)
        pd_mut_merge['检测结果'] = list_mut_info
        pd_mut_merge = pd_mut_merge[['基因', '检测结果']]
        pd_mut_merge = pd_mut_merge.groupby(['基因'])['检测结果'].apply(lambda x: Tools.cluster_columns(x, sep=';'))
        pd_mut_merge = pd_mut_merge.reset_index()
        pd_mut_merge['检测结果'] = '突变:;' + pd_mut_merge['检测结果']

        pd_cnv_merge = pd.merge(list_cnv[0], list_cnv[1], on='基因', how='outer')
        list_cnv_info = []
        for index, line in pd_cnv_merge.iterrows():
            '基因 拷贝数'
            _list = []
            if str(line[self.type1]) not in ['', ' ', 'nan']:
                _list.append(f'{line[self.type1]}({self.type1})')
            if str(line[self.type2]) not in ['', ' ', 'nan']:
                _list.append(f'{line[self.type2]}({self.type2})')
            info = line['基因'] + ' ' + '/'.join(_list)
            list_cnv_info.append(info)
        pd_cnv_merge['检测结果'] = list_cnv_info
        pd_cnv_merge = pd_cnv_merge[['基因', '检测结果']]
        pd_cnv_merge = pd_cnv_merge.groupby(['基因'])['检测结果'].apply(lambda x: Tools.cluster_columns(x, sep=';'))
        pd_cnv_merge = pd_cnv_merge.reset_index()
        pd_cnv_merge['检测结果'] = '拷贝数变异:;' + pd_cnv_merge['检测结果']

        pd_fus_merge = pd.merge(list_fus[0], list_fus[1], on=['基因', 'fusion_gene_name'], how='outer')

        dict_rna_fus = {}
        if os.path.isfile(self.rna_fusion_file):
            header_tmp = ['fusion_gene_name', 'value']
            pd_rna_fus = self.filter_rna(if_do_targeting_site_drug=True)
            pd_rna_fus = pd_rna_fus[header_tmp]
            dict_rna_fus = {pd_rna_fus['fusion_gene_name'][x]: pd_rna_fus['value'][x] for x in range(len(pd_rna_fus))}

        list_fus_info = []
        for index, line in pd_fus_merge.iterrows():
            _list = []
            if str(line[self.type1]) not in ['', ' ', 'nan']:
                _list.append(f"{line[self.type1]}({self.type1})")
            if str(line[self.type2]) not in ['', ' ', 'nan']:
                _list.append(f"{line[self.type2]}({self.type2})")
            if str(line['fusion_gene_name']) not in dict_rna_fus.keys():
                info = f"{line['fusion_gene_name']} 融合;{'/'.join(_list)}"
            else:
                info = f"{line['fusion_gene_name']} 融合;DNA测序:{'/'.join(_list)};RNA测序:{dict_rna_fus[line['fusion_gene_name']]}"
            list_fus_info.append(info)

        pd_fus_merge['检测结果'] = list_fus_info
        pd_fus_merge = pd_fus_merge[['基因', '检测结果']].drop_duplicates()
        pd_fus_merge = pd_fus_merge.groupby(['基因'])['检测结果'].apply(lambda x: Tools.cluster_columns(x, sep=';'))
        pd_fus_merge = pd_fus_merge.reset_index()
        pd_fus_merge['检测结果'] = '重排/融合:;' + pd_fus_merge['检测结果']
        pd_fus_merge = pd_fus_merge[['基因', '检测结果']]

        pd_all = pd.concat([pd_mut_merge, pd_cnv_merge, pd_fus_merge], axis=0)
        pd_all = pd_all.groupby(['基因'])['检测结果'].apply(lambda x: Tools.cluster_columns(x, sep=';'))
        pd_all = pd_all.reset_index()
        pd_solid_tumor_NCCN_gene = pd_solid_tumor_NCCN_gene.drop('检测结果', axis=1)
        pd_all = pd.merge(pd_solid_tumor_NCCN_gene, pd_all, on='基因', how='left')
        pd_all = pd_all[['基因', '变异类型', '检测结果', '临床意义']]
        pd_all['检测结果'] = ['未检出' if str(x) in ['', ' ', 'nan'] else x for x in list(pd_all['检测结果'])]
        pd_all.to_csv(f'{self.outdir}/NCCN.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/NCCN.txt')

    def unit_2_0(self):  # add by rhd
        """内分泌模块：1去势治疗 2内分泌药物
        去势治疗只关注ADT	HSD3B1	1245C	AA()	敏感性	敏感	临床研究 这个变异来自胚系SNP，ctDNA好点儿
        野生参考基因型为C，GRCh37:chr1:120057246 rs1047303
        规则为 AA 敏感； AC ->敏感； CC->耐药；
        未检出该基因时默认为CC 耐药；其他基因型结果敏感性为‘未知’
        """
        # castration
        pd_cas = pd.read_excel(self.castration)  #, sep ='\t')# 单独跑一下chr1:119514623 rs1047303 的注释
        if pd_cas.empty:
            pd_cas['敏感性'] = '未知'
        pd_cas.to_csv(f'{self.outdir}/去势治疗疗效评估.txt', sep='\t', index=False)

        # endotherapy
        pd_mut1, pd_mut2 = self.filter_mut(if_do_targeting_site_drug=True)
        pd_cnv1, pd_cnv2 = self.filter_cnv(if_do_targeting_site_drug=True)
        if self.single_flag:
            header ='@@@'.join(list(pd_mut1)).replace(self.tumor_id1, self.tumor_id2).split('@@@')
            pd_mut2 = pd.DataFrame(columns=header)
            pd_cnv2 = pd.DataFrame(columns=list(pd_cnv1))
        dict_pd = {}
        for_list = [[pd_mut1, pd_cnv1, self.type1], [pd_mut2, pd_cnv2, self.type2]]
        for pd_mutation_temp, pd_cnv_temp, sample_type in for_list:
            pd_mut = pd_mutation_temp.copy()
            pd_cnv = pd_cnv_temp.copy()
            pd_mut = pd_mut[~pd_mut['内分泌药物'].isin([''])]
            pd_cnv = pd_cnv[~pd_cnv['内分泌药物'].isin([''])]

            header = ['基因', 'mRNA变体名称(NM)', '突变类型', '外显子', '核酸改变', '氨基酸改变', '突变率', '内分泌药物', '内分泌药物敏感性', '证据来源']
            pd_mut = pd_mut[header].drop_duplicates()
            pd_mut[sample_type] = pd_mut['突变率']
            pd_mut['基因'] = pd_mut['基因'].str.cat(pd_mut['mRNA变体名称(NM)'], sep=' ')
            pd_mut['核酸改变'] = [re.sub('\([a-z]*[A-Z]*[0-9]+\)', '', x) for x in list(pd_mut['核酸改变'])]
            pd_mut['变异'] = (pd_mut['外显子'].map(str) + ' ' +
                            pd_mut['核酸改变'].map(str) + ' ' +
                            pd_mut['氨基酸改变'].map(str) + ' ' +
                            pd_mut['突变类型'].map(str))

            header = ['基因', '拷贝数', '内分泌药物', '内分泌药物敏感性', '证据来源']
            pd_cnv = pd_cnv[header].drop_duplicates()
            pd_cnv[sample_type] = pd_cnv['拷贝数']
            pd_cnv['变异'] = '拷贝数扩增'

            header = ['基因', '变异', sample_type, '内分泌药物', '内分泌药物敏感性', '证据来源']
            pd_mut = pd_mut[header].drop_duplicates()
            pd_cnv = pd_cnv[header].drop_duplicates()
            dict_pd[sample_type] = {'pd_mut' : pd_mut, 'pd_cnv' : pd_cnv}

        pd_endo_mut = pd.merge(dict_pd[self.type1]['pd_mut'], dict_pd[self.type2]['pd_mut'], on=['基因', '变异', '内分泌药物', '内分泌药物敏感性', '证据来源'], how = 'outer')
        pd_endo_cnv = pd.merge(dict_pd[self.type1]['pd_cnv'], dict_pd[self.type2]['pd_cnv'], on=['基因', '变异', '内分泌药物', '内分泌药物敏感性', '证据来源'], how = 'outer')

        pd_endo = pd.concat([pd_endo_mut, pd_endo_cnv], axis=0, ignore_index=True)
        pd_endo = pd_endo[['基因', '变异', self.type1, self.type2, '内分泌药物', '内分泌药物敏感性', '证据来源' ]]
        pd_endo['证据来源'] = pd_endo['证据来源'].apply(lambda x: x.split(':')[1])
        pd_endo = pd_endo.sort_values(
            by=['基因', self.type1, self.type2, '变异'], ascending=[True, False, False, True])

        if self.single_flag:
            header = ['基因', '变异', self.type1, '内分泌药物', '内分泌药物敏感性', '证据来源']
        else:
            header = ['基因', '变异', self.type1, self.type2, '内分泌药物', '内分泌药物敏感性', '证据来源']
        pd_endo = pd_endo[header].fillna(value = '-')
        pd_endo.to_csv(f'{self.outdir}/内分泌治疗检测结果.txt', sep='\t', index=False)

    def unit_2_1(self):
        """ 本癌种 靶向用药提示 """
        pd_mut1, pd_mut2 = self.filter_mut(if_do_targeting_site_drug=True)
        pd_cnv1, pd_cnv2 = self.filter_cnv(if_do_targeting_site_drug=True)  # '基因', '变异（缺失/扩增）', '拷贝数'
        pd_fus1, pd_fus2 = self.filter_fus(if_do_targeting_site_drug=True)  # '基因', 'value', 'fusion_gene_name'

        if self.single_flag:
            header ='@@@'.join(list(pd_mut1)).replace(self.tumor_id1, self.tumor_id2).split('@@@')
            pd_mut2 = pd.DataFrame(columns=header)
            pd_cnv2 = pd.DataFrame(columns=list(pd_cnv1))
            pd_fus2 = pd.DataFrame(columns=list(pd_fus1))

        dict_pd = {}
        for_list = [[pd_mut1, pd_cnv1, pd_fus1, self.tumor_id1, self.type1],
                    [pd_mut2, pd_cnv2, pd_fus2, self.tumor_id2, self.type2]]

        for pd_mutation_temp, pd_cnv_temp, pd_fusion_temp, tumor_id, sample_type in for_list:
            pd_mut = pd_mutation_temp.copy()
            pd_cnv = pd_cnv_temp.copy()
            pd_fus = pd_fusion_temp.copy()

            pd_mut = pd_mut[~pd_mut['药物'].isin([''])]
            pd_mut['基因'] = pd_mut['基因'].str.cat(pd_mut['mRNA变体名称(NM)'], sep=' ')
            pd_cnv = pd_cnv[~pd_cnv['药物'].isin([''])]
            pd_fus = pd_fus[~pd_fus['药物'].isin([''])]
            # ['基因', 'fusion_gene_name', 'value', '疾病关键字', '药物', '敏感性预测', '证据来源']

            if os.path.isfile(self.rna_fusion_file) and str(sample_type).lower() != 'ctdna':
                # ['基因', 'fusion_gene_name', 'value', '疾病关键字', '药物', '敏感性预测', '证据来源']
                header = list(pd_fus)
                pd_rna_fusion = self.filter_rna(if_do_targeting_site_drug=True)
                pd_rna_fusion = pd_rna_fusion[header].drop_duplicates()
                pd_fus = pd.concat([pd_fus, pd_rna_fusion], axis=0)
                if not pd_fus.empty:
                    pd_fus = pd_fus.groupby(['基因', 'fusion_gene_name', '疾病关键字', '药物', '敏感性预测', '证据来源'])['value'].apply(lambda x: Tools.cluster_columns(x, sep=';'))
                    pd_fus = pd_fus.reset_index()
                list_value_format = []
                for i in list(pd_fus['value']):
                    if ';' in i:
                        value_info = str(i).split(';')
                        value_rna = [x for x in value_info if 'FFPM' in x][0]
                        value_dna = [x for x in value_info if '%' in x][0]
                        list_value_format.append(f"DNA:{value_dna};RNA:{value_rna}")
                    else:
                        list_value_format.append(i)
                pd_fus['value'] = list_value_format
                pd_fus = pd_fus[header]

            pd_mut_this = pd_mut[pd_mut['疾病关键字'] == self.cancer_flag]
            pd_cnv_this = pd_cnv[pd_cnv['疾病关键字'] == self.cancer_flag]
            pd_fus_this = pd_fus[pd_fus['疾病关键字'] == self.cancer_flag]

            pd_mut_no_this = pd_mut[~(pd_mut['疾病关键字'] == self.cancer_flag)]
            pd_cnv_no_this = pd_cnv[~(pd_cnv['疾病关键字'] == self.cancer_flag)]
            pd_fus_no_this = pd_fus[~(pd_fus['疾病关键字'] == self.cancer_flag)]

            header = ['基因', 'mRNA变体名称(NM)', '突变类型', '外显子', '核酸改变', '氨基酸改变', '突变率',
                      '药物', '敏感性预测', '证据来源']
            pd_mut_this = pd_mut_this[header]
            pd_mut_this = pd_mut_this.drop_duplicates()
            pd_mut_no_this = pd_mut_no_this[header]
            pd_mut_no_this = pd_mut_no_this.drop_duplicates()

            header = ['基因', '拷贝数', '药物', '敏感性预测', '证据来源']
            pd_cnv_this = pd_cnv_this[header]
            pd_cnv_this = pd_cnv_this.drop_duplicates()
            pd_cnv_no_this = pd_cnv_no_this[header]
            pd_cnv_no_this = pd_cnv_no_this.drop_duplicates()

            header = ['基因', 'fusion_gene_name', 'value', '药物', '敏感性预测', '证据来源']
            pd_fus_this = pd_fus_this[header]
            pd_fus_this = pd_fus_this.drop_duplicates()
            pd_fus_no_this = pd_fus_no_this[header]
            pd_fus_no_this = pd_fus_no_this.drop_duplicates()

            pd_mut_this['核酸改变'] = [re.sub('\([a-z]*[A-Z]*[0-9]+\)', '', x) for x in list(pd_mut_this['核酸改变'])]
            pd_mut_this['变异'] = (pd_mut_this['外显子'].map(str) + ' ' +
                                 pd_mut_this['核酸改变'].map(str) + ' ' +
                                 pd_mut_this['氨基酸改变'].map(str) + ' ' +
                                 pd_mut_this ['突变类型'].map(str)
                                 )
            pd_mut_this[sample_type] = pd_mut_this['突变率']
            pd_cnv_this['变异'] = [f"{x} 拷贝数扩增" for x in list(pd_cnv_this['基因'])]
            pd_cnv_this[sample_type] = pd_cnv_this['拷贝数']
            pd_fus_this['变异'] = [f"{x} 融合" for x in list(pd_fus_this['fusion_gene_name'])]
            pd_fus_this[sample_type] = pd_fus_this['value']

            pd_mut_no_this['核酸改变'] = [re.sub('\([a-z]*[A-Z]*[0-9]+\)', '', x) for x in list(pd_mut_no_this['核酸改变'])]
            pd_mut_no_this['变异'] = pd_mut_no_this['外显子'].map(str) + ' ' + \
                                   pd_mut_no_this['核酸改变'].map(str) + ' ' + \
                                   pd_mut_no_this['氨基酸改变'].map(str) + ' ' + \
                                   pd_mut_no_this['突变类型'].map(str)
            pd_mut_no_this[sample_type] = pd_mut_no_this['突变率']
            pd_cnv_no_this['变异'] = [f"{x} 拷贝数扩增" for x in list(pd_cnv_no_this['基因'])]
            pd_cnv_no_this[sample_type] = pd_cnv_no_this['拷贝数']
            pd_fus_no_this['变异'] = [f"{x} 融合" for x in list(pd_fus_no_this['fusion_gene_name'])]
            pd_fus_no_this[sample_type] = pd_fus_no_this['value']

            header = ['基因', '变异', sample_type, '药物', '敏感性预测', '证据来源']
            pd_mut_this = pd_mut_this[header].drop_duplicates()
            pd_cnv_this = pd_cnv_this[header].drop_duplicates()
            pd_fus_this = pd_fus_this[header].drop_duplicates()
            pd_mut_no_this = pd_mut_no_this[header].drop_duplicates()
            pd_cnv_no_this = pd_cnv_no_this[header].drop_duplicates()
            pd_fus_no_this = pd_fus_no_this[header].drop_duplicates()

            dict_pd[sample_type] = {'pd_mut_this': pd_mut_this,
                                    'pd_cnv_this': pd_cnv_this,
                                    'pd_fus_this': pd_fus_this,
                                    'pd_mut_no_this': pd_mut_no_this,
                                    'pd_cnv_no_this': pd_cnv_no_this,
                                    'pd_fus_no_this': pd_fus_no_this}

        merge_mutation_this = pd.merge(dict_pd[self.type1]['pd_mut_this'],
                                       dict_pd[self.type2]['pd_mut_this'],
                                       on=['基因', '变异', '药物', '敏感性预测', '证据来源'],
                                       how='outer')
        pd_merge_cnv_this = pd.merge(dict_pd[self.type1]['pd_cnv_this'],
                                     dict_pd[self.type2]['pd_cnv_this'],
                                     on=['基因', '变异', '药物', '敏感性预测', '证据来源'],
                                     how='outer')
        merge_fusion_this = pd.merge(dict_pd[self.type1]['pd_fus_this'],
                                     dict_pd[self.type2]['pd_fus_this'],
                                     on=['基因', '变异', '药物', '敏感性预测', '证据来源'],
                                     how='outer')
        merge_mutation_no_this = pd.merge(dict_pd[self.type1]['pd_mut_no_this'],
                                          dict_pd[self.type2]['pd_mut_no_this'],
                                          on=['基因', '变异', '药物', '敏感性预测', '证据来源'],
                                          how='outer')
        merge_cnv_no_this = pd.merge(dict_pd[self.type1]['pd_cnv_no_this'],
                                     dict_pd[self.type2]['pd_cnv_no_this'],
                                     on=['基因', '变异', '药物', '敏感性预测', '证据来源'],
                                     how='outer')
        merge_fusion_no_this = pd.merge(dict_pd[self.type1]['pd_fus_no_this'],
                                        dict_pd[self.type2]['pd_fus_no_this'],
                                        on=['基因', '变异', '药物', '敏感性预测', '证据来源'],
                                        how='outer')

        out_unit_2_1_this = pd.concat([merge_mutation_this, pd_merge_cnv_this, merge_fusion_this], axis=0, ignore_index=True)
        out_unit_2_1_no_this = pd.concat([merge_mutation_no_this, merge_cnv_no_this, merge_fusion_no_this], axis=0, ignore_index=True)

        out_unit_2_1_this = out_unit_2_1_this[['基因', '变异', self.type1, self.type2, '药物', '敏感性预测', '证据来源'	]]
        out_unit_2_1_this['证据来源'] = out_unit_2_1_this['证据来源'].apply(lambda x: x.split(':')[1])
        out_unit_2_1_this = out_unit_2_1_this.sort_values(
            by=['基因', self.type1, self.type2, '变异'], ascending=[True, False, False, True])

        out_unit_2_1_this_flag = list(out_unit_2_1_this['变异'].str.cat(out_unit_2_1_this['药物'], sep=' '))
        out_unit_2_1_no_this['flag'] = out_unit_2_1_no_this['变异'].str.cat(out_unit_2_1_no_this['药物'], sep=' ')
        out_unit_2_1_no_this = out_unit_2_1_no_this[~out_unit_2_1_no_this['flag'].isin(out_unit_2_1_this_flag)]
        out_unit_2_1_no_this = out_unit_2_1_no_this[['基因', '变异', self.type1, self.type2, '药物', '敏感性预测', '证据来源']]
        out_unit_2_1_no_this = out_unit_2_1_no_this.sort_values(
            by=['基因', self.type1, self.type2, '变异'], ascending=[True, False, False, True])

        if self.single_flag:
            out_unit_2_1_this = out_unit_2_1_this[['基因', '变异', self.type1, '药物',  '敏感性预测', '证据来源']]
            out_unit_2_1_no_this = out_unit_2_1_no_this[['基因', '变异', self.type1, '药物', '敏感性预测', '证据来源']]

        # 20210908肠癌新增(增加肠癌中KRAS/NRAS/BRAF三基因全阴性处理方案)
        def specical_gene(gene):
            if gene.split(' ')[0] in ['KRAS', 'NRAS', 'BRAF']:
                return True
        special_gene_frame = out_unit_2_1_this[out_unit_2_1_this['基因'].apply(specical_gene) == True]
        if self.cancer_flag == '肠' and special_gene_frame.empty:
            special_info = list()
            for gene in ['KRAS', 'BRAF', 'NRAS']:
                if self.single_flag:
                    special_info.append([gene, gene + ' 野生型', '-', '西妥昔单抗', '敏感', 'FDA/NMPA/NCCN'])
                    special_info.append([gene, gene + ' 野生型', '-', '帕尼单抗', '敏感', 'FDA/NMPA/NCCN'])
                else:
                    special_info.append([gene, gene + ' 野生型', '-', '-', '西妥昔单抗', '敏感', 'FDA/NMPA/NCCN'])
                    special_info.append([gene, gene + ' 野生型', '-', '-', '帕尼单抗', '敏感', 'FDA/NMPA/NCCN'])
            out_unit_2_1_this = out_unit_2_1_this.append(pd.DataFrame(special_info, columns=list(out_unit_2_1_this)),
                                                         ignore_index=True)
        out_unit_2_1_this = out_unit_2_1_this.fillna(value='-')
        out_unit_2_1_this.to_csv(f'{self.outdir}/本癌种相关靶向药物.txt', sep='\t', index=False)
        out_unit_2_1_no_this = out_unit_2_1_no_this.fillna(value='-')
        out_unit_2_1_no_this.to_csv(f'{self.outdir}/其他癌种相关靶向药物.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/本癌种相关靶向药物.txt')
        self.result_file_list.append(f'{self.outdir}/其他癌种相关靶向药物.txt')

    def unit_2_2(self):
        pd_mut1, pd_mut2 = self.filter_mut(if_do_targeting_site_drug=True)
        pd_cnv1, pd_cnv2 = self.filter_cnv(if_do_targeting_site_drug=True)
        pd_fus1, pd_fus2 = self.filter_fus(if_do_targeting_site_drug=True)

        if self.single_flag:
            pd_mut2 = pd.DataFrame(columns=list(pd_mut1))
            pd_cnv2 = pd.DataFrame(columns=list(pd_cnv1))
            pd_fus2 = pd.DataFrame(columns=list(pd_fus1))

        list_mut = []
        list_cnv = []
        list_fus = []
        for_list = [[pd_mut1, pd_cnv1, pd_fus1, self.tumor_id1, self.type1],
                    [pd_mut2, pd_cnv2, pd_fus2, self.tumor_id2, self.type2]]

        for pd_mut, pd_cnv, pd_fus, tumor_id, sample_type in for_list:
            """
            2.2.1 重要基因变异检测结果 ---基因点突变、插入和缺失
                靶向药物 + 其它重要基因 + 不过癌种
            """
            header = ['基因', 'mRNA变体名称(NM)', '外显子', '核酸改变', '氨基酸改变', '突变率', '突变类型']
            pd_mut = pd_mut[header].drop_duplicates()
            pd_mut = pd_mut.apply(lambda x: Tools.rm_exon(x), axis=1)
            pd_mut.columns = ['基因', '转录本ID', '外显子/内含子', '核苷酸改变', '氨基酸改变', sample_type, '变异来源']
            list_mut.append(pd_mut)

            """ 
            2.2.2 重要基因变异检测结果 --- 基因重排/融合
                靶向药物 + 其它重要基因 + 不过癌种
            排序  基因 —— 分度
            """
            pd_fus = pd_fus[['基因', 'fusion_gene_name', 'value']].drop_duplicates()
            pd_fus.columns = ['基因', '重排/融合方式', sample_type]
            list_fus.append(pd_fus)

            """ 
            2.2.3 重要基因变异检测结果 --- 基因拷贝数变异
            靶向药物 + 其它重要基因 + 不过癌种
            """
            header = ['基因', '变异（缺失/扩增）', '拷贝数']
            pd_cnv = pd_cnv[header].drop_duplicates()
            pd_cnv.columns = ['基因', '缺失/扩增', sample_type]
            list_cnv.append(pd_cnv)

        pd_mut_merge = pd.merge(list_mut[0], list_mut[1], on=['基因', '转录本ID', '外显子/内含子', '核苷酸改变', '氨基酸改变', '变异来源'], how='outer')
        pd_mut_merge = pd_mut_merge[['基因', '转录本ID', '外显子/内含子', '核苷酸改变', '氨基酸改变', self.type1, self.type2, '变异来源']]
        pd_mut_merge = pd_mut_merge.sort_values(by=['基因', self.type1, self.type2], ascending=[True, False, False])

        pd_fus_merge = pd.merge(list_fus[0], list_fus[1], on=['基因', '重排/融合方式'], how='outer')
        pd_fus_merge = pd_fus_merge.sort_values(by=['基因', '重排/融合方式'], ascending=[True, True])
        pd_fus_merge = pd_fus_merge.fillna(value='-')

        dict_rna_fus = {}

        if os.path.isfile(self.rna_fusion_file):
            header_tmp = ['fusion_gene_name', 'value']
            pd_rna_fus = self.filter_rna(if_do_targeting_site_drug=True)
            pd_rna_fus = pd_rna_fus[header_tmp]
            dict_rna_fus = {pd_rna_fus['fusion_gene_name'][x]: pd_rna_fus['value'][x] for x in range(len(pd_rna_fus))}

        info_list = []
        for index, line in pd_fus_merge.iterrows():
            fusion_gene = str(line['基因'])
            fusion_pair = str(line['重排/融合方式'])
            value_type1 = str(line[self.type1])
            value_type2 = str(line[self.type2])
            dna_info = [fusion_gene, fusion_pair, 'DNA', value_type1, value_type2]
            rna_info = [fusion_gene, fusion_pair, 'RNA', dict_rna_fus.get(fusion_pair, '-'), '-']
            info_list.append(dna_info)
            if os.path.isfile(self.rna_fusion_file):
                info_list.append(rna_info)

        pd_fus_merge = pd.DataFrame(info_list, columns=['基因', '重排/融合方式', '样本', self.type1, self.type2])
        pd_cnv_merge = pd.merge(list_cnv[0], list_cnv[1], on=['基因', '缺失/扩增'], how='outer')
        pd_cnv_merge = pd_cnv_merge.sort_values(by=['基因', self.type1, self.type2], ascending=[True, False, False])

        if self.single_flag:
            pd_mut_merge = pd_mut_merge[['基因', '转录本ID', '外显子/内含子', '核苷酸改变', '氨基酸改变', self.type1, '变异来源']]
            pd_fus_merge = pd_fus_merge[['基因', '重排/融合方式', '样本', self.type1]]
            pd_cnv_merge = pd_cnv_merge[['基因', '缺失/扩增', self.type1]]

        pd_mut_merge = pd_mut_merge.fillna(value='-')
        pd_mut_merge.to_csv(f'{self.outdir}/重要基因变异检测结果_mutation.txt', sep='\t', index=False)
        pd_fus_merge = pd_fus_merge.fillna(value='-')
        pd_fus_merge.to_csv(f'{self.outdir}/重要基因变异检测结果_fusion.txt', sep='\t', index=False)
        pd_cnv_merge = pd_cnv_merge.fillna(value='-')
        pd_cnv_merge.to_csv(f'{self.outdir}/重要基因变异检测结果_cnv.txt', sep='\t', index=False)

        self.result_file_list.append(f'{self.outdir}/重要基因变异检测结果_mutation.txt')
        self.result_file_list.append(f'{self.outdir}/重要基因变异检测结果_fusion.txt')
        self.result_file_list.append(f'{self.outdir}/重要基因变异检测结果_cnv.txt')


        """
        self.tumor_id1 = self.argv.tumor_id1
        self.tumor_id2 = self.argv.tumor_id2
        self.tumor_wes_id1 = self.argv.tumor_wes_id1
        self.tumor_wes_id2 = self.argv.tumor_wes_id2
        """
        for sample_id in [x for x in [self.tumor_id1, self.tumor_id2, self.tumor_wes_id1, self.tumor_wes_id2] if x]:
            os.makedirs(f'{self.outdir}/IGV_html/{sample_id}', exist_ok=True)
            if self.single_flag:
                args = f'--merge_new_file {self.mutation_file1}'
            else:
                args = f'--merge_new_file {self.mutation_file1},{self.mutation_file2}'
            scirpt = fr"""
python /mnt/share02/zhangxc/pipeline/IGV_html_pipeline/IGV_html.py \
	--mutation_file {self.outdir}/重要基因变异检测结果_mutation.txt \
	{args} \
	--bam_file {self.outdir}/../../bwa_alignment/{sample_id}/{sample_id}.final.bam \
	--sample_id {sample_id} \
	--out_dir {self.outdir}/IGV_html/{sample_id}
"""
        # try:
        # 	os.system(scirpt)
        # except:
        # 	print(scirpt, '\nIGV_html.py 执行失败')

    def unit_3(self):
        """ 1.2 中已输出 """
        pass

    def unit_4(self):
        """ 药物基因组 """
        pd_normal_bam_readcount = pd.read_excel(self.normal_bam_readcount)
        pd_normal_bam_readcount.columns = ['chr', 'pos', 'end', 'rs号', 'Deth', '频率', '结果', '次频率', '次结果']

        pd_chemotherapy_annotation = pd.read_excel(self.solid_tumor_chemotherapy_annotation)

        # 靶向药物list
        solid_tumor_chemotherapy_drug_list = pd.read_excel(self.solid_tumor_chemotherapy_drug_list)
        solid_tumor_chemotherapy_drug_list = solid_tumor_chemotherapy_drug_list.fillna(value='')
        drug_list = []
        for each_drug in solid_tumor_chemotherapy_drug_list[self.cancer_flag]:
            if each_drug not in ['', ' ']:
                drug_list.append(each_drug)

        list_gene_type = []
        for x in [str(x) for x in range(1, 7)]:
            gene_type = pd.merge(pd_normal_bam_readcount, pd_chemotherapy_annotation, left_on=['rs号', '结果'],
                                 right_on=['检测位点', '基因型' + x], how='inner')
            gene_type['基因型all'] = gene_type['基因型' + x]
            gene_type['基因型说明all'] = gene_type['基因型%s说明' % x]
            list_gene_type.append(gene_type)
        pd_all_gene_type = pd.concat(list_gene_type, axis=0, ignore_index=True)

        def judge_drug(x):
            _my_drug_list = set(str(x['药物']).split(';'))
            _need_drug_list = set(drug_list)
            _drug_intersection = _my_drug_list & _need_drug_list
            if len(_drug_intersection) == 0:
                return '不在靶向药list'
            else:
                return '在靶向药list'

        pd_all_gene_type['是否在靶向药list'] = ''
        pd_all_gene_type['是否在靶向药list'] = pd_all_gene_type.apply(lambda x: judge_drug(x), axis=1)
        pd_all_gene_type = pd_all_gene_type[pd_all_gene_type['是否在靶向药list'] == '在靶向药list']

        pd_all_gene_type = pd_all_gene_type[['药物', '基因', '检测位点', '分析内容', '证据等级', '基因型all', '基因型说明all', '基因型说明(评分)']]
        pd_all_gene_type.columns = ['药物', '基因', '检测位点', '分析内容', '证据等级', '基因型', '基因型说明', '基因型说明(评分)']
        # 将药物列拆开
        if not pd_all_gene_type.empty:
            pd_all_gene_type = pd_all_gene_type.drop('药物', axis=1).join(
                pd_all_gene_type['药物'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('药物'))

        # 排序
        pd_all_gene_type['证据等级2'] = ''
        pd_all_gene_type['证据等级2'] = [str(x).replace('-1A', '0').replace('1B', '1').replace('2A', '21').replace('2B', '22').replace('3', '30').replace(
            '4', '40') for x in list(pd_all_gene_type['证据等级'])]
        # pd_all_gene_type['证据等级2'] = pd_all_gene_type['证据等级'].apply(lambda x: re.sub("[A-Za-z]", "", str(x)))
        pd_all_gene_type = pd_all_gene_type.sort_values(by=['药物', '证据等级2', '检测位点', '基因'], ascending=[True, True, True, True])
        pd_all_gene_type = pd_all_gene_type[['基因', '检测位点', '分析内容', '证据等级', '基因型', '基因型说明', '药物', '基因型说明(评分)']]

        list_level = []
        for index, line in pd_all_gene_type.iterrows():
            gene_type = line['基因型']
            dict_gene_level = dict([x.split('::') for x in line['基因型说明(评分)'].split(';')])
            gene_level = dict_gene_level.get(gene_type, '-')
            list_level.append(gene_level)
        pd_all_gene_type['基因型说明(评分)'] = list_level
        # 小panel基因列表过滤

        def filter_drug_list(x):
            little_panel = ['NCPP1', 'NCPP2', 'NCPP3', 'NCPU1', 'NCPU2', 'NCPU3']
            if self.project_id in little_panel:
                sheet_dict = {'前列腺' : '前列腺13化疗药基因', '尿路上皮' : '尿路上皮17化疗药基因'}
                pd_panel_gene = pd.read_excel(self.solid_gene_list, sheet_name=sheet_dict[self.cancer_flag], header = None)
                filter_gene_list = [x[0] for x in pd_panel_gene.values.tolist()]
                x = x[x['基因'].isin(filter_gene_list)]
            return x
        pd_all_gene_type = filter_drug_list(pd_all_gene_type)
        pd_all_gene_type = pd_all_gene_type[['基因', '检测位点', '分析内容', '证据等级', '基因型', '基因型说明', '药物', '基因型说明(评分)']]
        pd_all_gene_type.to_csv(f'{self.outdir}/药物基因组学相关检测结果.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/药物基因组学相关检测结果.txt')

        return pd_all_gene_type

    def unit_4_conclusion(self, pd_all_gene_type):
        """化疗药评分得到最终结论"""
        def score_level(chemo_frame, flag):
            """统计降低、中等、升高不同程度得分"""
            def score(frame, chemo_score):
                trans_lib = {'升高': '较高', '中等': '中等', '降低': '较低'}
                score_lib = {'1A': 100, '1B': 100, '2A': 6, '2B': 5, '3': 2, '4': 1}
                score = score_lib[str(frame['证据等级'])]
                state = ['升高', '降低', '中等']
                if frame['基因型说明(评分)'] in state:
                    if frame['药物'] not in chemo_score:
                        chemo_score[frame['药物']] = dict()
                        chemo_score[frame['药物']]['较低'] = 0
                        chemo_score[frame['药物']]['中等'] = 0
                        chemo_score[frame['药物']]['较高'] = 0
                    chemo_score[frame['药物']][trans_lib[frame['基因型说明(评分)']]] += score
                else:
                    warning = '{}\t{}\t{}\t{}\t{}'.format(frame['药物'], frame['基因'], frame['基因型'], frame['分析内容'],
                                                          frame['基因型说明(评分)'])
                    print('Please check 分析内容：\n' + warning)

            chemo_score = dict()
            filter_frame = chemo_frame[chemo_frame['分析内容'] == flag]
            filter_frame.apply(lambda filter_frame: score(filter_frame, chemo_score), axis=1)
            return chemo_score

        def conclusion_line(drug, chemo_lib):
            """
            得到每个药物输出行内容
            """
            drug_line = '尚未知\t尚未知\t尚未知\t尚未知'
            result = '尚未知'
            if drug in chemo_lib:
                sort_lib = sorted(chemo_lib[drug].items(), key=lambda x: x[1], reverse=True)
                if sort_lib[0][1] > sort_lib[1][1]:
                    result = sort_lib[0][0]
                else:
                    result = '-'
                drug_line = '{}\t{}\t{}\t{}'.format(chemo_lib[drug]['较低'], chemo_lib[drug]['中等'], chemo_lib[drug]['较高'],
                                                    result)
            return drug_line, result

        effect_chemo = score_level(pd_all_gene_type, '疗效')
        toxic_chemo = score_level(pd_all_gene_type, '毒性')
        drugs = list()
        [drugs.append(i) for i in effect_chemo if i not in drugs]
        [drugs.append(i) for i in toxic_chemo if i not in drugs]
        final_chemo = list()
        chemo_report = list()
        for drug in drugs:
            effect_line, effect_result = conclusion_line(drug, effect_chemo)
            toxic_line, toxic_result = conclusion_line(drug, toxic_chemo)
            final_chemo.append('{}\t{}\t{}\n'.format(drug, effect_line, toxic_line))
            chemo_report.append('{}\t{}\t{}\n'.format(drug, effect_result, toxic_result))

        with open(f'{self.outdir}/药物基因组结论.txt', 'w') as f:
            f.write('药物\t敏感降低\t敏感中等\t敏感升高\t疗效结论\t毒性降低\t毒性中等\t毒性升高\t毒性结论\n')
            f.write(''.join(final_chemo))

        with open(f'{self.outdir}/药物基因组敏感性及毒副作用.txt', 'w') as g:
            g.write('药物\t敏感性\t毒副作用风险\n')
            g.write(''.join(chemo_report))

    def unit_5(self):
        """ 遗传性肿瘤相关基因突变检测结果 """
        """ 基因	转录本ID	外显子/内含子	核苷酸改变	氨基酸改变	纯合/杂合	致病风险	遗传方式 """
        pd_tumor_genetic_predisposed_gene_list = pd.read_csv(self.tumor_genetic_predisposed_gene_list, sep='\t')
        pd_mut = pd.read_csv(self.mutation_file1, sep='\t', low_memory=False, keep_default_na=False)
        pd_mut['外显子'] = pd_mut['核酸改变'].apply(lambda x: Tools.extrac_exon(x))
        # pd_mut = pd_mut.fillna(value='') 不可以进行替换NA  因为 pathogenicity 这块需要提取NA的信息
        pd_mut = pd_mut.apply(lambda x: Tools.rm_exon(x), axis=1)

        # (1) 纯合、杂合、半合子 判断
        pd_mut['突变状态'] = ''
        pd_mut["突变状态"] = pd_mut.apply(lambda x: Tools.mut_status(x, self.gender, self.normal_id), axis=1)

        # (2) 致病性、可能致病性、NA 判断  去掉其他无意义注释
        pd_mut['pathogenicity'] = ''
        pd_mut["pathogenicity"] = pd_mut.apply(lambda x: Tools.get_pathogenicity(x), axis=1)

        '基因	转录本ID	外显子/内含子	核苷酸改变	氨基酸改变	纯合/杂合	致病风险	遗传方式'

        # 检不出结果可以解开下面两行这样测试
        # pd_mut['突变级别'] = '1'
        # pd_mut['pathogenicity'] = '致病'
        # pd_mut.to_csv('遗传性肿瘤基因检测结果temp.txt', sep='\t', index=False)

        # (3) 进行'遗传性突变'的过滤
        pd_mut['结论'] = pd_mut.apply(lambda x: Tools.get_conclusion(x, self.product, panel_id=self.tumor_id1, wes_id=self.tumor_wes_id1), axis=1)
        pd_mut = pd_mut[(pd_mut['结论'] == '遗传性突变') &
                        pd_mut['突变级别'].isin(['1', '2']) &
                        ~pd_mut['突变类型'].isin(['非编码区突变', '同义突变'])&
                        (pd_mut['pathogenicity'].isin(['致病', '可能致病', '-']))
                        ]
        pd_mut = pd_mut[['基因', 'mRNA变体名称(NM)', '外显子', '核酸改变', '氨基酸改变', '突变状态', 'pathogenicity', '遗传方式']]

        pd_merge = pd.merge(pd_mut, pd_tumor_genetic_predisposed_gene_list, on='基因', how='inner')

        # (5) 提取此模块需要的列
        pd_merge = pd_merge[['基因', 'mRNA变体名称(NM)', '外显子', '核酸改变', '氨基酸改变',
                             '突变状态', 'pathogenicity', '遗传方式', '相关癌种', '基因描述', '疾病表型词']]

        def filter_germline_gene(pd_merge):
            little_panel = ['NCPK1', 'NCPK2', 'NCPK3']
            sheet_dict = {'肾': '肾32遗传性基因'}
            if self.project_id in little_panel:
                pd_panel_gene = pd.read_excel(self.solid_gene_list, sheet_name=sheet_dict[self.cancer], header=None)  # 库表 癌症为小类
                filter_gene_list = [x[0] for x in pd_panel_gene.values.tolist()]
                pd_merge = pd_merge[pd_merge['基因'].isin(filter_gene_list)]
            return pd_merge
        pd_merge = filter_germline_gene(pd_merge)

        # (6) 进行列的合并
        # pd_merge.colums = ['基因', '转录本ID', '外显子/内含子', '核苷酸改变', '氨基酸改变', '核酸改变', '氨基酸改变',
        # 				   '纯合/杂合', '致病风险', '遗传方式', '相关癌种', '基因描述', '疾病表型词']

        # (7) 遗传方式 行进行折叠 暂未发现有相同的突变输出多行，认为不需要
        # pd_mut = pd_mut.groupby(['基因', '核酸改变', '氨基酸改变'])['遗传方式'].apply(lambda x: Tools.cluster_columns(x, sep=';'))
        # pd_mut = pd_mut.reset_index()

        panel39_code = ['NCPC1', 'NCPC2', 'NCPC3']
        if self.project_id in panel39_code:
            panel39_gene_frame = pd.read_excel(panel39_germline_genes, sheet_name = 'panel39')
            pd_merge = pd.merge(pd_merge, panel39_gene_frame)
        pd_merge.to_csv(f'{self.outdir}/遗传性肿瘤相关基因突变检测结果.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/遗传性肿瘤相关基因突变检测结果.txt')

    def unit_b(self):
        """ 本次基因检测结果总览（体细胞变异） """

        pd_mut1, pd_mut2 = self.filter_mut(if_do_targeting_site_drug=False)
        pd_mut1_db, pd_mut2_db = self.filter_mut(if_do_targeting_site_drug=True)
        pd_cnv1, pd_cnv2 = self.filter_cnv(if_do_targeting_site_drug=True)
        pd_fus1, pd_fus2 = self.filter_fus(if_do_targeting_site_drug=True)
        head_num = 50
        if self.single_flag:
            pd_mut2 = pd.DataFrame(columns=list(pd_mut1))
            pd_mut2_db = pd.DataFrame(columns=list(pd_mut1))
            pd_cnv2 = pd.DataFrame(columns=list(pd_cnv1))
            pd_fus2 = pd.DataFrame(columns=list(pd_fus1))
            head_num = 100

        for_list = [[pd_mut1, pd_cnv1, pd_fus1, self.tumor_id1, self.type1],
                    [pd_mut2, pd_cnv2, pd_fus2, self.tumor_id2, self.type2]]

        list_mut = []
        list_cnv = []
        list_fus = []

        for_list_db = [[pd_mut1_db, self.tumor_id1, self.type1], [pd_mut2_db, self.tumor_id2, self.type2]]
        list_mut_db = []
        for pd_mut_tmp, tumor_id, sample_type in for_list_db:
            pd_mut = pd_mut_tmp.copy()
            pd_mut[sample_type] = pd_mut['突变率']
            # pd_mut[sample_type] = pd_mut['突变率'].str.strip('%').astype(float) / 100
            def unit_b_convert(x):
                if x['氨基酸改变'] in ['', ' ', 'NA']:
                    return '(-)'
                else:
                    return x['氨基酸改变']

            pd_mut['氨基酸改变'] = pd_mut.apply(lambda x: unit_b_convert(x), axis=1)
            header = ['变异基因', '外显子/内含子', '变异位点', sample_type, '变异类型', '数据库收录(COSMIC)', 'SIFT',
                      'Polyphen2', 'REVEL', 'MutationTaster']
            if not pd_mut.empty:
                pd_mut['变异基因'] = pd_mut['基因']
                pd_mut['外显子/内含子'] = pd_mut['外显子']
                pd_mut['变异位点'] = pd_mut['mRNA变体名称(NM)'] + ':;' + pd_mut['核酸改变'] + ':;' + pd_mut['氨基酸改变']
                # pd_mut[sample_type] = pd_mut[sample_type].map(lambda x: format(x, '.1%'))
                pd_mut['变异类型'] = pd_mut['突变类型']
                pd_mut['数据库收录(COSMIC)'] = pd_mut.apply(lambda x: Tools.convet_cosmic(x), axis=1)
                pd_mut['SIFT'] = pd_mut['SIFT危害性预测'].apply(lambda x: x.split('(')[0])
                pd_mut['Polyphen2'] = pd_mut['Polyphen2_HDIV危害性预测'].apply(lambda x: x.split('(')[0])
                pd_mut['REVEL'] = pd_mut['REVEL危害性预测'].apply(lambda x: x.split('(')[0])
                pd_mut['MutationTaster'] = pd_mut['MutationTaster危害性预测'].apply(lambda x: x.split('(')[0])
                pd_mut = pd_mut[header]
                pd_mut = pd_mut.drop_duplicates()
            else:
                pd_mut = pd.DataFrame(columns=header)
            list_mut_db.append(pd_mut)

        for pd_mut_tmp, pd_cnv_tmp, pd_fus_tmp, tumor_id, sample_type in for_list:
            pd_mut = pd_mut_tmp.copy()
            pd_mut = pd_mut.head(head_num)
            pd_cnv = pd_cnv_tmp.copy()
            pd_fus = pd_fus_tmp.copy()

            pd_mut[sample_type] = pd_mut['突变率'].str.strip('%').astype(float) / 100

            def unit_b_convert(x):
                if x['氨基酸改变'] in ['', ' ', 'NA']:
                    return '(-)'
                else:
                    return x['氨基酸改变']

            pd_mut['氨基酸改变'] = pd_mut.apply(lambda x: unit_b_convert(x), axis=1)

            header = ['变异基因', '外显子/内含子', '变异位点', sample_type, '变异类型', '数据库收录(COSMIC)', 'SIFT',
                      'Polyphen2', 'REVEL', 'MutationTaster']

            if not pd_mut.empty:
                pd_mut['变异基因'] = pd_mut['基因']
                pd_mut['外显子/内含子'] = pd_mut['外显子']
                pd_mut['变异位点'] = pd_mut['mRNA变体名称(NM)'] + ':;' + pd_mut['核酸改变'] + ':;' + pd_mut['氨基酸改变']
                # pd_mut[sample_type] = pd_mut[sample_type].map(lambda x: format(x, '.1%'))
                pd_mut['变异类型'] = pd_mut['突变类型']
                pd_mut['数据库收录(COSMIC)'] = pd_mut.apply(lambda x: Tools.convet_cosmic(x), axis=1)
                pd_mut['SIFT'] = pd_mut['SIFT危害性预测'].apply(lambda x: x.split('(')[0])
                pd_mut['Polyphen2'] = pd_mut['Polyphen2_HDIV危害性预测'].apply(lambda x: x.split('(')[0])
                pd_mut['REVEL'] = pd_mut['REVEL危害性预测'].apply(lambda x: x.split('(')[0])
                pd_mut['MutationTaster'] = pd_mut['MutationTaster危害性预测'].apply(lambda x: x.split('(')[0])
                pd_mut = pd_mut[header]
                pd_mut = pd_mut.drop_duplicates()
            else:
                pd_mut = pd.DataFrame(columns=header)
            list_mut.append(pd_mut)

            if not pd_cnv.empty:
                pd_cnv['变异基因'] = pd_cnv['基因']
                pd_cnv['外显子/内含子'] = '-'
                pd_cnv['变异位点'] = '-'
                pd_cnv[sample_type] = pd_cnv['拷贝数']
                pd_cnv['变异类型'] = '扩增'
                pd_cnv['数据库收录(COSMIC)'] = '-'
                pd_cnv['SIFT'] = '-'
                pd_cnv['Polyphen2'] = '-'
                pd_cnv['REVEL'] = '-'
                pd_cnv['MutationTaster'] = '-'
                pd_cnv = pd_cnv[header]
                pd_cnv = pd_cnv.drop_duplicates()
            else:
                pd_cnv = pd.DataFrame(columns=header)
            list_cnv.append(pd_cnv)

            if not pd_fus.empty:
                if os.path.isfile(self.rna_fusion_file) and str(sample_type).lower() != 'ctdna':
                    header_tmp = list(pd_fus)
                    pd_rna_fusion = self.filter_rna(if_do_targeting_site_drug=True)
                    pd_rna_fusion = pd_rna_fusion[header_tmp].drop_duplicates()
                    pd_fus = pd.concat([pd_fus, pd_rna_fusion], axis=0)
                    pd_fus = pd_fus.groupby(['基因', 'fusion_gene_name'])['value'].apply(lambda x: Tools.cluster_columns(x, sep=';'))
                    pd_fus = pd_fus.reset_index()
                    header_tmp = list(pd_fus)

                    list_value_format = []
                    for i in list(pd_fus['value']):
                        if ';' in i:
                            value_info = str(i).split(';')
                            value_rna = [x for x in value_info if 'FFPM' in x][0]
                            value_dna = [x for x in value_info if '%' in x][0]
                            list_value_format.append(f"DNA:{value_dna};RNA:{value_rna}")
                        else:
                            list_value_format.append(i)
                    pd_fus['value'] = list_value_format
                    pd_fus = pd_fus[header_tmp]

                pd_fus['变异基因'] = pd_fus['基因']
                pd_fus['外显子/内含子'] = '-'
                pd_fus['变异位点'] = pd_fus['fusion_gene_name']
                pd_fus[sample_type] = pd_fus['value']
                pd_fus['变异类型'] = '基因重排'
                pd_fus['数据库收录(COSMIC)'] = '-'
                pd_fus['SIFT'] = '-'
                pd_fus['Polyphen2'] = '-'
                pd_fus['REVEL'] = '-'
                pd_fus['MutationTaster'] = '-'
                pd_fus = pd_fus[header]
                pd_fus = pd_fus.drop_duplicates()
            else:
                pd_fus = pd.DataFrame(columns=header)
            list_fus.append(pd_fus)

        pd_mut_merge_no_db = pd.merge(list_mut[0], list_mut[1],
                                      on=['变异基因', '外显子/内含子', '变异位点', '变异类型', '数据库收录(COSMIC)', 'SIFT',
                                          'Polyphen2', 'REVEL', 'MutationTaster'], how='outer')
        pd_mut_merge_db = pd.merge(list_mut_db[0], list_mut_db[1],
                                   on=['变异基因', '外显子/内含子', '变异位点', '变异类型', '数据库收录(COSMIC)', 'SIFT',
                                       'Polyphen2', 'REVEL', 'MutationTaster'], how='outer')

        if not pd_mut_merge_no_db.empty:
            list_type = list(set(pd_mut_merge_no_db) & {'FFPE', 'ctDNA'})
            if len(list_type) == 1:
                if str(list_type[0]).lower() != 'ctdna':
                    """ 组织（新鲜组织/石蜡切片/胸水）：cutoff≥2% """
                    pd_mut_merge_no_db = pd_mut_merge_no_db[pd_mut_merge_no_db['FFPE'] >= 0.02]
                else:
                    """ 血液（cfDNA/ctDNA）：cutoff≥0.5% """
                    pd_mut_merge_no_db = pd_mut_merge_no_db[pd_mut_merge_no_db['ctDNA'] >= 0.005]

            elif len(list_type) == 2:
                pd_mut_merge_no_db = pd_mut_merge_no_db[(pd_mut_merge_no_db['FFPE'] >= 0.02) | (pd_mut_merge_no_db['ctDNA'] >= 0.005)]

            for each_type in list_type:
                pd_mut_merge_no_db[each_type] = pd_mut_merge_no_db[each_type].map(lambda x: format(x, '.1%'))
                pd_mut_merge_no_db[each_type] = pd_mut_merge_no_db[each_type].replace('nan%', '-')

        pd_mut_merge = pd.concat([pd_mut_merge_no_db, pd_mut_merge_db], axis=0).drop_duplicates()

        pd_cnv_merge = pd.merge(list_cnv[0], list_cnv[1],
                                on=['变异基因', '外显子/内含子', '变异位点', '变异类型', '数据库收录(COSMIC)', 'SIFT', 'Polyphen2',
                                    'REVEL', 'MutationTaster'], how='outer')

        pd_fus_merge = pd.merge(list_fus[0], list_fus[1],
                                on=['变异基因', '外显子/内含子', '变异位点', '变异类型', '数据库收录(COSMIC)', 'SIFT', 'Polyphen2',
                                    'REVEL', 'MutationTaster'], how='outer')

        out_unit_b = pd.concat([pd_mut_merge, pd_cnv_merge, pd_fus_merge], axis=0)
        out_unit_b = out_unit_b[~(out_unit_b['变异基因'] == '-')]
        out_unit_b = out_unit_b.fillna(value='-')
        out_unit_b = out_unit_b.sort_values(by=['变异基因', self.type1, self.type2], ascending=[True, False, False])
        # out_unit_b = out_unit_b.head(100)
        out_unit_b = out_unit_b[['变异基因', '外显子/内含子', '变异位点', self.type1, self.type2, '变异类型',
                                 '数据库收录(COSMIC)', 'SIFT', 'Polyphen2', 'REVEL', 'MutationTaster']]
        out_unit_b = out_unit_b.drop_duplicates()
        out_unit_b.to_csv(f'{self.outdir}/附录B.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/附录B.txt')

    def unit_c(self):
        """ 本次基因检测结果总览（胚系变异） """
        """ 基因	转录本ID	外显子/内含子	核苷酸改变	氨基酸改变	纯合/杂合	致病风险	遗传方式 """
        pd_mut = pd.read_csv(self.mutation_file1, sep='\t', low_memory=False, keep_default_na=False)
        pd_mut = self.filter_little_panel(pd_mut)
        pd_mut['外显子'] = pd_mut['核酸改变'].apply(lambda x: Tools.extrac_exon(x))
        pd_mut = pd_mut.apply(lambda x: Tools.rm_exon(x), axis=1)
        pd_mut = pd_mut.fillna(value='')

        def unit_c_convert(x):
            if x['氨基酸改变'] in ['', ' ', 'NA']:
                return '(-)'
            else:
                return x['氨基酸改变']

        pd_mut['氨基酸改变'] = pd_mut.apply(lambda x: unit_c_convert(x), axis=1)

        # (1) 纯合、杂合、半合子 判断
        pd_mut['突变状态'] = ''
        pd_mut["突变状态"] = pd_mut.apply(lambda x: Tools.mut_status(x, self.gender, self.normal_id), axis=1)

        # (2) 致病性、可能致病性、NA 判断  去掉其他无意义注释
        pd_mut['pathogenicity'] = ''
        pd_mut["pathogenicity"] = pd_mut.apply(lambda x: Tools.get_benign(x), axis=1)

        # 检不出结果可以解开下面两行这样测试
        # pd_mut['突变级别'] = '1'
        # pd_mut['pathogenicity'] = '致病'
        # pd_mut.to_csv('遗传性肿瘤基因检测结果temp.txt', sep='\t', index=False)

        # (3) 进行'遗传性突变'的过滤
        pd_mut['结论'] = pd_mut.apply(
            lambda x: Tools.get_conclusion(x, self.product, panel_id=self.tumor_id1, wes_id=self.tumor_wes_id1), axis=1)
        # pd_mut.to_csv('过滤前良性注释.txt', sep='\t', index=False)
        pd_mut = pd_mut[
            (pd_mut['结论'] == '遗传性突变') &
            pd_mut['突变级别'].isin(['1', '2']) &
            ~pd_mut['突变类型'].isin(['非编码区突变', '同义突变']) &
            (~pd_mut['pathogenicity'].isin(['良性', '可能良性']))]

        # pd_mut.to_csv('过滤后.txt', sep='\t', index=False)
        header = ['变异基因', '外显子/内含子', '变异位点', '突变状态', '变异类型', '数据库收录(COSMIC)',
                  'SIFT', 'Polyphen2', 'REVEL', 'MutationTaster']

        if not pd_mut.empty:
            # pd_mut['外显子'] = pd_mut['核酸改变'].apply(lambda x: self.extrac_exon(x))
            # pd_mut = pd_mut.apply(lambda x: self.rm_exon(x), axis=1)
            pd_mut['变异基因'] = pd_mut['基因']
            pd_mut['外显子/内含子'] = pd_mut['外显子']
            pd_mut['变异位点'] = pd_mut['mRNA变体名称(NM)'] + ':;' + pd_mut['核酸改变'] + ':;' + pd_mut['氨基酸改变']
            pd_mut['变异类型'] = pd_mut['突变类型']
            pd_mut['数据库收录(COSMIC)'] = pd_mut.apply(lambda x: Tools.convet_cosmic(x), axis=1)
            pd_mut['SIFT'] = pd_mut['SIFT危害性预测'].apply(lambda x: x.split('(')[0])
            pd_mut['Polyphen2'] = pd_mut['Polyphen2_HDIV危害性预测'].apply(lambda x: x.split('(')[0])
            pd_mut['REVEL'] = pd_mut['REVEL危害性预测'].apply(lambda x: x.split('(')[0])
            pd_mut['MutationTaster'] = pd_mut['MutationTaster危害性预测'].apply(lambda x: x.split('(')[0])

            # pd_mut = pd.merge(pd_mut, pd_tumor_genetic_predisposed_gene_list, on='基因', how='inner')
            pd_mut = pd_mut[header]
            pd_mut = pd_mut.sort_values(by=['变异基因', '变异位点'], ascending=[True, True])
            pd_mut = pd_mut.drop_duplicates()
            panel39_code = ['NCPC1', 'NCPC2', 'NCPC3']
            if self.project_id in panel39_code:
                panel39_gene_frame = pd.read_excel(panel39_germline_genes, sheet_name='panel39')
                pd_mut = pd.merge(pd_mut, panel39_gene_frame, left_on='变异基因', right_on='基因', how='inner')
            pd_mut = pd_mut.head(100)
        else:
            pd_mut = pd.DataFrame(columns=header)

        pd_mut = pd_mut.replace('NA', '-')
        pd_mut.to_csv(f'{self.outdir}/附录C.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/附录C.txt')

    def unit_d(self):
        """
        基因描述及研究解读
        不过药物
        """
        pd_mut1, pd_mut2 = self.filter_mut(if_do_targeting_site_drug=True)
        pd_cnv1, pd_cnv2 = self.filter_cnv(if_do_targeting_site_drug=True)
        pd_fus1, pd_fus2 = self.filter_fus(if_do_targeting_site_drug=True)

        list_for_f = []
        if self.single_flag:
            for_list = [[pd_mut1, pd_cnv1, pd_fus1, self.type1]]
        else:
            for_list = [[pd_mut1, pd_cnv1, pd_fus1, self.type1],
                        [pd_mut2, pd_cnv2, pd_fus2, self.type2]]

        for pd_mut_tmp, pd_cnv_tmp, pd_fus_tmp, sample_type in for_list:
            pd_mut = pd_mut_tmp[['基因', '药物']].copy()
            pd_cnv = pd_cnv_tmp[['基因', '药物']].copy()
            pd_fus = pd_fus_tmp[['基因', '药物']].copy()

            if os.path.isfile(self.rna_fusion_file) and str(sample_type).lower() != 'ctdna':
                header_tmp = list(pd_fus)
                pd_rna_fusion = self.filter_rna(if_do_targeting_site_drug=True)
                pd_fus = pd.concat([pd_fus, pd_rna_fusion], axis=0)
                pd_fus = pd_fus[header_tmp].drop_duplicates()

            # pd_mut_HasDrug = pd_mut   # [~pd_mut['药物'].isin([''])]
            # pd_cnv_HasDrug = pd_cnv   # [~pd_cnv['药物'].isin([''])]
            # pd_fus_HasDrug = pd_fus   # [~pd_fus['药物'].isin([''])]

            out_unit_d = pd.concat([pd_mut, pd_cnv, pd_fus], axis=0, ignore_index=True)
            out_unit_d = out_unit_d[['基因', '药物']].drop_duplicates()
            list_for_f.append(out_unit_d)

        pd_solid_tumor_gene_annotation = pd.read_excel(self.solid_tumor_gene_annotation)
        pd_solid_tumor_clinical_trials = pd.read_excel(self.solid_tumor_ClinicalTrials)
        pd_solid_tumor_clinical_trials = pd_solid_tumor_clinical_trials.fillna(value='')

        out_unit_for_e = pd.concat(list_for_f, axis=0)
        out_unit_for_e = out_unit_for_e.drop_duplicates()

        # 20210908: 肠癌新增(KRAS/NRAS/BRAF三基因阴性时在肠癌中需展示)
        special_gene_frame = out_unit_for_e[
            ((out_unit_for_e['基因'] == 'KRAS') | (out_unit_for_e['基因'] == 'NRAS') | (out_unit_for_e['基因'] == 'BRAF'))]
        if self.cancer_flag == '肠' and special_gene_frame.empty:
            special_info = list()
            for gene in ['KRAS', 'BRAF', 'NRAS']:
                special_info.append([gene, '西妥昔单抗'])
                special_info.append([gene, '帕尼单抗'])
            out_unit_for_e = out_unit_for_e.append(pd.DataFrame(special_info, columns=list(out_unit_for_e)),
                                                   ignore_index=True)
        out_unit_for_e.to_csv('out_unit_for_e.txt', sep='\t', index=False)
        pd_gene_anno = pd.merge(out_unit_for_e, pd_solid_tumor_gene_annotation, on='基因', how='left')

        pd_gene_anno_drug_anno = pd.merge(pd_gene_anno, pd_solid_tumor_clinical_trials, left_on='药物', right_on='相关药物', how='left')
        pd_mut_gene_anno_drug_anno = pd_gene_anno_drug_anno.fillna(value='')

        header = ['基因', '基因简介', f'临床价值({self.cancer_flag})', '研究类型', '研究癌种', '相关药物', '药物类型', '研究说明']
        pd_mut_gene_anno_drug_anno = pd_mut_gene_anno_drug_anno[header]
        # pd_mut_gene_anno_drug_anno = pd_mut_gene_anno_drug_anno[~(pd_mut_gene_anno_drug_anno['研究类型'].isin(['', ' ']))]
        pd_mut_gene_anno_drug_anno = pd_mut_gene_anno_drug_anno.drop_duplicates()
        pd_mut_gene_anno_drug_anno.to_csv(f'{self.outdir}/附录D.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/附录D.txt')
        return out_unit_for_e

    def unit_e(self, out_unit_for_e):
        """ 相关临床试验 """
        gene_list = list(set(out_unit_for_e['基因']))
        db_linchuangshiyan = pd.read_excel(self.solid_tumor_linchuangshiyan)
        out_unit_e = db_linchuangshiyan[db_linchuangshiyan['基因'].isin(gene_list)]
        out_unit_e.to_csv(self.outdir + '附录E.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/附录E.txt')

    def unit_f(self, out_unit_for_e):
        """ 有药和无药都要 """
        gene_list = set(out_unit_for_e['基因'])
        "信号通路名称（英文）	信号通路名称（中文）	通路介绍	靶向通路药物	检测基因"
        if self.project_id in ['NCPM5', 'NCPM6', 'NCPM7']:
            db_solid_tumor_pathway = pd.read_excel(self.solid_tumor_pathway_393)
        else:
            db_solid_tumor_pathway = pd.read_excel(self.solid_tumor_pathway)

        db_solid_tumor_pathway['通路介绍'] = db_solid_tumor_pathway['通路介绍'].apply(lambda x: str(x).replace('\n', ''))

        def convert_unit_f(x):
            my_gene_list = set(str(x['检测基因']).split())
            intersection = my_gene_list & gene_list
            if len(intersection) == 0:
                return ''
            else:
                return ' '.join(intersection)

        db_solid_tumor_pathway['实际检测'] = ''
        db_solid_tumor_pathway['实际检测'] = db_solid_tumor_pathway.apply(lambda x: convert_unit_f(x), axis=1)
        db_solid_tumor_pathway.to_csv(f'{self.outdir}/附录F.txt', sep='\t', index=False)
        self.result_file_list.append(f'{self.outdir}/附录F.txt')


if __name__ == '__main__':
    print('-----------------------------------------')
    print('\tStart...')
    argv = my_parser()
    ob_run = CancerFilter(argv)
    ob_run.unit_1_1()
    ob_run.unit_1_2()
    ob_run.unit_1_3()

    # NCCN

    if argv.cancer == '脑胶质':
        ob_run.unit_nccn_glioma()
    else:
        pd_solid_tumor_NCCN_gene = pd.read_excel(argv.solid_tumor_NCCN_gene, sheet_name=None)
        if argv.cancer in list(pd_solid_tumor_NCCN_gene):
            ob_run.unit_nccn_lung()
        else:
            nccn = pd.DataFrame(columns=['基因', '变异类型', '检测结果', '临床意义'])
            nccn.to_csv(f'{argv.outdir}/NCCN.txt', sep='\t', index=False)

    if '前列腺' in argv.cancer and argv.project_id not in ['NCPM5', 'NCPM6', 'NCPM7', 'NCP393', 'NCP393c']:
        ob_run.unit_2_0()

    ob_run.unit_2_1()
    ob_run.unit_2_2()
    pd_all_gene_type = ob_run.unit_4()
    ob_run.unit_4_conclusion(pd_all_gene_type)
    ob_run.unit_5()
    ob_run.unit_b()
    ob_run.unit_c()
    out_unit_for_e = ob_run.unit_d()
    ob_run.unit_e(out_unit_for_e)
    ob_run.unit_f(out_unit_for_e)

    dict_result = {os.path.basename(x).split('.')[0]: os.path.realpath(x) for x in ob_run.result_file_list}
    response_dict = json.dumps(dict_result, indent=4, ensure_ascii=False, sort_keys=False)
    open(f'{os.path.realpath(argv.outdir)}/filter_result.json', 'w').write(response_dict)
    print('\tEnd...')
    print('-----------------------------------------')

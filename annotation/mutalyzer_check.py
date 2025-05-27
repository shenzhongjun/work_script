#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
对生产流程报出的结果调用mutalyzer进行校正
"""
import argparse
import requests
import os
import re
import json
import pandas as pd

with open('/mnt/share02/zhouyj/script/annotation/make_mutation_table.json') as f:
    json_content = json.load(f)


def process_row(row):
    # 去掉括号如c.2740(exon17)C>T；去掉多余的c.如c.3192_c.3193delGCinsAT
    row['核酸改变'] = re.sub(r'\(.*?\)', '', row['核酸改变'])
    row['核酸改变'] = row['核酸改变'].replace('_c.', '_') if '_c.' in row['核酸改变'] else row['核酸改变']
    description = f"{row['mRNA变体名称(NM)']}:{row['核酸改变']}"
    normalized_id, normalized_c, normalized_p = mutalyzer_normalize(description)
    return pd.Series([normalized_id, normalized_c, normalized_p])


def mutalyzer_normalize(description):
    # 基本的APIURL
    base_url = "https://mutalyzer.nl/api"
    # 变异校正计算端点
    endpoint = f"/normalize/{description}"
    # 发送GET请求
    url = base_url + endpoint
    response = requests.get(url)

    # 检查响应
    if '+' in description or '-' in description:
        print(f"Error:不支持剪切位点校正 {description}")
        return '', '', ''
    elif response.status_code == 200:
        print("Success:", description)
        # print(response.json())
        protein_description = response.json()['protein']['description']
        p_nm = protein_description.split('(')[0]
        hgvsp = 'p.' + re.search(r"p\.\((\w+\d+.+)\)", protein_description).group(1) if '?' not in protein_description else 'p.?'
        if 'chromosomal_descriptions' in response.json().keys():
            chr_description = response.json()['chromosomal_descriptions'][0]['c']
            hgvsc = chr_description.split(':')[1]
            c_nm = re.search(r"\((NM_\d+\.\d+)\)", chr_description).group(1)
        else:
            chr_description = response.json()['normalized_description']
            c_nm, hgvsc = chr_description.split(':')
        if p_nm == c_nm:
            # print(hgvsc, hgvsp, sep='\n')
            return c_nm, hgvsc, hgvsp
        else:
            print(f'{description}校正失败，c.和p.的转录本不一致')
            return '', f'{c_nm}:{hgvsc}', f'{p_nm}:{hgvsp}'
    else:
        print(f"Error:{response.status_code} {description}")
        return '', '', ''


def normalize_df(df, out):
    if df.empty:
        with open(out, 'w') as f:
            f.write("阴性报告")
    else:
        df[["HGVS转录本", "HGVS核酸改变", "HGVS氨基酸改变"]] = df.apply(process_row, axis=1)
        df['HGVS氨基酸改变2'] = df['HGVS氨基酸改变'].replace(json_content['amino_acids_dict'], regex=True)
        df.to_csv(out, sep='\t', index=False)


def normalize_blood(args):
    print(f"血肿报出结果校正：\n输入 {args.known} {args.manual} {args.germline}\n输出 {args.out}")
    known_df = pd.read_csv(args.known, sep='\t').iloc[:, :6]
    manual_df = pd.read_csv(args.manual, sep='\t').iloc[:, :6]
    germline_df = pd.read_csv(args.germline, sep='\t').iloc[:, :6]
    merge_df = pd.concat([known_df, manual_df, germline_df], ignore_index=True)
    merge_df = merge_df.rename(columns={"转录本": "mRNA变体名称(NM)"})
    normalize_df(merge_df, args.out)


def normalize_solid(args):
    print(f"实体瘤报出结果校正：\n输入 {args.mut}\n输出 {args.out}")
    df = pd.read_csv(args.mut, sep='\t').iloc[:, :9]
    normalize_df(df, args.out)


def normalize_child(args):
    print(f"儿童实体瘤报出结果校正：\n输入 {args.ks} {args.ms} {args.kg} {args.mg}\n输出 {args.out}")
    ks_df = pd.read_csv(args.ks, sep='\t').iloc[:, :8]
    ms_df = pd.read_csv(args.ms, sep='\t').iloc[:, :8]
    kg_df = pd.read_csv(args.kg, sep='\t').iloc[:, list(range(0, 7)) + [8]]
    mg_df = pd.read_csv(args.mg, sep='\t').iloc[:, list(range(0, 7)) + [8]]
    merge_df = pd.concat([ks_df, ms_df, kg_df, mg_df], ignore_index=True)

    normalize_df(merge_df, args.out)


def get_args():
    parser = argparse.ArgumentParser(description="对生产流程报出的结果调用mutalyzer进行校正")
    # 定义子命令
    subparsers = parser.add_subparsers(dest='mode', help='选择血肿或实体瘤模式')
    # 血肿
    parser_a = subparsers.add_parser('blood', help='对血肿结果进行校正')
    parser_a.add_argument('--known', '-k', type=str, required=True, help='血肿known_somatic_mutation文件')
    parser_a.add_argument('--manual', '-m', type=str, required=True, help='血肿manual_somatic_mutation文件')
    parser_a.add_argument('--germline', '-g', type=str, required=True, help='血肿know_germline_mutation文件')
    parser_a.add_argument('--out', '-o', type=str, default='HGVS_normalize.xls', help='输出文件')
    # 实体瘤
    parser_b = subparsers.add_parser('solid', help='对实体瘤结果进行校正')
    parser_b.add_argument('--mut', '-m', type=str, required=True, help='实体瘤merge_new.txt_filter_mut文件')
    parser_b.add_argument('--out', '-o', type=str, default='HGVS_normalize.xls', help='输出文件')
    # 儿童实体瘤
    parser_a = subparsers.add_parser('child', help='对儿童实体瘤blood_a_base结果进行校正')
    parser_a.add_argument('--ks', type=str, required=True, help='儿童实体瘤known_somatic_mutation文件')
    parser_a.add_argument('--ms', type=str, required=True, help='儿童实体瘤manual_somatic_mutation文件')
    parser_a.add_argument('--kg', type=str, required=True, help='儿童实体瘤germline_known_mutation文件')
    parser_a.add_argument('--mg', type=str, required=True, help='儿童实体瘤germline_manual_mutation文件')
    parser_a.add_argument('--out', '-o', type=str, default='HGVS_normalize.xls', help='输出文件')

    return parser.parse_args()


if __name__ == "__main__":
    # 设置环境变量使节点可联网
    os.environ['http_proxy'] = "http://192.168.1.10:3128"
    os.environ['https_proxy'] = "http://192.168.1.10:3128"

    args = get_args()
    if args.mode == 'blood':
        normalize_blood(args)
    elif args.mode == 'solid':
        normalize_solid(args)
    elif args.mode == 'child':
        normalize_child(args)
    else:
        print("参数错误.")

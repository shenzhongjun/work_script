"""
以订单号+样本名为输入，获得保存在集群的merge表并进整理，因为有新旧版差异，需要写脚本处理
"""
import os
import glob
import re
import argparse
import subprocess
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='以订单号+样本名为输入，获得保存在集群的merge表并进整理')
    parser.add_argument('--orders', help='订单号-样本名列表，查找对应的merge表', required=True)
    parser.add_argument('--groups', help='订单号-分组列表，设置分组', required=True)
    parser.add_argument('--bak', help='merge表备份路径', default='/mnt/share03/tumor_pipeline/Somatic/DNA/Other/New_report/')
    parser.add_argument('--out', help='输出文件路径', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    argv = get_args()
    bak_path = argv.bak
    out_path = argv.out
    with open(argv.orders) as f1, open(argv.groups) as f2, open(f'run.sh', 'w') as w:
        f2.readline()
        order_dict = {}
        for i in f2:
            order, name, group, gender, age = i.strip().split('\t')
            order_dict[order] = group
        for i in f1:
            order, sample = i.strip().split('\t')
            if glob.glob(f'{bak_path}/{order}_{sample}_NBCW_*/{sample}.merge.txt'):
                file_path = glob.glob(f'{bak_path}/{order}_{sample}_NBCW_*/{sample}.merge.txt')[0]
                print(f'新---{order}', file_path)
                subprocess.call(f'ln -sf {file_path} {out_path}/{order}_{order_dict[order]}_merge.txt', shell=True)
            elif glob.glob(f'{bak_path}/{order}_{sample}_NBCW_*/{order}_{sample}_NBCW_*_Somatic.xlsx'):
                file_path = glob.glob(f'{bak_path}/{order}_{sample}_NBCW_*/{order}_{sample}_NBCW_*_Somatic.xlsx')[0]
                print(f'旧---{order}', file_path)
                if not os.path.exists(f'{out_path}/{order}_{sample}_merge.txt'):
                    df = pd.read_excel(file_path, sheet_name=0).rename(columns={'突变可靠性\n': '突变可靠性'})
                    df.to_csv(f'{out_path}/{order}_{order_dict[order]}_merge.txt', sep='\t', index=False)
            else:
                print(f'未发现！{order}')

            w.write(f'python merge_table_to_maf.py --wes_merge {out_path}/{order}_{order_dict[order]}_merge.txt '
                    f'--prefix {order}_{order_dict[order]} --wes_id {sample} --out mafs &\n')














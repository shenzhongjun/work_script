#!/mnt/share02/zhouyj/software/miniconda3/envs/python3.9/bin/python
# -*- coding: utf-8 -*-

"""
快速查看fastq.gz文件前200条序列的Q30
用法: q30 xx.fastq.gz
"""

import gzip
import sys
import numpy as np


def extract_quality_values(fastq_file, num_sequences=200):
    """从FASTQ文件中提取前num_sequences条序列的质量值"""
    quality_values = []  # 用于存储质量值行
    with gzip.open(fastq_file, 'rt') as f:  # 以文本模式读取.gz文件
        count = 0
        for line_number, line in enumerate(f, start=1):
            if line_number % 4 == 0:  # FASTQ文件的每第4行是质量值
                quality_values.append(line.strip())
                count += 1
                if count == num_sequences:  # 提取到指定数量的序列后停止
                    break
    return ''.join(quality_values)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python script.py <fastq.gz文件路径>")
        sys.exit(1)

    fastq_file = sys.argv[1]
    try:
        # 提取质量值
        quality_string = extract_quality_values(fastq_file, num_sequences=10)

        # 计算
        quality_values = [ord(c) - 33 for c in quality_string]  # 将字符转换为质量值
        q30_count = sum(1 for q in quality_values if q >= 30)  # 统计Q30数量
        q30_percentage = (q30_count / len(quality_values)) * 100  # 计算Q30百分比
        print(f"Q30百分比: {q30_percentage:.2f}%")

        min_quality = min(quality_values)  # 最小质量值
        max_quality = max(quality_values)  # 最大质量值
        mean_value = np.mean(quality_values)    # 计算平均值
        median_value = np.median(quality_values)    # 计算中位数
        print(f"最小质量值: {min_quality}")
        print(f"最大质量值: {max_quality}")
        print(f"平均值: {mean_value:.2f}")
        print(f"中位数: {median_value}")
    except FileNotFoundError:
        print(f"错误: 文件 {fastq_file} 不存在！")
    except Exception as e:
        print(f"发生错误: {e}")
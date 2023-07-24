#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
CMD:
/home/zhouyj/anaconda3/envs/rseqc/bin/geneBody_coverage.py \
    -i Aligned.out.new.bam \
    -r /mnt/share02/zhouyj/database/RNA/RSeQC_Files/hg19_RefSeq.bed \
    -f png -o rna_B_4_3.coverage

/home/zhouyj/anaconda3/envs/rseqc/bin/read_distribution.py \
    -i Aligned.out.new.bam \
    -r /mnt/share02/zhouyj/database/RNA/RSeQC_Files/hg19_RefSeq.bed

/home/zhouyj/anaconda3/envs/samtools_1_15/bin/samtools reheader \
    -c "sed  's,^@RG.*,@RG\tID:GRPundef\tSM:rna_P_0_1\tPL:Illumina,g'" Aligned.out.bam | \
sambamba sort -t 32 -m 40G -o Aligned.out.new.bam /dev/stdin && \
gatk CollectInsertSizeMetrics -I Aligned.out.new.bam -O AFJ525.insert_size_histogram.txt -H AFJ525.insert_size_histogram.pdf -M 0.5
"""

import argparse
import os
import sys

if __name__ == "__main__":
    input_path = os.path.abspath(sys.argv[1])
    output_path = f'{os.path.splitext(input_path)[0]}.txt'
    with open(input_path) as f, open(output_path, 'w') as w:
        lines = []
        total_ratio = 0
        for line in f:
            if not line.startswith('Total') and not line.startswith('='):
                if line.startswith('Group'):
                    header = '\t'.join(line.strip().split()) + '\tTag_percent\n'
                    w.write(header)
                else:
                    lines.append(line)
                    # print(line.split())
                    total_ratio += float(line.split()[3])
        for line in lines:
            percent = f'{float(line.split()[3])/total_ratio*100:.2f}%'
            line = '\t'.join(line.strip().split()) + f'\t{percent}\n'
            w.write(line)

    # args = get_args()
    # out = os.path.abspath(args.out)



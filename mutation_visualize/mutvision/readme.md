## 说明

```bash
# 使用示例
python mutationVisualize.py \
--bam /bioinfo/commercial/cancer/BN201229_1/200304973T1/Mapping/200304973T1.tbr.bam \
--mutlist /bioinfo/commercial/cancer/BN201229_1/200304973T1/rpt/200304973T1.annot.vcf \
--ref /bioinfo/home/acorndx07/DB/ucsc.hg19.fasta \
--tlenfile tlen.txt --html out.html
```

anno.vcf文件内容：
```txt
chr10   89720633    89720633    -   T
chr7    140498360   140498360   T   -
```

- 默认bpQ阈值为：30，在html中低于该阈值的碱基会以灰色呈现
- 默认mapQ阈值为：40，主要用于从bam中提出目标reads

在html可视化结果中:

- <font style="background:#00ff00">绿底色</font>表示错配，
- <font style="background:#8B6914" color='#C4C4C4'>土底色</font>表示softclip碱基，
- <font style="background:#EEEE00">黄底色</font>表示插入碱基，
- <font style="background:#C0D9D9">浅蓝底色</font>表示该reads的变异与输入变异有部分不一致
- <font color="#C0C0C0">灰色字体</font>表示碱基质量值低，
- <font color="#0165FC">蓝色字体</font>表示该reads比对到正链，
- <font color="#EF4026">橘红色字体</font>表示该reads比对到负链，
- ‘*’表示碱基缺失，
- ‘-’表示reads1，
- ‘=’表示reads2。


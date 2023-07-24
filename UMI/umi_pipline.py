#!/mnt/share01/tools/bin/python
# -*- coding: UTF-8 -*-

"""
UMI纳昂达+Sentieon分析流程
python umi_pipline.py --sample_list sample_list.txt
"""
__author__ = "Fei Xue"
__email__ = "xuefei@zhuanhuayixue.org"
__version__ = "1.0.2"
__date__ = "8/19/21 17:08 "
__copyright__ = "Copyright (C) 2021 All rights reserved"

import os
import time
import argparse
import yaml

localtime = time.strftime("%Y_%m_%d_%H%M%S")


def my_parser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--config', help="配置文件路径", default="/mnt/share02/wangwp/scripts/project/UMI/UMI_pipline.yaml")
    parser.add_argument('--sample_list', help="样本信息，包含样本名、fq1、fq2", required=True)
    parser.add_argument('--sample_ng', help="上样量（ng），默认30ng", type=int, default=30)
    parser.add_argument('--path', help="工作目录", default=f'{os.getcwd()}')
    parser.add_argument('--PLATFORM_UNIT', help="测序仪", choices=['HiseqX10', 'T7'], default='T7')
    parser.add_argument('--PLATFORM', help="测序平台", choices=['Illumina', 'DNBSEQ'], default='DNBSEQ')
    parser.add_argument('--read_structure', help="UMI格式，Illumina：'3M2S+T', BGI：'3M4S+T'",
                        choices=['3M2S+T', '3M4S+T'], default='3M4S+T')
    parser.add_argument('--duplex_umi', help='UMI添加形式，默认为双端UMI', default=True)
    parser.add_argument('--min_reads', help="形成一致性双链的过滤标准，默认为6 3 3",
                        choices=['2 1 1', '6 3 3'], default='6 3 3')
    parser.add_argument('--min_depth', help='sentieon过滤最低总深度，默认为100', type=int, default=100)
    parser.add_argument('--normal_sample', help="对照样本名")
    parser.add_argument('--normal_bam', help="对照样本bam文件路径")
    return parser.parse_args()


class Run:

    def __init__(self):
        self.argv = my_parser()
        self.sample_name, self.fq1, self.fq2 = open(self.argv.sample_list).readline().strip().split('\t')
        self.conf = yaml.load(open(self.argv.config), Loader=yaml.FullLoader)
        self.bin = self.conf["bin"]
        self.target = self.conf["database"]["bed"]
        self.ref = self.conf["database"]["ref"]
        self.dbsnp = self.conf["database"]["dbsnp"]
        self.sample_ng = self.argv.sample_ng
        self.path = self.argv.path
        self.PLATFORM_UNIT = self.argv.PLATFORM_UNIT
        self.PLATFORM = self.argv.PLATFORM
        self.read_structure = self.argv.read_structure
        self.duplex_umi = self.argv.duplex_umi
        self.min_reads = self.argv.min_reads
        self.min_depth = self.argv.min_depth
        self.normal_sample = self.argv.normal_sample
        self.normal_bam = self.argv.normal_bam
        self.job_file = os.path.join(self.path, f'umi_{localtime}.job')
        self.job_file_open = open(self.job_file, 'w')

    def step0(self):
        sample_path = os.path.join(self.path, 'step0.cutadapt')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/bin/cutadapt \
    -j 8 \
    -a CTCTCAGTACGTCAGCAGTTNNNNNNNNNNCAACTCCTTGGCTCACAGAAC \
    -A GCATGGCGACCTTATCAGNNNNNNNNNNTTGTCTTCCTAAGACCGCTTGG \
    -o {sample_path}/{self.sample_name}.clean.R1.fastq.gz \
    -p {sample_path}/{self.sample_name}.clean.R2.fastq.gz \
    -m 50 \
    --match-read-wildcards \
    {self.fq1} \
    {self.fq2}
"""
        job_name = f'step0_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = 'None'
        job_write(self.job_file_open, script, job_name, '8', job_script, pre_jobs, sample_path, 'ready')
        return f'{sample_path}/{self.sample_name}.clean.R1.fastq.gz', f'{sample_path}/{self.sample_name}.clean.R2.fastq.gz'

    def fastqc(self):
        sample_path = os.path.join(self.path, 'step0.cutadapt')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/bin/fastqc \
    -o  {sample_path} \
    -t 8 -q \
    {sample_path}/{self.sample_name}.clean.R1.fastq.gz \
    {sample_path}/{self.sample_name}.clean.R2.fastq.gz
        """
        job_name = f'step0_fastqc_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step0_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '8', job_script, pre_jobs, sample_path, 'ready')

    def fastqQC(self):
        """
        用于统计fastq的数据量，gc含量，Q30%
        """
        sample_path = os.path.join(self.path, 'QC')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/software/gzGCqual_Q33 \
    {self.fq1} \
    {self.fq2} \
    {sample_path}/{self.sample_name}.GC_content \
    {sample_path}/{self.sample_name}.qulity \
    {sample_path}/{self.sample_name}.totstat
        """
        job_name = f'step0_fastqQC_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step0_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '8', job_script, pre_jobs, sample_path, 'ready')
        return f"{sample_path}/{self.sample_name}.totstat"

    def step1(self, fq1, fq2):
        sample_path = os.path.join(self.path, 'step1.FastqToSam')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar \
    FastqToSam \
    FASTQ={fq1} \
    FASTQ2={fq2} \
    OUTPUT={sample_path}/{self.sample_name}.uBAM \
    READ_GROUP_NAME={self.sample_name} \
    SAMPLE_NAME={self.sample_name} \
    LIBRARY_NAME={self.sample_name} \
    PLATFORM_UNIT={self.PLATFORM_UNIT} \
    PLATFORM={self.PLATFORM} \
    STRIP_UNPAIRED_MATE_NUMBER=true \
    ALLOW_AND_IGNORE_EMPTY_LINES=true
"""
        job_name = f'step1_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step0_{self.sample_name}, step0_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '15', job_script, pre_jobs, sample_path, 'ready')
        step1_bam = os.path.join(sample_path, f'{self.sample_name}.uBAM')
        return step1_bam

    def step2(self, step1_bam):
        sample_path = os.path.join(self.path, 'step2.ExtractUmisFromBam')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java  \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/fgbio/fgbio.jar \
    ExtractUmisFromBam \
    --input={step1_bam} \
    --output={sample_path}/{self.sample_name}.umi.uBAM \
    --read-structure={self.read_structure} {self.read_structure} \
    --single-tag=RX \
    --molecular-index-tags=ZA ZB
rm {step1_bam}
"""
        job_name = f'step2_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step1_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '15', job_script, pre_jobs, sample_path, 'ready')
        step2_bam = os.path.join(sample_path, f'{self.sample_name}.umi.uBAM')
        return step2_bam

    def step3(self, step2_bam):
        sample_path = os.path.join(self.path, 'step3.SamToFastq_bwa_MergeBamAlignment')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar \
    SamToFastq \
    I={step2_bam} \
    F=/dev/stdout \
    INTERLEAVE=true \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/bin/bwa mem -p -t 20 {self.ref} /dev/stdin \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar \
    MergeBamAlignment \
    UNMAPPED={step2_bam} \
    ALIGNED=/dev/stdin \
    O={sample_path}/{self.sample_name}.umi.merged.uBAM \
    R={self.ref} \
    SO=coordinate \
    ALIGNER_PROPER_PAIR_FLAGS=true \
    MAX_GAPS=-1 \
    ORIENTATIONS=FR \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true
rm {step2_bam}
"""
        job_name = f'step3_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step2_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '15', job_script, pre_jobs, sample_path, 'ready')
        step3_bam = os.path.join(sample_path, f'{self.sample_name}.umi.merged.uBAM')
        return step3_bam

    def bamdst(self, bam, outdir, bed, job_name, pre_job):
        os.makedirs(f"{outdir}/depthQC", exist_ok=True)
        script = fr"""
/mnt/share01/tools/bamdst/bamdst-master/bamdst \
    -p {bed} \
    -o {outdir}/depthQC \
    --flank 150 \
    {bam}
    """
        job_name = f'{job_name}_{self.sample_name}'
        job_script = f'{outdir}/{job_name}.sh'
        job_write(self.job_file_open, script, job_name, '8', job_script, pre_job, outdir, 'ready')

    def markdup(self, step3_bam):
        sample_path = os.path.join(self.path, 'step3.SamToFastq_bwa_MergeBamAlignment')
        os.makedirs(sample_path, exist_ok=True)
        sorted_bam = os.path.join(sample_path, f'{self.sample_name}.umi.sorted.bam')
        markdup_bam = os.path.join(sample_path, f'{self.sample_name}.umi.markdup.bam')
        script = fr"""
/mnt/share01/tools/bin/samtools sort \
    -T {sample_path} \
    -@ 8 -m 4G \
    -o {sorted_bam} \
    {step3_bam}
/mnt/share01/tools/bin/samtools index {sorted_bam}
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -XX:+UseParallelGC -XX:ParallelGCThreads=8 \
    -Xmx4g \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar MarkDuplicates \
    I={sorted_bam} \
    O={markdup_bam} \
    M={sample_path}/{self.sample_name}.umi.dupStat \
    VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=false \
    MAX_RECORDS_IN_RAM=1750000 TMP_DIR=./alignment
/mnt/share01/tools/bin/samtools index {markdup_bam}
rm {sorted_bam}
        """
        job_name = f'step3_markdup_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step3_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '8', job_script, pre_jobs, sample_path, 'ready')
        self.bamdst(markdup_bam, sample_path, self.target, "step3_bamdst", f'step3_markdup_{self.sample_name}')
        return markdup_bam

    def panel_qc(self, panel_QC_bin, indir):
        sample_path = os.path.join(self.path, 'QC')
        os.makedirs(sample_path, exist_ok=True)
        depth_file = os.path.join(indir, "depthQC", "depth.tsv.gz")
        stats_file = os.path.join(indir, "depthQC", "coverage.report")
        script = fr"""
python {panel_QC_bin} \
    --depth {depth_file} \
    --stats {stats_file} \
    --outfile {sample_path}/{self.sample_name}.qc.xls \
    --gzGCqual {sample_path}/{self.sample_name}.totstat
        """
        job_name = f'step3_panel_qc_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step3_markdup_{self.sample_name}, step0_fastqQC_{self.sample_name}, step3_bamdst_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '8', job_script, pre_jobs, sample_path, 'ready')

    def insert_size(self, insert_size_bin, markdup_bam):
        sample_path = os.path.join(self.path, 'QC')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
python {insert_size_bin} \
    --bam {markdup_bam} \
    --bed {self.target} \
    --InsertSize {sample_path}/{self.sample_name}.insertSize.txt \
    --outfile {sample_path}/{self.sample_name}.insertSize_QC.txt \
    --png {sample_path}/{self.sample_name}.insertSize.png
        """
        job_name = f'step3_insert_size_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step3_markdup_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '8', job_script, pre_jobs, sample_path, 'ready')

    def UMI_QC(self, UMI_QC_bin, step3_bam):
        sample_path = os.path.join(self.path, 'QC')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
python {UMI_QC_bin} \
    --input {step3_bam} \
    --output {sample_path}/{self.sample_name}.UMI_count.txt \
    --bed {self.target} \
    --cpu 8
        """
        job_name = f'step3_UMI_qc_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step3_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '8', job_script, pre_jobs, sample_path, 'ready')

    def step4(self, step3_bam):
        sample_path = os.path.join(self.path, 'step4.GroupReadsByUmi')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java  \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/fgbio/fgbio.jar \
    GroupReadsByUmi \
    --input={step3_bam} \
    --output={sample_path}/{self.sample_name}.umi.group.uBAM \
    --strategy=paired \
    --min-map-q=20 \
    --edits=1
"""
        job_name = f'step4_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step3_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '15', job_script, pre_jobs, sample_path, 'ready')
        step4_bam = os.path.join(sample_path, f'{self.sample_name}.umi.group.uBAM')
        return step4_bam

    def step5(self, step4_bam):
        sample_path = os.path.join(self.path, 'step5.CallDuplexConsensusReads')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java  \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/fgbio/fgbio.jar \
    CallDuplexConsensusReads \
    --min-reads=1 \
    --min-input-base-quality=30 \
    --error-rate-pre-umi=45 \
    --error-rate-post-umi=30 \
    --input={step4_bam} \
    --output={sample_path}/{self.sample_name}.consensus.uBAM
"""
        job_name = f'step5_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step4_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '15', job_script, pre_jobs, sample_path, 'ready')
        step5_bam = os.path.join(sample_path, f'{self.sample_name}.consensus.uBAM')
        return step5_bam

    def pick_step4_bam(self, pick_bam_bin, step4_bam):
        """
        根据设置的梯度挑选bam中的reads，生成一系列的bam
        """
        sample_path = os.path.join(self.path, 'step4.GroupReadsByUmi')
        os.makedirs(sample_path, exist_ok=True)
        sorted_bam = os.path.join(sample_path, f'{self.sample_name}.umi.group.sorted.bam')
        script = fr"""
/mnt/share01/tools/bin/samtools sort \
    -T {sample_path} \
    -@ 8 -m 4G \
    -o {sorted_bam} \
    {step4_bam}
/mnt/share01/tools/bin/samtools index {sorted_bam}
python {pick_bam_bin} \
    --input {sorted_bam} \
    --outDir {sample_path} \
    --samtools /mnt/share01/tools/bin/samtools \
    --cpu 8
        """
        job_name = f'step4_pick_bam_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step4_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '8', job_script, pre_jobs, sample_path, 'ready')
        step4_bam = os.path.join(sample_path, f'{self.sample_name}.umi.group.uBAM')
        return step4_bam

    def rerun_step5(self):
        """
        按比例挑选的bam进行sort，并生成双链一致性序列
        """
        bam_path = os.path.join(self.path, 'step4.GroupReadsByUmi')
        out_path = os.path.join(self.path, 'step5.CallDuplexConsensusReads')
        for i in range(9):
            bam = os.path.join(bam_path, f"pick_{i}.bam")
            sort_bam = os.path.join(bam_path, f"pick_{i}.sort.bam")
            consensus_bam = os.path.join(out_path, f"pick_{i}.consensus.uBAM")
            script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/fgbio/fgbio.jar SortBam \
    --input={bam} \
    --output={sort_bam} \
    -s TemplateCoordinate
rm {bam}

/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/fgbio/fgbio.jar CallDuplexConsensusReads \
    --min-reads=1 \
    --min-input-base-quality=30 \
    --error-rate-pre-umi=45 \
    --error-rate-post-umi=30 \
    --input={sort_bam} \
    --output={consensus_bam}
"""
            job_name = f'rerun_step5_{i}_{self.sample_name}'
            job_script = f'{out_path}/{job_name}.sh'
            pre_jobs = f'step4_pick_bam_{self.sample_name}'
            job_write(self.job_file_open, script, job_name, '8', job_script, pre_jobs, out_path, 'ready')

    def PCT(self, PCT_count_bin, step4_bam, step5_bam):
        """
        计算一致性双链饱和率
        """
        sample_path = os.path.join(self.path, 'QC')
        os.makedirs(sample_path, exist_ok=True)
        r_bam_path = os.path.join(self.path, 'step4.GroupReadsByUmi')
        c_bam_path = os.path.join(self.path, 'step5.CallDuplexConsensusReads')
        script = fr"""
python {PCT_count_bin} \
    --cbam {step5_bam} \
    --rbam {step4_bam} \
    --rbamDir  {r_bam_path} \
    --cbamDir {c_bam_path} \
    --outfile {sample_path}/{self.sample_name}.PCT_status.txt \
    --png {sample_path}/{self.sample_name}.PCT_plot.png \
    --samtools /mnt/share01/tools/bin/samtools \
    --cpu 8
rm {r_bam_path}/pick_*.sort.bam
"""
        job_name = f'step5_PCT_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = ""
        for i in range(9):
            pre_jobs += f'rerun_step5_{i}_{self.sample_name}, '
        pre_jobs += f'step5_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '1', job_script, pre_jobs, sample_path, 'ready')

    def step6(self, step5_bam):
        sample_path = os.path.join(self.path, 'step6.remap')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar \
    SamToFastq \
    I={step5_bam} \
    F=/dev/stdout \
    INTERLEAVE=true \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/bin/bwa mem \
    -p -t 20 {self.ref} /dev/stdin \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar \
    MergeBamAlignment \
    UNMAPPED={step5_bam} \
    ALIGNED=/dev/stdin \
    O={sample_path}/{self.sample_name}.consensus.BAM \
    R={self.ref} \
    SO=coordinate \
    ALIGNER_PROPER_PAIR_FLAGS=true \
    MAX_GAPS=-1 \
    ORIENTATIONS=FR \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true
"""
        job_name = f'step6_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step5_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '10', job_script, pre_jobs, sample_path, 'ready')
        step6_bam = os.path.join(sample_path, f'{self.sample_name}.consensus.BAM')
        return step6_bam

    def step7(self, step6_bam):
        sample_path = os.path.join(self.path, 'step7.FilterConsensusReads')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java  \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/fgbio/fgbio.jar \
    FilterConsensusReads \
    --input={step6_bam} \
    --output={sample_path}/{self.sample_name}.consensus.filter.BAM \
    --ref={self.ref} \
    --min-reads={self.min_reads} \
    --max-read-error-rate=0.05 \
    --max-base-error-rate=0.1 \
    --min-base-quality=50 \
    --max-no-call-fraction=0.05
"""
        job_name = f'step7_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step6_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '10', job_script, pre_jobs, sample_path, 'ready')
        step7_bam = os.path.join(sample_path, f'{self.sample_name}.consensus.filter.BAM')

        self.bamdst(step7_bam, sample_path, self.target, 'step7_bamdst', f'step7_{self.sample_name}')
        return step7_bam

    def step7_libraryQC(self):
        sample_path = os.path.join(self.path, 'QC')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
a=`grep "Average depth" {self.path}/step8.ClipBam/depthQC/coverage.report|grep "Target" |grep "rmdup"| cut -f 2`
b=`grep  "Average depth" {self.path}/step3.SamToFastq_bwa_MergeBamAlignment/depthQC/coverage.report|grep "Target" |grep "rmdup" |cut -f 2`
c=`awk 'BEGIN{{printf "%.2f\n", '$a'/(303*{self.sample_ng})*100}}'`
echo "靶向捕获平均深度(×)    DCS平均深度(×)    文库转化效率(%)" > {sample_path}/{self.sample_name}.libraryQC.txt
echo "$b    $a    $c" >> {sample_path}/{self.sample_name}.libraryQC.txt
"""
        job_name = f'step7_libraryQC'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step7_bamdst_{self.sample_name}, step3_bamdst_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '1', job_script, pre_jobs, sample_path, 'ready')

    def step8(self, step7_bam):
        sample_path = os.path.join(self.path, 'step8.ClipBam')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java  \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/fgbio/fgbio.jar \
    ClipBam \
    --input={step7_bam} \
    --output={sample_path}/{self.sample_name}.consensus.filter.clip.BAM \
    --ref={self.ref} \
    --soft-clip=false \
    --clip-overlapping-reads=true
"""
        job_name = f'step8_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step7_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '5', job_script, pre_jobs, sample_path, 'ready')
        step8_bam = os.path.join(sample_path, f'{self.sample_name}.consensus.filter.clip.BAM')

        self.bamdst(step8_bam, sample_path, self.target, 'step8_bamdst', f'step8_{self.sample_name}')
        return step8_bam

    def step9(self, step8_bam):
        sample_path = os.path.join(self.path, 'step9.vardict')
        os.makedirs(sample_path, exist_ok=True)
        if not self.normal_sample:
            script = fr"""
AF_THR="0.00000001"

/mnt/share01/tools/miniconda/envs/UMI_naangda/share/VarDict-1.8.2/bin/VarDict \
    -G {self.ref} \
    -f $AF_THR \
    -N {self.sample_name} \
    -b {step8_bam} \
    -z -c 1 -S 2 -E 3 -g 4 -r 1 -B 1 -th 4 \
    {self.target} \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/share/VarDict-1.8.2/bin/teststrandbias.R \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/share/VarDict-1.8.2/bin/var2vcf_valid.pl -N {self.sample_name} -E -f $AF_THR \
    | awk '{{if($1 ~ /^#/) print; else if ($4 != $5) print}}' > {sample_path}/tmp.vcf\

/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar \
    SortVcf \
    I={sample_path}/tmp.vcf \
    O={sample_path}/{self.sample_name}.target.vcf \
    SD=/mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/hg19_reference_with_NC.dict
"""
        else:
            script = fr"""
AF_THR="0.00000001"

/mnt/share01/tools/miniconda/envs/UMI_naangda/share/VarDict-1.8.2/bin/VarDict \
    -G {self.ref} \
    -f $AF_THR \
    -N {self.sample_name} \
    -b "{step8_bam}|{self.normal_bam}" \
    -z -c 1 -S 2 -E 3 -g 4 -r 1 -B 1 -th 4 \
    {self.target} \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/share/VarDict-1.8.2/bin/testsomatic.R \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/share/VarDict-1.8.2/bin/var2vcf_paired.pl -N "{self.sample_name}|{self.normal_sample}" -f $AF_THR \
    | awk '{{if($1 ~ /^#/) print; else if ($4 != $5) print}}' > {sample_path}/tmp.somatic.vcf

/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar \
    SortVcf \
    I={sample_path}/tmp.somatic.vcf \
    O={sample_path}/{self.sample_name}.target.somatic.vcf \
    SD=/mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/hg19_reference_with_NC.dict
"""
        job_name = f'step9_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step8_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '5', job_script, pre_jobs, sample_path, 'ready')
        return f"{sample_path}/{self.sample_name}.target.vcf"

    def step10(self, vcf, freq_computebin):
        sample_path = os.path.join(self.path, 'step10.annotation')
        os.makedirs(sample_path, exist_ok=True)
        script = fr"""
perl /mnt/share01/tools/annovar/table_annovar.pl \
    {vcf} \
    /mnt/share01/tools/annovar/humandb/ -buildver hg19 \
    -out {sample_path}/{self.sample_name}.anno \
    -remove -protocol refGene,cytoBand,clinvar_20200316,cosmic68,1000g2015aug_all,snp138,esp6500si_all,ljb26_all -operation g,r,f,f,f,f,f,f -vcfinput
python {freq_computebin} \
    --annot {sample_path}/{self.sample_name}.anno.hg19_multianno.txt \
    --annot_re {sample_path}/{self.sample_name}.hg19_multianno.reform.xls
        """
        job_name = f'step10_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step9_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '1', job_script, pre_jobs, sample_path, 'ready')

    def step_sentieon(self, fq1, fq2):
        sample_path = os.path.join(self.path, 'sentieon')
        os.makedirs(sample_path, exist_ok=True)
        read_structure = self.read_structure
        if self.duplex_umi:
            read_structure = f'-d {read_structure},{read_structure}'
        script = fr"""
WD={sample_path}
TUMOR_SM={self.sample_name}
TUMOR_RGID="rg_$TUMOR_SM"
PL={self.PLATFORM}
TUMOR_FASTQ_1={fq1}
TUMOR_FASTQ_2={fq2}
FASTA={self.ref}
KNOWN_DBSNP={self.dbsnp}
INTERVAL_FILE={self.target}
READ_STRUCTURE="{read_structure}"
MIN_DEPTH={self.min_depth}
NT=24

SENTIEON_INSTALL_DIR={self.bin['sentieon']}
BCFTOOLS_BINARY={self.bin['bcftools']}
export SENTIEON_LICENSE=10.90.1.100:8990

# ******************************************
# 1. Pre-processing of FASTQ containing UMIs
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon umi extract $READ_STRUCTURE $TUMOR_FASTQ_1 $TUMOR_FASTQ_2 | \
  $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -p -C -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" \
  -t $NT -K 10000000 $FASTA - | \
  $SENTIEON_INSTALL_DIR/bin/sentieon umi consensus -o $WD/$TUMOR_SM.umi_consensus.fastq.gz

$SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -p -C -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" \
  -t $NT -K 10000000 $FASTA $WD/$TUMOR_SM.umi_consensus.fastq.gz | \
  $SENTIEON_INSTALL_DIR/bin/sentieon util sort --umi_post_process --sam2bam -i - -o $WD/$TUMOR_SM.umi_consensus.bam

# ******************************************
# 2. Somatic and Structural variant calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i $WD/$TUMOR_SM.umi_consensus.bam \
    --interval $INTERVAL_FILE \
    --algo TNscope \
    --tumor_sample $TUMOR_SM \
    --dbsnp $KNOWN_DBSNP \
    --pcr_indel_model NONE \
    --disable_detector sv \
    --min_tumor_allele_frac 0.0003 \
    --filter_t_alt_frac 0.0003 \
    --max_error_per_read 3 \
    --min_tumor_lod 3.0 \
    --min_base_qual 40 \
    --resample_depth 100000 \
    --assemble_mode 4 \
    $WD/$TUMOR_SM.output_tnscope.pre_filter.vcf.gz

$BCFTOOLS_BINARY annotate -R $INTERVAL_FILE -x "FILTER/triallelic_site" $WD/$TUMOR_SM.output_tnscope.pre_filter.vcf.gz | \
  $BCFTOOLS_BINARY filter -m + -s "low_neighboring_bq" -e "FMT/NBQPS[0] < 47.7" | \
  $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -2.14" | \
  $BCFTOOLS_BINARY filter -m + -s "base_qual_bias" -e "FMT/BaseQRankSumPS[0] < -5.02" | \
  $BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < 11.2" | \
  $BCFTOOLS_BINARY filter -m + -s "low_alt_reads" -e "FMT/ALTHC[0] < 4.5" | \
  $BCFTOOLS_BINARY filter -m + -s "low_depth" -e "SUM(FMT/AD[0:]) < $MIN_DEPTH" | \
  $BCFTOOLS_BINARY sort | \
  $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - $WD/$TUMOR_SM.output_tnscope.filtered.vcf.gz

$BCFTOOLS_BINARY filter -i 'FILTER=="PASS"' $WD/$TUMOR_SM.output_tnscope.filtered.vcf.gz | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - $WD/$TUMOR_SM.output_tnscope.filtered.PASS.vcf.gz
"""
        job_name = f'step_sentieon_{self.sample_name}'
        job_script = f'{sample_path}/{job_name}.sh'
        pre_jobs = f'step0_{self.sample_name}'
        job_write(self.job_file_open, script, job_name, '24', job_script, pre_jobs, sample_path, 'ready')
        sentieon_bam = os.path.join(sample_path, f'{self.sample_name}.umi_consensus.bam')
        self.bamdst(sentieon_bam, sample_path, self.target, 'sentieon_bamdst', f'step_sentieon_{self.sample_name}')

    def main(self):
        fq1, fq2 = self.step0()
        self.fastqc()
        self.fastqQC()
        step1_bam = self.step1(fq1, fq2)
        step2_bam = self.step2(step1_bam)
        step3_bam = self.step3(step2_bam)
        markdup_bam = self.markdup(step3_bam)
        self.panel_qc(self.bin["panel_QC"], os.path.dirname(markdup_bam))
        self.insert_size(self.bin["insert_size"], markdup_bam)
        self.UMI_QC(self.bin["UMI_QC"], step3_bam)
        step4_bam = self.step4(step3_bam)
        self.pick_step4_bam(self.bin["pick_bam"], step4_bam)
        step5_bam = self.step5(step4_bam)
        self.rerun_step5()
        self.PCT(self.bin["PCT_count"], step4_bam, step5_bam)
        step6_bam = self.step6(step5_bam)
        step7_bam = self.step7(step6_bam)
        self.step7_libraryQC()
        step8_bam = self.step8(step7_bam)
        vcf = self.step9(step8_bam)
        self.step10(vcf, self.bin["freq_compute"])
        self.step_sentieon(fq1, fq2)
        self.job_file_open.close()
        os.makedirs(os.path.join(self.path, 'log'), exist_ok=True)
        os.system(
            '/mnt/share01/tools/pipeline/clinical_pipeline/dna_pipeline/script/Job_Dag.py %s' % self.job_file)
        os.system(f'/mnt/share01/tools/pipeline/clinical_pipeline/dna_pipeline/script/sjm_job.py {self.job_file} {self.path}/log')


def job_write(job_file_open, shell_info, job_name, ncores, job_script, pre_jobs, work_dir, status):
    sge_job_script = fr"""
[{job_name}]
job_id = 0
ncores = {ncores}
status = {status}
queue = all.q
job_script = {job_script}
work_dir = {work_dir}
pre_jobs = {pre_jobs}

"""
    info = """set -euo pipefail\necho ==== start at `date "+%F  %H:%M:%S"` ====\n"""
    info += shell_info
    info += """\necho ==== end at `date "+%F  %H:%M:%S"` ====\n"""
    open(job_script, 'w').write(info)
    job_file_open.write(sge_job_script)


if __name__ == '__main__':
    umi = Run()
    umi.main()

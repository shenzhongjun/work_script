#!/mnt/share01/tools/bin/python
# -*- coding: UTF-8 -*-

__author__ = "Fei Xue"
__email__ = "xuefei@zhuanhuayixue.org"
__version__ = "1a"
__date__ = "8/19/21 17:08 "
__copyright__ = "Copyright (C) 2021 All rights reserved"

import os
import time
import argparse
localtime = time.strftime("%Y_%m_%d_%H%M%S")

"""
python pipline.py \
	--sample ng10 \
	--path /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/other \
	--fq1 /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/NL190929-1C.R1.fastq.gz \
	--fq2 /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/NL190929-1C.R2.fastq.gz \
	--ref /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \
	--PLATFORM_UNIT HiseqX10 \
	--PLATFORM Illumina \
	--read_structure 3M2S+T \
	--min_reads "6 3 3" \
	--target /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/snp.probe.bed
"""

def my_parser():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('--sample', help="", required=True)
	parser.add_argument('--path', help="", required=True)
	parser.add_argument('--fq1', help="", required=True)
	parser.add_argument('--fq2', help="", required=True)
	parser.add_argument('--ref', help="", default='/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta', required=True)
	parser.add_argument('--PLATFORM_UNIT', help="", choices=['HiseqX10', 'T7'], required=True)
	parser.add_argument('--PLATFORM', help="", choices=['Illumina', 'MGISEQ'], required=True)
	parser.add_argument('--read_structure', help="Illumina, BGI", choices=['3M2S+T', '3M4S+T'], required=True)
	parser.add_argument('--min_reads', help="I", choices=['3 3 0', '6 3 3'], default='6 3 3', required=True)
	parser.add_argument('--target',
						help="/mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/snp.probe.bed",
						required=True)
	return parser.parse_args()


class Run:

	def __init__(self):
		self.argv = my_parser()
		self.sample_name = self.argv.sample
		self.path = self.argv.path
		self.target = self.argv.target
		self.fq1 = self.argv.fq1
		self.fq2 = self.argv.fq2
		self.PLATFORM_UNIT = self.argv.PLATFORM_UNIT
		self.PLATFORM = self.argv.PLATFORM
		self.read_structure = self.argv.read_structure
		self.ref = self.argv.ref
		self.min_reads = self.argv.min_reads
		self.job_file = os.path.join(self.path, localtime + ".job")
		self.job_file_open = open(self.job_file, 'w')

	def step0(self):
		sample_path = os.path.join(self.path, 'step0.fastp')
		os.makedirs(sample_path, exist_ok=True)
		script = fr"""
/mnt/share01/tools/bin/fastp \
	-i {self.fq1} \
	-I {self.fq2} \
	-o {sample_path}/{self.sample_name}.chean.R1.fastq.gz \
	-O {sample_path}/{self.sample_name}.chean.R2.fastq.gz \
	--thread 16 \
	--json {sample_path}/{self.sample_name}.fastp.json \
	--html {sample_path}/{self.sample_name}.fastp.html\
	
python /mnt/share02/xuefei/clinical_WGBS/script/QC/fastp_plot.py \
        --fastp-json {sample_path}/{self.sample_name}.fastp.json \
        --outdir {sample_path} \
        --samplename {self.sample_name}
"""
		job_name = f'step0_{self.sample_name}'
		job_script = f'{sample_path}/{job_name}.sh'
		pre_jobs = 'None'
		job_write(self.job_file_open, script, job_name, '16', job_script, pre_jobs, sample_path, 'ready')

	def step1(self):
		sample_path = os.path.join(self.path, 'step1.FastqToSam')
		os.makedirs(sample_path, exist_ok=True)
		script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar \
    FastqToSam \
    FASTQ={self.fq1} \
    FASTQ2={self.fq2} \
    OUTPUT={sample_path}/{self.sample_name}.uBAM \
    READ_GROUP_NAME={self.sample_name} \
    SAMPLE_NAME={self.sample_name} \
    LIBRARY_NAME={self.sample_name} \
    PLATFORM_UNIT={self.PLATFORM_UNIT} \
    PLATFORM={self.PLATFORM}
"""
		job_name = f'step1_{self.sample_name}'
		job_script = f'{sample_path}/{job_name}.sh'
		pre_jobs = 'None'
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
"""
		job_name = f'step3_{self.sample_name}'
		job_script = f'{sample_path}/{job_name}.sh'
		pre_jobs = f'step2_{self.sample_name}'
		job_write(self.job_file_open, script, job_name, '15', job_script, pre_jobs, sample_path, 'ready')
		step3_bam = os.path.join(sample_path, f'{self.sample_name}.umi.merged.uBAM')

		script = fr"""
/mnt/share01/tools/bamdst/bamdst-master/bamdst \
	-p {self.target} \
	-o {sample_path} \
	{step3_bam}
"""
		job_name = f'step3_bamdst_{self.sample_name}'
		job_script = f'{sample_path}/{job_name}.sh'
		pre_jobs = f'step3_{self.sample_name}'
		job_write(self.job_file_open, script, job_name, '10', job_script, pre_jobs, sample_path, 'ready')

		return step3_bam

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

		script = fr"""
/mnt/share01/tools/bamdst/bamdst-master/bamdst \
	-p {self.target} \
	-o {sample_path} \
	{step6_bam}
"""
		job_name = f'step6_bamdst_{self.sample_name}'
		job_script = f'{sample_path}/{job_name}.sh'
		pre_jobs = f'step6_{self.sample_name}'
		job_write(self.job_file_open, script, job_name, '10', job_script, pre_jobs, sample_path, 'ready')

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

		script = fr"""
/mnt/share01/tools/bamdst/bamdst-master/bamdst \
	-p {self.target} \
	-o {sample_path} \
	{step7_bam}
"""
		job_name = f'step7_bamdst_{self.sample_name}'
		job_script = f'{sample_path}/{job_name}.sh'
		pre_jobs = f'step7_{self.sample_name}'
		job_write(self.job_file_open, script, job_name, '10', job_script, pre_jobs, sample_path, 'ready')

		return step7_bam

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

		script = fr"""
/mnt/share01/tools/bamdst/bamdst-master/bamdst \
	-p {self.target} \
	-o {sample_path} \
	{step8_bam}
"""
		job_name = f'step8_bamdst_{self.sample_name}'
		job_script = f'{sample_path}/{job_name}.sh'
		pre_jobs = f'step8_{self.sample_name}'
		job_write(self.job_file_open, script, job_name, '10', job_script, pre_jobs, sample_path, 'ready')

		return step8_bam

	def step9(self, step8_bam):
		sample_path = os.path.join(self.path, 'step9.vardict')
		os.makedirs(sample_path, exist_ok=True)
		script = fr"""
AF_THR="0.00000001"

/mnt/share01/tools/miniconda/envs/UMI_naangda/share/VarDictJava-1.5.1/build/install/VarDict/bin/VarDict\
    -G {self.ref} \
    -f $AF_THR \
    -N {self.sample_name} \
    -b {step8_bam} \
    -z -c 1 -S 2 -E 3 -g 4 -r 1 -B 1 -th 4 \
    {self.target} \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/share/VarDictJava-1.5.1/VarDict/teststrandbias.R \
    | /mnt/share01/tools/miniconda/envs/UMI_naangda/share/VarDictJava-1.5.1/VarDict/var2vcf_valid.pl -N {self.sample_name} -E -f $AF_THR \
    | awk '{{if($1 ~ /^#/) print; else if ($4 != $5) print}}' > {sample_path}/tmp.vcf\
    
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/picard-2.19.2-0/picard.jar \
    SortVcf \
    I={sample_path}/tmp.vcf \
    O={sample_path}/{self.sample_name}.target.vcf \
    SD=/mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/hg19_reference_with_NC.dict
"""
		job_name = f'step9_{self.sample_name}'
		job_script = f'{sample_path}/{job_name}.sh'
		pre_jobs = f'step8_{self.sample_name}'
		job_write(self.job_file_open, script, job_name, '5', job_script, pre_jobs, sample_path, 'ready')

	def main(self):
		self.step0()
		step1_bam = self.step1()
		step2_bam = self.step2(step1_bam)
		step3_bam = self.step3(step2_bam)
		step4_bam = self.step4(step3_bam)
		step5_bam = self.step5(step4_bam)
		step6_bam = self.step6(step5_bam)
		step7_bam = self.step7(step6_bam)
		step8_bam = self.step8(step7_bam)
		self.step9(step8_bam)

		self.job_file_open.close()
		os.system(
			'/mnt/share01/tools/pipeline/clinical_pipeline/dna_pipeline/script/Job_Dag.py %s' % self.job_file)
		os.system('/mnt/share01/tools/pipeline/clinical_pipeline/dna_pipeline/script/sjm_job.py %s %s/log' % (
			self.job_file, self.path))


def job_write(job_file_open, shell_info, job_name, ncores, job_script, pre_jobs, work_dir, status):
	sge_job_script = '''
	[{job_name}]
	job_id = 0
	ncores = {ncores}
	status = {status}
	queue = all.q
	job_script = {job_script}
	work_dir = {work_dir}
	pre_jobs = {pre_jobs}

	'''.replace('\t', '').format(job_name=job_name,
								 ncores=ncores,
								 job_script=job_script,
								 pre_jobs=pre_jobs,
								 work_dir=work_dir,
								 status=status)

	info = """set -e\necho ==== start at `date "+%F  %H:%M:%S"` ====\n"""
	info += shell_info
	info += """\necho ==== end__ at `date "+%F  %H:%M:%S"` ====\n"""
	open(job_script, 'w').write(info)
	job_file_open.write(sge_job_script)


if __name__ == '__main__':
	object = Run()
	object.main()





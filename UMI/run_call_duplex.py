# /usr/bin/env python
# conding:utf-8
#生成call duplex步骤的shell脚本
import os
import argparse

class call_duplex:
    def __init__(self, sample_name, path, bam):
        self.sample_name = sample_name
        self.path = path
        self.bam = bam
        self.sample_path = os.path.join(self.path, 'step5.CallDuplexConsensusReads')
        os.makedirs(self.sample_path, exist_ok=True)

    def run_sort(self):
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/fgbio/fgbio.jar SortBam \
    --input={self.bam} \
    --output={self.sample_path}/{self.sample_name}.sort.bam \
    -s TemplateCoordinate
        """
        return script
    def run_call_duplex(self):
        script = fr"""
/mnt/share01/tools/miniconda/envs/UMI_naangda/bin/java  \
    -Xmx20G \
    -jar /mnt/share01/tools/miniconda/envs/UMI_naangda/share/fgbio/fgbio.jar \
    CallDuplexConsensusReads \
    --min-reads=1 \
    --min-input-base-quality=30 \
    --error-rate-pre-umi=45 \
    --error-rate-post-umi=30 \
    --input={self.sample_path}/{self.sample_name}.sort.bam \
    --output={self.sample_path}/{self.sample_name}.consensus.uBAM
        """
        return script

    def deploy(self):
        out_shell = os.path.join(self.sample_path, "run_{sample_name}.sh".format(sample_name=self.sample_name))
        out = open(out_shell, "w")
        sort_script = self.run_sort()
        call_duplex_script = self.run_call_duplex()
        out.write("set -e")
        out.write(sort_script + call_duplex_script)
        out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input', help="input bam file", required=True)
    parser.add_argument('--sample', help="sample name", required=True)
    parser.add_argument('--outdir', help="output dir", required=True)
    argv = parser.parse_args()
    p = call_duplex(argv.sample, argv.outdir, argv.input)
    p.deploy()
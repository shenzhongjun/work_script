"""
产生分析流程
"""
import os
import argparse


def run(fq1, fq2, outfile, sample, outdir):
    script = fr"""
python /mnt/share02/wangwp/scripts/project/UMI/pipline.py \
    --sample {sample} \
    --path {outdir}/{sample} \
    --fq1 {fq1} \
    --fq2 {fq2} \
    --ref /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \
    --PLATFORM_UNIT T7 \
    --PLATFORM MGISEQ \
    --read_structure 3M4S+T \
    --min_reads "3 3 0" \
    --config /mnt/share02/wangwp/scripts/project/UMI/UMI_pipline.yaml
    """
    out = open(outfile, "w")
    out.write(script)

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--fqlist', help="", required=True)
    parser.add_argument('--outdir', help="", required=True)
    argv = parser.parse_args()
    for line in open(argv.fqlist):
        sample, fq1, fq2 = line.strip().split("\t")
        outshell = os.path.join(argv.outdir, f'run_{sample}.sh')
        run(fq1, fq2, outshell, sample, argv.outdir)

if __name__ == '__main__':
    main()


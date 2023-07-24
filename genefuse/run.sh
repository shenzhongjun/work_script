path=`pwd`

/mnt/share02/xuefei/dna/genefusion/packaging_test/make_genefuse.py \
    --name ACB472 \
    --out_path ${path} \
    --ref /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \
    --refflat /mnt/share02/xuefei/dna/genefusion/db/hg19/refFlat.txt \
    --read1 ACB472_1.fastq.gz \
    --read2 ACB472_2.fastq.gz \
    --fusion_csv /mnt/share02/xuefei/dna/genefusion/db/hg19/cancer.hg19.csv \
    --unique 6 \
    --thread 20 \
    --deletion 50 \
    --gene_list ${path}/list

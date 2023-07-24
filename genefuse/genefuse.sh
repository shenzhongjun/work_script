
# 第一步：配置需要验证的基因和转录本号\
# 		 方法：在 /mnt/share02/xuefei/dna/genefusion/db/hg19/refFlat.txt 中 grep你需要的基因和转录本号 写入 gene_list.txt 文件中
#        gene_list.txt 文件 第一列是基因，第二列是转录本号

# 第二步：gene_list.txt 作为输入，执行下面代码，生成 my_fusion_AEB017.csv

fq1=/mnt/share05/clinical_project/projects/solid_tumor/DDN21022850_PN2110220005/DNA_NCPLt_NCPLt/clean_data/AEI144/AEI144.clean_1.fq.gz
fq2=/mnt/share05/clinical_project/projects/solid_tumor/DDN21022850_PN2110220005/DNA_NCPLt_NCPLt/clean_data/AEI144/AEI144.clean_2.fq.gz
out=`pwd`
#/mnt/share01/tools/miniconda/envs/cnvkit/bin/python \
#    /mnt/share02/xuefei/dna/genefusion/test/genefuse/scripts/make_fusion_genes.py \
#    $out/gene_list.txt \
#    -r /mnt/share02/xuefei/dna/genefusion/db/hg19/refFlat.txt \
#    -o $out/my_fusion.csv && \
    

# 第三步：my_fusion.csv 作为输入，执行以下代码

/mnt/share01/tools/miniconda/envs/genefusion/bin/genefuse \
    --read1 $fq1 \
    --read2 $fq2 \
    --fusion $out/my_fusion.csv \
    --ref /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta --thread 20 --deletion 50 --unique 2 \
    --html $out/genefuse.html \
    --json $out/genefuse.json \
    --output_deletions $out/long_deletions.txt \
    --output_untranslated_fusions $out/untranslated_fusions.txt	

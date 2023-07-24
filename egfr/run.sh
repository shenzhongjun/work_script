out=`pwd`
prefix=DDN21019340

/mnt/share01/tools/bin/pindel \
    -f /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \
    -c chr7:55086970-55273310 \
    -i pindel.conf \
    -o $out/$prefix \

/mnt/share01/tools/miniconda/envs/pindel/bin/pindel2vcf \
    -p $out/${prefix}_D \
    -r /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \
    -R hg19-UCSC -d 200902 \
    -v $out/$prefix.vcf

perl /mnt/share01/tools/annovar/table_annovar.pl \
	$prefix.vcf \
	/mnt/share01/tools/annovar/humandb/ -buildver hg19 \
	-out $prefix.anno -remove \
	-remove -protocol refGene,cytoBand,clinvar_20200316 -operation g,r,f -vcfinput

python /mnt/share02/lixx/04.project_test/06.EGFR19del/muts_table.py \
	--original_anno $out/$prefix.anno.hg19_multianno.txt


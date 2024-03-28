# normal
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep \
	--input_file AIJ481.merged.vcf.gz \
	--output_file AIJ481.merged.anno_vep.final.vcf \
	--dir_cache /mnt/share02/zhouyj/database/annotation \
	--dir_plugins /home/zhouyj/anaconda3/envs/vep/Plugins \
	--fasta /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta \
	--format vcf \
	--vcf \
	--cache \
	--offline \
	--merged \
	--buffer_size 10000 \
	--fork 16 \
	--force_overwrite \
	--variant_class \
	--sift b \
	--polyphen b \
	--numbers \
	--gene_phenotype \
	--hgvs \
	--protein \
	--symbol \
	--canonical \
	--biotype \
	--domains \
	--transcript_version \
	--ccds \
	--check_existing \
	--af \
	--af_1kg \
	--af_gnomade \
	--pubmed

# 使用config文件
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep \
	--config /home/zhouyj/anaconda3/envs/vep/basic.config.ini \
	--input_file AIJ481.merged.vcf.gz \
	--output_file AIJ481.merged.anno_vep.final.vcf \
	--fork 16

# 使用config文件+Plugin+Custom
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep \
	--config /home/zhouyj/anaconda3/envs/vep/plugin_custom.config.ini \
	--input_file AIJ481.merged.vcf.gz \
	--output_file AIJ481.merged.anno_vep.final.vcf \
	--fork 8

# other
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep --input_file AIJ481.merged.vcf.gz --output_file AIJ481.merged.anno_vep.nsfp.txt --tab --fork 16 --dir_cache /mnt/share02/zhouyj/database/annotation --dir_plugins /home/zhouyj/anaconda3/envs/vep/Plugins --fasta /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta --format vcf --cache --offline --merged --force_overwrite --show_ref_allele --plugin dbNSFP,/mnt/share02/zhouyj/database/annotation/dbNSFP/dbNSFP4.3a_grch37.gz,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,REVEL_score,MutationTaster_score,MutationTaster_pred,CADD_phred,FATHMM_score,FATHMM_pred --sift b --polyphen b --custom /mnt/share02/zhouyj/database/annotation/COSMICv97/hg19/CosmicCodingMuts.normal.vcf.gz,COSMIC_num,vcf,exact,,ID,LEGACY_ID,CDS,AA,CNT --pick

# other-config
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep \
	--config /mnt/share02/zhouyj/database/annotation/plugin_custom.config.ini \
	--input_file AIJ481.merged.vcf.gz \
	--output_file AIJ481.merged.anno_vep.test.txt \
	--tab \
	--fork 64 \
	--plugin dbNSFP,/mnt/share02/zhouyj/database/annotation/dbNSFP/dbNSFP4.3a_grch37.gz,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,REVEL_score,MutationTaster_score,MutationTaster_pred,CADD_phred,FATHMM_score,FATHMM_pred \
	--plugin SpliceAI,snv=/mnt/share02/zhouyj/database/annotation/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/mnt/share02/zhouyj/database/annotation/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
	--plugin dbscSNV,/mnt/share02/zhouyj/database/annotation/dbscSNV/dbscSNV1.1_GRCh37.txt.gz \
	--custom /mnt/share02/zhouyj/database/annotation/COSMICv97/hg19/CosmicCodingMuts.normal.vcf.gz,COSMIC,vcf,exact,,LEGACY_ID,CNT

# other-config-vcf
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep \
	--config /mnt/share02/zhouyj/database/annotation/plugin_custom.config.ini \
	--input_file AIJ481.merged.vcf.gz \
	--output_file AIJ481.merged.anno_vep.test.vcf \
	--vcf \
	--fork 64 \
	--plugin dbNSFP,/mnt/share02/zhouyj/database/annotation/dbNSFP/dbNSFP4.3a_grch37.gz,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,REVEL_score,MutationTaster_score,MutationTaster_pred,CADD_phred,FATHMM_score,FATHMM_pred \
	--plugin SpliceAI,snv=/mnt/share02/zhouyj/database/annotation/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/mnt/share02/zhouyj/database/annotation/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
	--plugin dbscSNV,/mnt/share02/zhouyj/database/annotation/dbscSNV/dbscSNV1.1_GRCh37.txt.gz \
	--custom /mnt/share02/zhouyj/database/annotation/COSMICv97/hg19/CosmicCodingMuts.normal.vcf.gz,COSMIC,vcf,exact,,LEGACY_ID,CNT

# 结果整理
sed -n '1,50000p' AIJ481.merged.anno_vep.test.txt > AIJ481.merged.anno_vep.test.sub.txt
python handle_vep_result.py --anno AIJ481.merged.anno_vep.test.sub.txt --vcf AIJ481.merged.vcf.gz --out test


#!/bin/sh

set -auo pipefail
echo ==== start at $(date "+%F  %H:%M:%S") ====

db=/mnt/share02/zhouyj/database/annotation
dbnsfp_str="HGVSc_VEP,CADD_phred,FATHMM_pred,SIFT_pred,Polyphen2_HDIV_pred,PROVEAN_pred,MutationTaster_pred,M-CAP_pred,REVEL_score"
dbnsfp_str=$dbnsfp_str",clinvar_id,clinvar_clnsig,clinvar_trait,Interpro_domain"

fields_str="Uploaded_variation,Location,REF_ALLELE,Allele,SYMBOL,Gene,Feature,HGVSc,HGVSp,Consequence,IMPACT,EXON,STRAND,Existing_variation,PUBMED,"
fields_str=$fields_str$dbnsfp_str",SpliceAI_pred,ada_score,rf_score,COSMIC,COSMIC_CNT"
fields_str=$fields_str",AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_EAS_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,WBBC_AF,WBBC_North_AF,WBBC_Central_AF,WBBC_South_AF,WBBC_Lingnan_AF"

# VEP注释。
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep \
	--config $db/plugin_custom.config.ini \
	--input_file AIJ481.merged.vcf.gz \
	--output_file AIJ481.anno_vep.txt \
	--distance 10 \
	--tab \
	--fields $fields_str \
	--fork 16 \
	--plugin dbNSFP,$db/dbNSFP/dbNSFP4.3a_grch37.gz,$dbnsfp_str \
	--plugin SpliceAI,snv=$db/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=$db/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
	--plugin dbscSNV,$db/dbscSNV/dbscSNV1.1_GRCh37.txt.gz \
	--custom $db/COSMICv97/hg19/CosmicCodingMuts.normal.vcf.gz,COSMIC,vcf,exact,,CNT \
	--custom $db/WBBC/hg19/WBBC.GRCh37.vcf.gz,WBBC,vcf,exact,,AF,North_AF,Central_AF,South_AF,Lingnan_AF

# 结果整理
#sed -n '1,10000p' AIJ481.merged.anno_vep.test.txt > AIJ481.merged.anno_vep.test.sub.txt
python handle_vep_result.py --anno AIJ481.anno_vep.txt --vcf AIJ481.merged.vcf.gz --out AIJ481.anno_vep.merged.txt

echo ==== end__ at $(date "+%F  %H:%M:%S") ====

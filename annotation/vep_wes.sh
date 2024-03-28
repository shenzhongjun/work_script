#!/bin/sh

set -auo pipefail
echo ==== start at $(date "+%F  %H:%M:%S") ====

db=/mnt/share02/zhouyj/database/annotation
dbnsfp_str="HGVSc_VEP,CADD_phred_hg19,FATHMM_pred,SIFT_pred,Polyphen2_HDIV_pred,PROVEAN_pred,MutationTaster_pred,VEST4_score,REVEL_score"
dbnsfp_str=$dbnsfp_str",clinvar_id,clinvar_clnsig,clinvar_trait,Interpro_domain"
dbnsfp_str=$dbnsfp_str",gnomAD_exomes_controls_AF,gnomAD_exomes_controls_AFR_AF,gnomAD_exomes_controls_AMR_AF,gnomAD_exomes_controls_EAS_AF,gnomAD_exomes_controls_NFE_AF,gnomAD_exomes_controls_SAS_AF"

fields_str="Uploaded_variation,Location,REF_ALLELE,Allele,SYMBOL,Gene,Feature,HGVSg,HGVSc,HGVSp,Consequence,IMPACT,EXON,STRAND,Existing_variation,PUBMED,"
fields_str=$fields_str$dbnsfp_str",AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,WBBC_AF,WBBC_North_AF,WBBC_Central_AF,WBBC_South_AF,WBBC_Lingnan_AF"
fields_str=$fields_str",SpliceAI_pred,ada_score,rf_score,COSMIC,COSMIC_LEGACY_ID,COSMIC_CNT"

# VEP注释
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep \
	--config $db/plugin_custom.config.ini \
	--input_file /mnt/share05/clinical_project/projects/blood_tumor/test/xuefei_test/v2.0test/v2.0/snvindel/merge/AIG138/AIG138.merged.vcf.gz \
	--output_file AIG138.anno_vep.txt \
	--distance 10 \
	--tab \
	--fields $fields_str \
	--fork 64 \
	--plugin dbNSFP,$db/dbNSFP/dbNSFP4.3a_grch37.gz,$dbnsfp_str \
	--plugin SpliceAI,snv=$db/spliceAI/spliceai_scores.raw.snv.hg19.vcf.gz,indel=$db/spliceAI/spliceai_scores.raw.indel.hg19.vcf.gz,cutoff=0.5 \
	--plugin dbscSNV,$db/dbscSNV/dbscSNV1.1_GRCh37.txt.gz \
	--custom $db/COSMICv97/hg19/CosmicCodingMuts.normal.vcf.gz,COSMIC,vcf,exact,,LEGACY_ID,CNT \
	--custom $db/WBBC/hg19/WBBC.GRCh37.vcf.gz,WBBC,vcf,exact,,AF,North_AF,Central_AF,South_AF,Lingnan_AF

# 结果整理
#sed -n '1,10000p' AIJ481.merged.anno_vep.other.txt > AIJ481.merged.anno_vep.other.sub.txt
python handle_vep_result.v0.2.py \
	--anno AIG138.anno_vep.txt \
	--vcf /mnt/share05/clinical_project/projects/blood_tumor/test/xuefei_test/v2.0test/v2.0/snvindel/merge/AIG138/AIG138.merged.vcf.gz \
	--out AIG138.anno_vep.merged.txt

echo ==== end__ at $(date "+%F  %H:%M:%S") ====

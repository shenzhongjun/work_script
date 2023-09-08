#!/bin/sh

set -auo pipefail
echo ==== start at $(date "+%F  %H:%M:%S") ====

sample=AJL106		# Panel
#sample=AIG138		# WES
db=/mnt/share02/zhouyj/database/annotation
dbnsfp_str="HGVSc_VEP,CADD_phred_hg19,FATHMM_score,SIFT_score,Polyphen2_HDIV_score,MutationTaster_pred,VEST4_score,REVEL_score,Interpro_domain"
# dbNSFP用的gnomAD太旧，有的突变未收录导致报出为“致病”，改为使用VEP自带gnomAD
# dbnsfp_str=$dbnsfp_str",Interpro_domain,gnomAD_exomes_controls_AF,gnomAD_exomes_controls_AFR_AF,gnomAD_exomes_controls_AMR_AF,gnomAD_exomes_controls_EAS_AF,gnomAD_exomes_controls_NFE_AF,gnomAD_exomes_controls_SAS_AF"

fields_str="Uploaded_variation,Location,REF_ALLELE,Allele,SYMBOL,Gene,Feature,HGVSg,HGVSc,HGVSp,Consequence,IMPACT,EXON,STRAND,Existing_variation,PUBMED,"
fields_str=$fields_str"gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_EAS_AF,gnomADe_NFE_AF,gnomADe_SAS_AF,"
fields_str=$fields_str$dbnsfp_str",AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,WBBC_AF,WBBC_North_AF,WBBC_Central_AF,WBBC_South_AF,WBBC_Lingnan_AF"
fields_str=$fields_str",SpliceAI_pred,ada_score,rf_score,COSMIC,COSMIC_LEGACY_ID,COSMIC_CNT,ClinVar,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_CLNVI,ClinVar_CLNREVSTAT,UniProt"

# VEP注释
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep \
	--config $db/plugin_custom.config.ini \
	--input_file ${sample}.merged.vcf.gz \
	--output_file ${sample}.anno_vep.txt \
	--distance 10 \
	--tab \
	--fields $fields_str \
	--fork 32 \
	--plugin dbNSFP,$db/dbNSFP/dbNSFP4.3a_grch37.gz,$dbnsfp_str \
	--plugin SpliceAI,snv=$db/spliceAI/spliceai_scores.raw.snv.hg19.vcf.gz,indel=$db/spliceAI/spliceai_scores.raw.indel.hg19.vcf.gz,cutoff=0.5 \
	--plugin dbscSNV,$db/dbscSNV/dbscSNV1.1_GRCh37.txt.gz \
	--custom $db/clinvar/clinvar_20230527.vcf.gz,ClinVar,vcf,exact,,CLNSIG,CLNDN,CLNVI,CLNREVSTAT \
	--custom $db/COSMICv97/hg19/CosmicCodingMuts.normal.vcf.gz,COSMIC,vcf,exact,,LEGACY_ID,CNT \
	--custom $db/WBBC/hg19/WBBC.GRCh37.vcf.gz,WBBC,vcf,exact,,AF,North_AF,Central_AF,South_AF,Lingnan_AF \
	--custom $db/custom/UniProt/uniprot_for_OM1.hg19.20230706.bed.gz,UniProt,bed,overlap,

# 结果整理
sed -n '1,10000p' ${sample}.anno_vep.txt > ${sample}.anno_vep.tmp.txt
#cp ${sample}.anno_vep.txt ${sample}.anno_vep.tmp.txt

/home/zhouyj/anaconda3/envs/python3.9/bin/python make_mutation_table_test.py \
	--anno ${sample}.anno_vep.tmp.txt \
	--vcf ${sample}.merged.vcf.gz \
	--out ${sample}.anno_vep.merged

echo ==== end__ at $(date "+%F  %H:%M:%S") ====

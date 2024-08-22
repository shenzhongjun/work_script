#!/bin/sh

set -auo pipefail
echo ==== start at $(date "+%F  %H:%M:%S") ====

db=/mnt/share02/zhouyj/database/annotation
dbnsfp_str="HGVSc_VEP,CADD_phred_hg19,FATHMM_score,SIFT_score,Polyphen2_HDIV_score,MutationTaster_pred,VEST4_score,REVEL_score,Interpro_domain"

fields_str="Uploaded_variation,Location,REF_ALLELE,Allele,SYMBOL,Gene,Feature,HGVSg,HGVSc,HGVSp,Consequence,IMPACT,EXON,INTRON,STRAND,Existing_variation,PUBMED,"
fields_str=$fields_str"gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_EAS_AF,gnomADe_NFE_AF,gnomADe_SAS_AF,"
fields_str=$fields_str$dbnsfp_str",AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,WBBC_AF,WBBC_North_AF,WBBC_Central_AF,WBBC_South_AF,WBBC_Lingnan_AF"
fields_str=$fields_str",SpliceAI_pred,ada_score,rf_score,COSMIC,COSMIC_LEGACY_ID,COSMIC_CNT,ClinVar,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_CLNVI,ClinVar_CLNREVSTAT,UniProt"

# VEP注释
/home/zhouyj/anaconda3/envs/vep/bin/perl /home/zhouyj/anaconda3/envs/vep/bin/vep \
	--config $db/plugin_custom.config.ini \
	--input_file /mnt/share05/clinical_project/projects/blood_tumor/test/xuefei_test/v2.0test/v2.0/snvindel/merge/AIG138/AIG138.merged.vcf.gz \
	--output_file AIG138.anno_vep.txt \
	--distance 10 \
	--tab \
	--fields $fields_str \
	--fork 16 \
	--plugin dbNSFP,$db/dbNSFP/dbNSFP4.3a_grch37.gz,$dbnsfp_str \
	--plugin SpliceAI,snv=$db/spliceAI/spliceai_scores.raw.snv.hg19.vcf.gz,indel=$db/spliceAI/spliceai_scores.raw.indel.hg19.vcf.gz,cutoff=0.5 \
	--plugin dbscSNV,$db/dbscSNV/dbscSNV1.1_GRCh37.txt.gz \
	--custom $db/COSMICv97/hg19/CosmicCodingMuts.normal.vcf.gz,COSMIC,vcf,exact,,LEGACY_ID,CNT \
	--custom $db/clinvar/clinvar_20230527.vcf.gz,ClinVar,vcf,exact,,CLNSIG,CLNDN,CLNVI,CLNREVSTAT \
	--custom $db/WBBC/hg19/WBBC.GRCh37.vcf.gz,WBBC,vcf,exact,,AF,North_AF,Central_AF,South_AF,Lingnan_AF \
	--custom $db/custom/UniProt/uniprot_for_OM1.hg19.20230706.bed.gz,UniProt,bed,overlap,

# 结果整理
#sed -n '1,10000p' AIJ481.merged.anno_vep.other.txt > AIJ481.merged.anno_vep.other.sub.txt
/mnt/share02/zhouyj/software/miniconda3/envs/python3.9/bin/python /mnt/share02/zhouyj/script/annotation/make_mutation_table.py \
	--anno AIG138.anno_vep.txt \
	--vcf /mnt/share05/clinical_project/projects/blood_tumor/test/xuefei_test/v2.0test/v2.0/snvindel/merge/AIG138/AIG138.merged.vcf.gz \
	--out AIG138.anno_vep.merged

echo ==== end__ at $(date "+%F  %H:%M:%S") ====

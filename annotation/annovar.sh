perl /mnt/share02/zhouyj/software/annovar/table_annovar.pl \
	test-3.vcf.gz \
	/mnt/share02/zhouyj/software/annovar/humandb/ -buildver hg19 \
	-out test-3.annovar \
	-vcfinput -remove --thread 10 \
	-protocol refGene,cytoBand,avsnp150,snp138,clinvar_20220320,cosmic70,1000g2015aug_all,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur,exac03,gnomad211_genome,gnomad211_exome,esp6500siv2_all,cadd13gt20,gerp++gt2,ljb26_all,intervar_20180118,wgRna,targetScanS,gwasCatalog,rmsk \
	-operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r

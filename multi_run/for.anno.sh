out=$PWD
for i  in `cat list`
do
	cd $out/$i/Sentieon/$i
	echo -e "
perl /mnt/share02/zhangxc/software/annovar/table_annovar.pl \\
		$out/$i/Sentieon/$i/$i.output_tnscope.filtered.PASS.vcf.gz \\
		/mnt/share01/tools/annovar/humandb/ -buildver hg19 \\
		-out $out/$i/Sentieon/$i/$i.filtered.PASS.anno \\
		-remove -protocol refGene,cytoBand,clinvar_20200316,cosmic68,1000g2015aug_all,snp138,esp6500si_all,ljb26_all \\
		-operation g,r,f,f,f,f,f,f -vcfinput
" > anno.sh
qsub -V -cwd -l p=10 anno.sh
cd -
done

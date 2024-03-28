less /mnt/share05/data/product_raw_data/rawdata/script/nas_data/外送测序信息汇总-肿瘤全部信息.txt|grep IDT|grep -E 'NBCW|Ncet|NPCW|Ncec'|cut -f 2,3,4,5,7,15|grep DZ|tail -300 > IDTsamples.txt
python get_bams.py -i step0.get_sample_bams/list -o sample_bam_list.txt
python step1_call_snp_indel.py
cd step1.call_snp_indel/shell
for i in `ls`;do qsub $i;done

# step2.create_genomics_db
mkdir -p step2.create_genomics_db/tmp && cd step2.create_genomics_db
cat /mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/IDT.bed \
	/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NPC980.bed \
	/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NLP.bed \
	/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NVT.bed \
	| sort -k1,1 -k2,2n > merged_sorted.bed
bedtools merge -i merged_sorted.bed  > final.bed
for i in `ls ../step1.call_snp_indel/*.vcf.gz`;do sample=`echo $i|cut -f3 -d/|cut -f1 -d.`;echo -e "$sample\t$i" >> sample.list;done
qsub create_genomics_db.sh	# --genomicsdb-workspace-path为程序自动创建，请勿手动创建
cd ..

# step3.combine_PoN
mkdir -p step3.combine_PoN && cd step3.combine_PoN
qsub combine_PoN.sh


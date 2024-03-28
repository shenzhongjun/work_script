###需要修改,从/mnt/share02/lixx/05.database/05.GTex/1.correlation/01.data/tissue.list中选择
control=Stomach

###无需修改
outdir=`pwd`
tumor=`cat sample_list.txt|cut -f1|tail -1`
fq1=$outdir/clean_data/$tumor/$tumor_1.clean.fastq.gz
fq2=$outdir/clean_data/$tumor/$tumor_2.clean.fastq.gz

mkdir -p $outdir/STAR_EBseq/01.STAR
cd $outdir/STAR_EBseq/01.STAR

/mnt/share01/tools/analysis_module/rna_analysis/STAR/STAR-2.5.2b/bin/Linux_x86_64/STAR \
	--runThreadN 16 \
	--genomeDir /mnt/share01/tools/analysis_module/rna_analysis/star_fusion/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/ref_genome.fa.star.idx \
	--readFilesCommand zcat \
	--readFilesIn $fq1 \
	$fq2 \
	--outFileNamePrefix  $tumor \
	--outSAMtype BAM SortedByCoordinate \
	--limitBAMsortRAM 70000000000 \
	--outBAMsortingThreadN 16 \
	--quantMode GeneCounts

mkdir -p $outdir/STAR_EBseq/02.DEG
cd $outdir/STAR_EBseq/02.DEG

star_result=$outdir/STAR_EBseq/01.STAR/*ReadsPerGene.out.tab

#python /mnt/share02/lixx/05.database/05.GTex/3.other/STAR_EBseq/1.script/star_gtex_merge_EBseq.py \
python /mnt/share02/lixx/05.database/05.GTex/3.test/STAR_DEseq2/1.script/star_gtex_merge_EBseq_DEseq2.py \
	--star_ReadsPerGene $star_result \
	--gtex_path /mnt/share02/lixx/05.database/05.GTex/2.GTex_database \
	--tissue_type $control \
	--outdir $outdir/STAR_EBseq/02.DEG \
	--prefix $tumor \
	--method EBseq

/mnt/share01/tools/miniconda/bin/python /mnt/share01/tools/analysis_module/reports/mutation_report/script/zhixuan/rna_expression_process.py \
	--expression_file $outdir/STAR_EBseq/02.DEG/R${tumor}_vs_Control.final.txt \
	--outdir $outdir/STAR_EBseq/02.DEG \
	--sample_id R$tumor \
	--project_id NBBR


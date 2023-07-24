sample=RAFD456
dir=/mnt/share04/clinical_project/projects/blood_tumor/DDN22001210_PN2201180005_auto/RNA_NBBR/
project=DDN22001210

/mnt/share01/tools/miniconda/bin/python \
	/mnt/share01/tools/analysis_module/reports/mutation_report/script/zhixuan/rna_expression_process.py \
	--expression_file $dir/STAR_EBseq/02.DEG/${sample}_vs_Control.final.txt \
	--outdir $dir/STAR_EBseq/02.DEG \
	--project_id NBBR \
	--sample_id $sample
	
perl /mnt/share02/xuefei/clinical_RNA/script/DGEplot_heatmap_volcano_point.pl \
	-p $dir \
	-c Sample:Control \

ln -sf $dir/STAR_EBseq/02.DEG/$sample.known_diff_expression_gene_filter.tsv $dir/STAR_EBseq/02.DEG/$project.known_diff_expression_gene_filter.tsv 


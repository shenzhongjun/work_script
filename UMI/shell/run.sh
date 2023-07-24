echo ==== start at `date "+%F  %H:%M:%S"` ====
#分子标签质控
#python UMI_qc.py --input /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step3.SamToFastq_bwa_MergeBamAlignment/ng10.umi.merged.uBAM --output UMI_count.txt

#双链饱和度质控
#1. 以原始reads的0.1，0.2...0.9倍共9个梯度挑选reads，生成多个bam
#python pick_bam.py --input /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step4.GroupReadsByUmi/ng10.umi.group.sorted.bam --outDir /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step4.GroupReadsByUmi 
#2. 生成的bam重新分析第4步
#for i in $(seq 0 8);do python /mnt/share02/wangwp/scripts/project/UMI/run_call_duplex.py --input /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step4.GroupReadsByUmi/pick_${i}.bam --sample pick_${i} --outdir /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo;done
#for i in $(seq 0 8);do qsub -cwd -q development -l p=10 /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step5.CallDuplexConsensusReads/run_pick_${i}.sh;done
#3. 统计梯度提取的原始reads数，生成的duplex数，并使用lowess计算 duplex0.5，根据设置的cutoff进行质控状态的判定
#python PCTcount.py --cbam /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step5.CallDuplexConsensusReads/ng10.consensus.uBAM --rbam /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step4.GroupReadsByUmi/ng10.umi.group.sorted.bam --rbamDir /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step4.GroupReadsByUmi --cbamDir /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step5.CallDuplexConsensusReads --outfile PCT_status.txt

#插入片段质控
python /mnt/share02/wangwp/scripts/project/UMI/instert_size_qc.py --bam /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step3.SamToFastq_bwa_MergeBamAlignment/ng10.umi.merged.uBAM --bed  /mnt/share02/wangwp/scripts/project/UMI/bed/merge.bed --InsertSize ng10.insertSize.txt --outfile insertSize_QC.txt --png insertSize.png
#使用脚本测试结果
#python quality_control_of_panel.py --sample ng10 --bam /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step3.SamToFastq_bwa_MergeBamAlignment/ng10.umi.markdup.bam --bed /mnt/share02/wangwp/scripts/project/UMI/bed/merge.bed --json /mnt/share05/clinical_project/projects/blood_tumor/xuefei_test/UMI_naangda/demo/UMI_NanoPrep_10ng_illumina_Demo/step0.fastp/test.json --qcFile qc.xls 
echo ==== end__ at `date "+%F  %H:%M:%S"` ====

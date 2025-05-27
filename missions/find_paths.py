"""科研任务，林俊忠肠癌项目计算CIN ITH评分，查找需要的文件"""
import os
import csv

# 输入文件路径
weslist_file = "weslist"
base_dir = "/mnt/share05/clinical_project/projects/blood_tumor/test/zhouyj_test/research/LJZ_CIN_ITH/wesdata"
output_wes = "wes.list.txt"
output_panel = "panel.list.txt"
output_wes_cnv2 = 'wes.cnv2.list.txt'
output_wes_merge = 'wes.merge.list.txt'
output_wes_huizong = 'wes.huiizong.list.txt'
missing_file = "missing_files.txt"
wes_cnv_data = "wes.cnv.data.txt"

# 初始化输出文件
with open(output_wes, "w") as wes_out, open(output_panel, "w") as panel_out, open(missing_file, "w") as miss, open(wes_cnv_data, "w", newline="") as wes_cnv_out, open(output_wes_merge, 'w') as wes_merge, open(output_wes_huizong, 'w') as wes_huizong, open(output_wes_cnv2, 'w') as wes_cnv2:
    wes_cnv_writer = csv.writer(wes_cnv_out, delimiter="\t")
    wes_cnv_writer.writerow(["chromosome", "start", "end", "cn", "samplename"])  # 写入表头

    # 遍历 weslist 文件
    with open(weslist_file, "r") as f:
        for line in f:
            # 读取每一行的 ID 和姓名
            parts = line.strip().split("\t")
            if len(parts) != 2:
                miss.write(f"格式错误: {line.strip()}\n")
                continue

            order, name = parts
            folder = os.path.join(base_dir, f"{order}_{name}")

            # 检查文件夹是否存在
            if not os.path.isdir(folder):
                miss.write(f"文件夹不存在: {folder}\n")
                continue

            # 查找文件路径
            cnvkit_path = os.path.join(folder, "DNA_Ncet_NCP980", "cnvkit")
            found_files = []
            huizong_path = os.path.join(folder, 'DNA_Ncet_NCP980', 'single_wes', 'new', '检测结果汇总.txt')

            if os.path.isfile(huizong_path):
                wes_huizong.write(f'{huizong_path}\n')

            if os.path.isdir(cnvkit_path):
                for root, _, files in os.walk(cnvkit_path):
                    for file in files:
                        if file.endswith(".genemetrics.addCN.xls"):
                            found_files.append(os.path.join(root, file))

            if found_files:
                for file_path in found_files:
                    try:
                        # 检查文件行数
                        with open(file_path, "r") as file:
                            line_count = sum(1 for _ in file)

                        if line_count > 1000:
                            wes_out.write(f"{file_path}\n")  # 输出到 WES 文件

                            # 获取样本名，是路径中 cnvkit 文件夹后的一级目录
                            sample_name = os.path.basename(os.path.dirname(file_path))
                            print(sample_name, file_path, os.path.dirname(file_path))

                            # 读取 WES 文件内容并提取指定列
                            with open(file_path, "r") as cnv_file:
                                cnv_reader = csv.reader(cnv_file, delimiter="\t")
                                next(cnv_reader)  # 跳过表头
                                for row in cnv_reader:
                                    # 提取第 2, 3, 4, 9 列，加上文件名
                                    try:
                                        wes_cnv_writer.writerow([row[1], row[2], row[3], row[8], sample_name])
                                    except IndexError:
                                        continue  # 跳过不完整行
                            # cnv结果和merge表路径写入文件
                            cnv2_path = os.path.join(folder, 'DNA_Ncet_NCP980', 'cnvkit', sample_name, f'{sample_name}.final.call.cns')
                            wes_cnv2.write(f'{cnv2_path}\n')
                            merge_path = os.path.join(folder, 'DNA_Ncet_NCP980', 'single_wes', 'new', f'{order}_merge_new.txt')
                            wes_merge.write(f"{merge_path}\n")

                            # 创建文件夹，写入命令
                            if not os.path.exists(f'ITH/{order}_{name}'):
                                os.makedirs(f'ITH/{order}_{name}')
                            cmd = f"""
python /mnt/share02/zhouyj/script/mrd/prepare.PyClone-VI.py --merge {merge_path} --cnvkit {cnv2_path} --tissue_id {sample_name} --purity 0.8
/mnt/share02/lixx/02.software/envs_python3.7/PyClone-VI/bin/pyclone-vi fit -i {sample_name}.pyclone-vi.txt -o {sample_name}.pyclone-vi.h5  -c 40 -d beta-binomial -r 10
/mnt/share02/lixx/02.software/envs_python3.7/PyClone-VI/bin/pyclone-vi write-results-file -i {sample_name}.pyclone-vi.h5 -o {sample_name}.pyclone-vi.output.txt
python /mnt/share02/zhouyj/script/missions/cal_ith.py {sample_name}.pyclone-vi.output.txt > {sample_name}.ith.txt
"""
                            with open(f'ITH/{order}_{name}/{sample_name}.cmd.sh', 'w') as wcmd:
                                wcmd.write(cmd)
                        else:
                            panel_out.write(f"{file_path}\n")  # 输出到 Panel 文件
                    except Exception as e:
                        miss.write(f"无法读取文件: {file_path} 错误: {str(e)}\n")
            else:
                miss.write(f"目标文件缺失: {folder}\n")


print(f"文件查找完成。\nWES 样本保存在 {output_wes}。\nWES cnv cns保存在{output_wes_cnv2}。\nPanel 样本保存在 {output_panel}。\n缺失信息保存在 {missing_file}。\nMerge表保存在{output_wes_merge}。\n检测结果汇总保存在{output_wes_huizong}。")



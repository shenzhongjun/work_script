"""科研任务，林俊忠肠癌项目计算CIN ITH评分，查找需要的文件"""
import sys
import pandas as pd

# 1. 读取 PyClone-VI 输出文件
file_path = sys.argv[1]   # "ALI946.pyclone-vi.output.txt"
data = pd.read_csv(file_path, sep="\t")

# 2. 按 cluster_id 分组，计算每个克隆的平均 cellular_prevalence 和突变数量
cluster_stats = data.groupby("cluster_id").agg(
    cellular_prevalence_mean=("cellular_prevalence", "mean"),
    mutation_count=("mutation_id", "count")
).reset_index()

# 3. 查找最大和次大 CCF 的克隆
sorted_clusters = cluster_stats.sort_values(by="cellular_prevalence_mean", ascending=False)
largest_cluster = sorted_clusters.iloc[0]  # 最大 CCF 克隆
second_largest_cluster = sorted_clusters.iloc[1]  # 次大 CCF 克隆

# 4. 判断是否需要合并克隆
if largest_cluster["mutation_count"] == 1:
    # 如果最大 CCF 的克隆仅有一个突变，合并最大和次大克隆
    main_clone_ids = [largest_cluster["cluster_id"], second_largest_cluster["cluster_id"]]
else:
    # 否则，主克隆为最大 CCF 的克隆
    main_clone_ids = [largest_cluster["cluster_id"]]

# 5. 统计主克隆和子克隆的突变数量
n_main = data[data["cluster_id"].isin(main_clone_ids)]["mutation_id"].count()
n_sub = data[~data["cluster_id"].isin(main_clone_ids)]["mutation_id"].count()

# 6. 计算 ITH 分数
ith_score = n_sub / (n_main + n_sub)

# 7. 输出结果
print(f"主克隆 ID: {main_clone_ids}")
print(f"主克隆突变数量 (n_main): {n_main}")
print(f"子克隆突变数量 (n_sub): {n_sub}")
print(f"ITH 分数: {ith_score:.4f}")

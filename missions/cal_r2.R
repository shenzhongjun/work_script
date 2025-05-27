# 科研任务，林俊忠肠癌项目计算R2
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggExtra)

prefix = 'Ep'
seurat_obj = readRDS('Ep.rds')
percent = read.table('样本_编号.对应关系2')
percent$Percent <- as.numeric(sub("%", "", percent$V3)) / 100
percent$Percent[percent$Percent == 0] <- NA
state1 <- data.frame(barcode = rownames(seurat_obj@meta.data), Sample = seurat_obj$orig.ident, Cluster = seurat_obj$seurat_clusters)
sum <- tapply(state1$Sample, state1$Cluster, table)
matrix1 <- data.frame()
for (i in names(sum)){
  tmp <- data.frame(sum[[i]])
  tmp$cluster <- i
  matrix1 <- rbind(matrix1, tmp)
}
colnames(matrix1) <- c("Sample", "nCells", "Cluster")
matrix1$nCells <- as.numeric(as.character(matrix1$nCells))
matrix1 <- ddply(matrix1, "Sample", transform, Percent = round(nCells/sum(nCells) * 100, 2))
write.table(matrix1, file =paste(prefix,"nCells-percent_samplestats.txt",sep="."), row.names = F, col.names = T, quote = F, sep = "\t")

result <- merge(matrix1, percent[, c("V2", "Percent")], by.x = "Sample", by.y = "V2", all.x = TRUE)
result$Percent.x <- result$Percent.x / 100


# ======循环=====
# 获取所有cluster的编号
all_clusters <- unique(result$Cluster)

# 定义一个函数来进行相关性分析和绘图，并返回相关性结果
plot_correlation <- function(result_subset, cluster_names, filename_prefix) {
  # 计算相关性和 p 值

  cor_test <- cor.test(result_subset$Percent.x, result_subset$Percent.y, use = "complete.obs")
  r_value <- round(cor_test$estimate, 2)
  p_value <- round(cor_test$p.value, 3)
  # # 创建散点图，添加回归线和置信区间
  # p <- ggplot(result_subset, aes(x = Percent.x, y = Percent.y)) +
    # geom_point(color = "brown", size = 2) +
    # geom_smooth(method = "lm", color = "black", fill = "gray", level = 0.95) +
    # labs(x = paste0("Percent.x (", paste(cluster_names, collapse = ", "), ")"),
         # y = "Percent.y (Tumor size reduction)") +
    # theme_minimal(base_size = 15) +
    # annotate("text",
             # x = max(result_subset$Percent.x, na.rm = TRUE) * 0.8,
             # y = max(result_subset$Percent.y, na.rm = TRUE) * 0.95,
             # label = paste0("p = ", p_value, "\n", "r = ", r_value),
             # size = 5, hjust = 0, color = "blue")

  # # 添加边缘直方图
  # p_with_marginal <- ggMarginal(
    # p,
    # type = "histogram",
    # fill = "gold",
    # bins = 30,
    # margins = "both",
    # size = 7
  # )

  # # 保存图像
  # ggsave(paste0(filename_prefix, ".png"), p_with_marginal)

  # 返回相关性分析结果
  return(c(p_value, r_value))
}

# 创建一个空的结果数据框来存储所有组合的p_value和r_value
correlation_results <- data.frame(Cluster_Combination = character(0), p_value = numeric(0), r_value = numeric(0))

# 循环处理1到19个cluster的组合
for (n in 1:19) {
  # 生成所有n个Cluster的组合
  cluster_combinations <- combn(all_clusters, n, simplify = FALSE)

  # 遍历每个组合进行分析
  for (cluster_comb in cluster_combinations) {
    result_subset <- subset(result, Cluster %in% cluster_comb)

    # 定义文件名前缀
    filename_prefix <- paste0("Cluster_", paste(cluster_comb, collapse = "_"))

    # 调用绘图函数并获取相关性结果
    cor_values <- plot_correlation(result_subset, cluster_comb, filename_prefix)

    # 将结果添加到数据框
    correlation_results <- rbind(correlation_results,
                                 data.frame(Cluster_Combination = paste(cluster_comb, collapse = "_"),
                                            p_value = cor_values[1], r_value = cor_values[2]))
  }
}

# 将所有的p_value和r_value输出到文件
write.table(correlation_results[, c("Cluster_Combination", "p_value", "r_value")],
            file = "correlation_results.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)



# --------无作图，实时输出------------
# 获取所有cluster的编号
all_clusters <- unique(result$Cluster)

# 创建一个空的数据框来存储所有组合的p_value和r_value
correlation_results <- data.frame(p_value = numeric(0), r_value = numeric(0), Cluster_Combination = character(0))

# 打开文件连接，以便实时写入
output_file <- file("correlation_results_new.txt", open = "w")
# 写入文件的表头
writeLines("p_value\tr_value\tCluster_Combination", output_file)

# 循环处理1到19个cluster的组合
for (n in 1:19) {
  # 生成所有n个Cluster的组合
  cluster_combinations <- combn(all_clusters, n, simplify = FALSE)

  # 遍历每个组合进行分析
  for (cluster_comb in cluster_combinations) {
    result_subset <- subset(result, Cluster %in% cluster_comb)
    result_summary <- result %>%
    group_by(Sample) %>%
    summarise(
      total_nCells = sum(nCells),                                     # 计算每个样本的总细胞数
      Cluster_nCells = sum(nCells[Cluster %in% cluster_comb]),           # 计算 Cluster 1+2 的细胞数
      Cluster_ratio = Cluster_nCells / total_nCells,              # 计算比例
      Percent = unique(Percent.y)
    )


    # 计算相关性和 p 值
    cor_test <- cor.test(result_summary$Cluster_ratio, result_summary$Percent, method="spearman")
    # cor_test <- cor.test(result_subset$Percent.x, result_subset$Percent.y, use = "complete.obs")
    r_value <- round(cor_test$estimate, 2)
    p_value <- round(cor_test$p.value, 3)

    # 输出到屏幕
    cat("Cluster Combination:", paste(cluster_comb, collapse = "_"),
        "r_value:", r_value,
        "p_value:", p_value, "\n")

    # 将结果添加到数据框
    correlation_results <- rbind(correlation_results,
                                 data.frame(p_value = p_value, r_value = r_value,
                                            Cluster_Combination = paste(cluster_comb, collapse = "_")))

    # 实时输出到文件
    writeLines(paste(p_value, r_value, paste(cluster_comb, collapse = "_"), sep = "\t"), output_file)
  }
}

# 关闭文件连接
close(output_file)



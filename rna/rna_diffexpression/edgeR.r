#!/home/zhouyj/anaconda3/envs/R4/bin/Rscript

library(getopt)

spec <- matrix(c(
    'help',      'h',  0,  'logical',    'Show this help.',
    'sampledir', 's',  1,  'character',  'The sample dir.',
    'controldir',   'c',  1,  'character',  'The control dir.',
    'sample',    'a',  1,  'character',  'The sample name.',
    'control',   'b',  1,  'character',  'The control name.',
    'outdir',    'o',  1,  'character',  'The out dir.'
), byrow=T, ncol=5)
opt = getopt(spec=spec)

# 输入参数检查
if( !is.null(opt$help) || is.null(opt$sampledir) || is.null(opt$controldir) || is.null(opt$sample) || is.null(opt$control) || is.null(opt$out)) {
    cat(getopt(spec=spec, usage = T), '\n')
    quit()
}

# 参数设置
opt$outdir = sub("/$", "", opt$outdir)    # 去掉末尾的'/'防止报错
opt$sampledir = sub("/$", "", opt$sampledir)
opt$controldir = sub("/$", "", opt$controldir)

# START
library(glue)
library(edgeR)
library(tidyverse)

outdir = normalizePath(opt$outdir)
sampledir = normalizePath(opt$sampledir)
controldir = normalizePath(opt$controldir)

sample_name = opt$sample
control_name = opt$control
dir.create(glue('{outdir}/{sample_name}_vs_{control_name}'), showWarnings=F)
setwd(glue('{outdir}/{sample_name}_vs_{control_name}'))

print(glue('开始读取表达矩阵等文件并进行差异基因分析...{Sys.time()}'))
sample_readcount = read.table(glue('{sampledir}/merged_readcount'), header = T, check.names = F, col.names = c('基因名', sample_name))
control_readcount = read.table(glue('{controldir}/merged_readcount'), header = T, check.names = F, col.names = c('基因名', control_name))
sample_fpkm = read.table(glue('{sampledir}/merged_fpkm'), header = T, check.names = F, col.names = c('基因名', sample_name))
control_fpkm = read.table(glue('{controldir}/merged_fpkm'), header = T, check.names = F, col.names = c('基因名', control_name))
id2name_table = read.table(glue('{controldir}/geneid2genename.txt'), header = F, col.names=c('id', '基因名'))
id2name_table$基因id = do.call(rbind, strsplit(id2name_table$id, ':'))[,2]    # 拆分字符串获得基因id
id2name_table$id = NULL
append_info_table = reduce(list(sample_fpkm, control_fpkm, id2name_table), full_join, by='基因名')

raw_count_table = full_join(sample_readcount, control_readcount, by='基因名')
rownames(raw_count_table) = raw_count_table[, 1]
raw_count_table[,1] = NULL

dge = DGEList(counts=data.matrix(raw_count_table), group = 1:2)
keep = rowSums(cpm(dge)>1) >= 1  # 保留在至少在一个样本里有表达的基因(CPM > 1)
dge = dge[keep, , keep.lib.sizes=FALSE]
dge = calcNormFactors(dge)

normed_count_table = data.frame(matrix(ncol = 3, nrow = nrow(dge$counts)))
colnames(normed_count_table) = c('基因名', 'sampleNCounts', 'controlNCounts')
normed_count_table$基因名 = rownames(dge$counts)
normed_count_table$sampleNCounts = dge$samples$norm.factors[1]*dge$counts[, 1]
normed_count_table$controlNCounts = dge$samples$norm.factors[2]*dge$counts[, 2]
append_info_table = inner_join(append_info_table, normed_count_table, by='基因名')

bcv = 0.2
et = exactTest(dge, dispersion = bcv ^ 2)
sign_genes = decideTestsDGE(et, p.value = 0.05, lfc = 2)
result = topTags(et, n = nrow(et$table))$table
result$基因名 = rownames(result)
result$基因上下调 = ifelse(result$logFC > 0, "上调", "下调")
result$是否差异显著基因 = ifelse(result$FDR < 0.05, "是", "否")
result = inner_join(result, append_info_table, by='基因名')
result = subset(result, select=c(5, 10, 11, 12, 8, 9, 1, 3, 4, 6, 7))
colnames(result) = c('基因名', '基因id',
                     glue('标准化后的序列支持数:{sample_name}'), glue('标准化后的序列支持数:{control_name}'),
                     glue('FPKM值:{sample_name}'), glue('FPKM值:{control_name}'),
                     'logFC', 'PValue', 'FDR', '基因上下调', '是否差异显著基因')

deg = result[which(result$是否差异显著基因=='是'), -11]
deg_up = deg[which(deg$基因上下调=='上调'), -10]
deg_down = deg[which(deg$基因上下调=='下调'), -10]

write.table(result, file=glue('{sample_name}_vs_{control_name}.final.txt'), sep="\t", quote=F, row.names=F)
write.table(deg, file=glue('{sample_name}_vs_{control_name}.final.DEG.txt'), sep="\t", quote=F, row.names=F)
write.table(deg_up, file=glue('{sample_name}_vs_{control_name}.final.DEG.Up.txt'), sep="\t", quote=F, row.names=F)
write.table(deg_down, file=glue('{sample_name}_vs_{control_name}.final.DEG.Down.txt'), sep="\t", quote=F, row.names=F)

# 富集分析
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

dir.create("enrichment", showWarnings=F)
setwd("enrichment")
print(glue('开始富集分析...{Sys.time()}'))

enrichment = function(enrichData, type) {
    go = enrichGO(gene=enrichData$基因id, OrgDb="org.Hs.eg.db", ont="ALL", readable=T)
    write.table(go@result,
                file=glue('{sample_name}_vs_{control_name}.enrichGO.{type}.txt'),
                sep="\t", quote=F, row.names=F)
    dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    ggsave(glue('{sample_name}_vs_{control_name}.enrichGO.{type}.pdf'), width=12, height=16)
    ggsave(glue('{sample_name}_vs_{control_name}.enrichGO.{type}.png'), width=12, height=16)

    options(clusterProfiler.download.method = "wget")
    kegg = enrichKEGG(gene=enrichData$基因id, keyType = 'ncbi-geneid')    # 如果运行失败将阈值放宽到1:, pvalueCutoff = 1, qvalueCutoff = 1
    kegg = setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
    write.table(kegg@result,
                file=glue('{sample_name}_vs_{control_name}.enrichKEGG.{type}.txt'),
                sep="\t", quote=F, row.names=F)
    dotplot(kegg, showCategory=20)
    ggsave(glue('{sample_name}_vs_{control_name}.enrichKEGG.{type}.pdf'), width=12, height=12)
    ggsave(glue('{sample_name}_vs_{control_name}.enrichKEGG.{type}.png'), width=12, height=12)
}

enrichment(deg_up, 'Up')
enrichment(deg_down, 'Down')

# kegg GSEA分析（Gene Set Enrichment Analysis），用上调+下调全部基因进行分析
log2fc = deg$logFC  # 提取log2fc值为向量
names(log2fc) = as.character(deg$基因id)    # add vector name
sortedFc = sort(log2fc, decreasing=T)     # decreasing order，按值大小降序排列
gseakegg = gseKEGG(geneList=sortedFc, keyType="ncbi-geneid")
gseakegg = setReadable(gseakegg, 'org.Hs.eg.db', 'ENTREZID')
write.table(gseakegg@result, file=glue("{sample_name}_vs_{control_name}.gseaKEGG.txt"), sep="\t", quote=F, row.names=F)

if ((nrow(gseakegg) >= 2) && (nrow(gseakegg) <= 10)) {
    gseaplot2(gseakegg, seq_len(nrow(gseakegg)), color="red", base_size=15, rel_heights = c(1.5, 0.2, 0.2))
    ggsave(glue("{sample_name}_vs_{control_name}.GSEA_KEGG.all.pdf"), width=12, height=12)
    ggsave(glue("{sample_name}_vs_{control_name}.GSEA_KEGG.all.png"), width=12, height=12)
}
for (path in gseakegg$ID) {
    title = gseakegg[gseakegg$ID == path, 2]
    gseaplot2(gseakegg, title = title, path, color="red", base_size=15, pvalue_table=T)
    ggsave(glue("{sample_name}_vs_{control_name}.GSEA_KEGG.{path}.pdf"), width=12, height=10)
    ggsave(glue("{sample_name}_vs_{control_name}.GSEA_KEGG.{path}.png"), width=12, height=10)
}

# 画热图
library(pheatmap)

print(glue('开始画热图...{Sys.time()}'))
plot_heatmap = function(heatData, type) {
    heatData = heatData[order(heatData$FDR), ]
    rownames(heatData) = heatData[ ,1]
    top20 = heatData[1:20, c(5, 6)]
    colnames(top20) = c(sample_name, control_name)
    p = pheatmap(top20, scale="row",
                 main=glue('Top20 {type} Regulated Genes'),
                 border="white", # 设置边框为白色
                 cluster_cols=F,  # 去掉列聚类
                 cluster_rows=F,  # 去掉行聚类
                 angle_col=0, # 设置显示角度
                 fontsize=15
                 )
    ggsave(glue('{sample_name}_vs_{control_name}.heatmap.{type}.pdf'), p, width=12, height=12)
    ggsave(glue('{sample_name}_vs_{control_name}.heatmap.{type}.png'), p, width=12, height=12)
}
setwd(glue('{outdir}/{sample_name}_vs_{control_name}'))
plot_heatmap(deg_up, 'Up')
plot_heatmap(deg_down, 'Down')

print(glue('任务结束...{Sys.time()}'))
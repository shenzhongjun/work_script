#!/home/zhouyj/anaconda3/envs/R4/bin/Rscript

library(getopt)

spec <- matrix(c(
    'help',     'h',  0,  'logical',    'Show this help.',
    'degfile',  'i',  1,  'character',  'The diff expression file like "R1REPZT11_vs_R1REPZTK-.final.txt" file path.',
    'id2name',  't',  1,  'character',  'The geneid2genename file path.',
    'sample',   'a',  1,  'character',  'The sample name.',
    'control',  'b',  1,  'character',  'The control name.',
    'outdir',   'o',  1,  'character',  'The out dir.'
), byrow=T, ncol=5)
opt = getopt(spec=spec)

# 输入参数检查
if( !is.null(opt$help) || is.null(opt$degfile) || is.null(opt$sample) || is.null(opt$control) || is.null(opt$out)) {
    cat(getopt(spec=spec, usage = T), '\n')
    quit()
}

# 参数设置
opt$outdir = sub("/$", "", opt$outdir)    # 去掉末尾的'/'防止报错
degfile = normalizePath(opt$degfile)    # 获取绝对路径
id2name = normalizePath(opt$id2name)
outdir = normalizePath(opt$outdir)
sample_name = opt$sample
control_name = opt$control

# START
library(glue)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(pheatmap)

setwd(outdir)

deg_raw_table = read.table(degfile, header=T, check.names=F)
id2name_table = read.table(id2name, col.names=c('id', '基因名'))
deg_raw_table = merge(deg_raw_table, id2name_table, by="基因名", all.x = TRUE)   # 合并获得基因id字符串
deg_raw_table$基因id = do.call(rbind, strsplit(deg_raw_table$id, ':'))[,2]    # 拆分字符串获得基因id
deg_raw_table = subset(deg_raw_table, select=c(1, 11, 2, 3, 4, 5, 6, 7, 8, 9))     # 列重新排列
deg_table = deg_raw_table[which(deg_raw_table$是否差异显著基因 == '是'), ]
deg_up = deg_table[deg_table['基因上下调']=='上调', -c(9, 10)]
deg_down = deg_table[deg_table['基因上下调']=='下调', -c(9, 10)]

#
# write.table(deg_up,
#             file=glue("{sample_name}_vs_{control_name}.DEG.final.Up.txt"),
#             sep="\t", quote=F, row.names=F)
# write.table(deg_down,
#             file=glue("{sample_name}_vs_{control_name}.DEG.final.Down.txt"),
#             sep="\t", quote=F, row.names=F)
#
# # 富集分析
# dir.create("enrichment", showWarnings=F)
# setwd("enrichment")
#
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
#
# # kegg GSEA分析（Gene Set Enrichment Analysis），用上调+下调全部基因进行分析
log2fc = deg_raw_table$"log2(FoldChange)值"  # 提取log2fc值为向量
names(log2fc) = as.character(deg_raw_table$基因id)    # add vector name
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
plot_heatmap = function(heatData, type) {
    heatData = heatData[order(heatData$FDR值), ]
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

plot_heatmap(deg_up, 'Up')
plot_heatmap(deg_down, 'Down')

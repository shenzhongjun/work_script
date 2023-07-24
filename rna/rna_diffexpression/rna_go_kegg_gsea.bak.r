#!/home/zhouyj/anaconda3/envs/R4/bin/Rscript

# R路径：/home/zhouyj/anaconda3/envs/R4/bin/R

library(glue)
library(getopt)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(pheatmap)

spec <- matrix(c(
    'help',     'h',  0,  'logical',    'Show this help.',
    'degfile',  'i',  1,  'character',  'The DEG file path.',
    'id2name',  't',  1,  'character',  'The geneid2genename file path.',
    'sample',   'a',  1,  'character',  'The sample name.',
    'control',  'b',  1,  'character',  'The control name.',
    'outdir',   'o',  1,  'character',  'The out dir.'
), byrow=T, ncol=5)
opt = getopt(spec=spec)

# 将usage参数打开，getpot()就会返回一个特定写法的帮助文件
if( !is.null(opt$help) || is.null(opt$degfile) || is.null(opt$sample) || is.null(opt$control) || is.null(opt$out)) {
    cat(getopt(spec=spec, usage = T), '\n')
    quit()
}

# 参数路径去掉末尾的'/'防止报错
opt$degfile = sub("/$","",opt$degfile)
opt$outdir = sub("/$","",opt$outdir)

degfile = normalizePath(opt$degfile)
id2name = normalizePath(opt$id2name)
outdir = normalizePath(opt$outdir)
setwd(outdir)

sample = opt$sample
control = opt$control

deg_table = read.table(degfile, header=T, check.names=F)
id2name_table = read.table(id2name, col.names=c('id', '基因名'))

deg_table = merge(deg_table, id2name_table, by="基因名", all.x = TRUE)   # 合并获得基因id字符串
deg_table$基因id = do.call(rbind, strsplit(deg_table$id, ':'))[,2]    # 拆分字符串获得基因id
deg_table = subset(deg_table, select=c(1, 10, 2, 3, 4, 5, 6, 7, 8))     # 列重新排列
deg_up = deg_table[deg_table['基因上下调']=='上调', -9]
deg_down = deg_table[deg_table['基因上下调']=='下调', -9]

write.table(deg_up, file=glue("{sample}_vs_{control}.DEG.final.up.txt"), sep="\t", quote=F, row.names=F)
write.table(deg_down, file=glue("{sample}_vs_{control}.DEG.final.down.txt"), sep="\t", quote=F, row.names=F)

dir.create("GOandPathway", showWarnings=F)
setwd("GOandPathway")

# go = enrichGO(myDMG$Entre_ID, OrgDb="org.Hs.eg.db", keyType="ENTREZID",
              # ont="ALL",pvalueCutoff=0.05, qvalueCutoff=0.2, readable=T)
# kegg = enrichKEGG(organism='hsa', gene=myDMG$Entre_ID, pvalueCutoff=0.05)

# 函数默认参数如下：
# enrichGO(gene, OrgDb, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe,
#   qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE)
# enrichKEGG(gene, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe, minGSSize = 10,
#   maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
# gseKEGG(geneList, organism = "hsa", keyType = "kegg", exponent = 1, minGSSize = 10, maxGSSize = 500, eps = 1e-10, pvalueCutoff = 0.05,
#   pAdjustMethod = "BH", verbose = TRUE, use_internal_data = FALSE, seed = FALSE, by = "fgsea", ...)
go = enrichGO(gene=deg_up$基因id, OrgDb="org.Hs.eg.db", ont="ALL", readable=T)
options(clusterProfiler.download.method = "wget")
kegg = enrichKEGG(gene=deg_up$基因id, keyType = 'ncbi-geneid')
kegg = setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')

# kegg GSEA分析（Gene Set Enrichment Analysis），用上调+下调差异显著基因进行
log2fc = deg_table$"log2(FoldChange)值"  # 提取log2fc值为向量
names(log2fc) = as.character(deg_table$基因id)    # add vector name
sortedFc = sort(log2fc, decreasing=T)     # decreasing order，按值大小降序排列
gseakegg = gseKEGG(geneList=sortedFc, keyType="ncbi-geneid")
gseakegg = setReadable(gseakegg, 'org.Hs.eg.db', 'ENTREZID')

write.table(go@result, file="enrichGO.txt", sep="\t", quote=F, row.names=F)
write.table(kegg@result, file="enrichKEGG.txt", sep="\t", quote=F, row.names=F)
write.table(gseakegg@result, file="gseaKEGG.txt", sep="\t", quote=F, row.names=F)

# barplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")  # 暂无需条形图
# ggsave('go_bar_plot.pdf', width=12, height=16)
# cnetplot(go, showCategory=5)  # 暂无需Gene-Concept Network图
# cnetplot(go, showCategory=10, circular=TRUE,colorEdge=TRUE) # 圆形布局，给线条上色
# ggsave("go_cnet_plot.pdf",  width=16, height=16)
dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
ggsave("go_dot_plot.pdf",  width=12, height=16)
# barplot(kegg, showCategory=20, drop=T)
# ggsave("kegg_bar_plot.pdf", width=12, height=12)
dotplot(kegg, showCategory=20)
ggsave("kegg_dot_plot.pdf", width=12, height=12)

# dotplot(gseakegg)     # 暂无需气泡图
# ggsave("gsea_kegg_dot_plot.pdf", width=12, height=12)

if (nrow(gseakegg) <= 10) {
    gseaplot2(gseakegg, 1:nrow(gseakegg), color="red", base_size=15, rel_heights = c(1.5, 0.2, 0.2))
    ggsave(glue("gsea_kegg_all_plot.pdf"), width=12, height=12)
    ggsave(glue("gsea_kegg_all_plot.png"), width=12, height=12)
}
for (path in gseakegg$ID) {
    title = gseakegg[gseakegg$ID == path, 2]
    gseaplot2(gseakegg, title = title, path, color="red", base_size=15, pvalue_table=T)
    ggsave(glue("gsea_kegg_{path}_plot.pdf"), width=12, height=10)
    ggsave(glue("gsea_kegg_{path}_plot.png"), width=12, height=10)      # , scale = 1.2
}
# gseaplot2默认参数
# gseaplot2(x, geneSetID, title = "", color = "green", base_size = 11,
#     rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE,
#     ES_geom = "line")
# gseaplot2参数解析
# gseaplot2(gseakegg,
#           title = "Salivary secretion",  #设置title
#           "hsa04970", #绘制hsa04740通路的结果
#           color="red", #线条颜色
#           base_size = 20, #基础字体的大小
#           subplots = 1:2, #展示上2部分
#           pvalue_table = T) # 显示p值

# 画热图
plot_heatmap = function(heatdata, type) {
    heatdata = heatdata[order(heatdata$FDR值), ]
    rownames(heatdata) = heatdata[ ,1]
    top20 = heatdata[1:20, c(5, 6)]
    colnames(top20) = c(sample, control)
    p = pheatmap(top20,
                 main=glue('Top20 {type} Regulated Genes'),
                 border="white", # 设置边框为白色
                 cluster_cols=F,  # 去掉列聚类
                 cluster_rows=F,  # 去掉行聚类
                 angle_col=0, # 设置显示角度
                 fontsize=15
                 )
    ggsave(glue('{sample}_vs_{control}.heatmap.{type}.pdf'), p, width=12, height=12)
    ggsave(glue('{sample}_vs_{control}.heatmap.{type}.png'), p, width=12, height=12)
}

plot_heatmap(deg_up, 'Up')
plot_heatmap(deg_down, 'Down')

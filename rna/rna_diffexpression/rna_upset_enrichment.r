#!/home/zhouyj/anaconda3/envs/R4/bin/Rscript

library(getopt)

spec <- matrix(c(
  'help', 'h', 0, 'logical', 'Show this help.',
  'vennlist', 'n', 1, 'character', 'The vennlist file path.',
  'outdir', 'o', 1, 'character', 'The out dir.'
), byrow = T, ncol = 5)
opt = getopt(spec = spec)

# 输入参数检查
if (!is.null(opt$help) ||
  is.null(opt$vennlist) ||
  is.null(opt$outdir)) {
  cat(getopt(spec = spec, usage = T), '\n')
  quit()
}

# 参数设置
opt$outdir = sub("/$", "", opt$outdir)    # 去掉末尾的'/'防止报错
vennlist = normalizePath(opt$vennlist)    # 获取绝对路径
wd = normalizePath(opt$outdir)

# START
library(glue)
library(UpSetR)

setwd(wd)

vl = read.table(vennlist, check.names = F, col.names=c('sample', 'control'))
samples = vl$sample
control = unique(vl$control)
stopifnot('对照样本只能有一个！' = length(control) == 1)

get_deg_table = function (x, ud) {
    sample_name = x[1]
    control_name = x[2]
    if (ud == 'up') {
        deg_table = read.table(
          glue('{wd}/{sample_name}_vs_{control_name}/{sample_name}_vs_{control_name}.final.DEG.Up.txt'),
          header = T, check.names = F)
    } else {
        deg_table = read.table(
          glue('{wd}/{sample_name}_vs_{control_name}/{sample_name}_vs_{control_name}.final.DEG.Down.txt'),
          header = T, check.names = F)
    }
    return(deg_table)
}

deg_up_list = apply(vl, 1, get_deg_table, 'up')
deg_down_list = apply(vl, 1, get_deg_table, 'down')
names(deg_up_list) = names(deg_down_list) = samples

deg_up_genenames_list = list()
deg_down_genenames_list = list()
for (sample in samples) {
    deg_up_genenames_list[[sample]] = deg_up_list[[sample]]$基因名
    deg_down_genenames_list[[sample]] = deg_down_list[[sample]]$基因名
}
upplot = upset(fromList(deg_up_genenames_list), order.by = "freq", decreasing = TRUE)
downplot = upset(fromList(deg_down_genenames_list), order.by = "freq", decreasing = TRUE)

dir.create(glue('{wd}/venn_enrichment'), showWarnings = F)
setwd(glue('{wd}/venn_enrichment'))

pdf('UpGenes.UpSet.pdf')
print(upplot)
dev.off()
png('UpGenes.UpSet.png', width = 2400, height = 2400, res=300)
print(upplot)
dev.off()
pdf('DownGenes.UpSet.pdf')
print(downplot)
dev.off()
png('DownGenes.UpSet.png', width = 2400, height = 2400, res=300)
print(downplot)
dev.off()

# union_group = rbind(sub_group1_up[, 1:2], sub_group1_down[, 1:2], sub_group2_up[, 1:2], sub_group2_down[, 1:2], sub_group3_up[, 1:2], sub_group3_down[, 1:2])
# up_intersect = intersect(intersect(sub_group1_up$基因名, sub_group2_up$基因名), sub_group3_up$基因名)
# down_intersect = intersect(intersect(sub_group1_down$基因名, sub_group2_down$基因名), sub_group3_down$基因名)
# sub_group1_up_uniq = setdiff(setdiff(sub_group1_up$基因名, sub_group2_up$基因名), sub_group3_up$基因名)
# sub_group2_up_uniq = setdiff(setdiff(sub_group2_up$基因名, sub_group1_up$基因名), sub_group3_up$基因名)
# sub_group3_up_uniq = setdiff(setdiff(sub_group3_up$基因名, sub_group1_up$基因名), sub_group2_up$基因名)
# sub_group1_down_uniq = setdiff(setdiff(sub_group1_down$基因名, sub_group2_down$基因名), sub_group3_down$基因名)
# sub_group2_down_uniq = setdiff(setdiff(sub_group2_down$基因名, sub_group1_down$基因名), sub_group3_down$基因名)
# sub_group3_down_uniq = setdiff(setdiff(sub_group3_down$基因名, sub_group1_down$基因名), sub_group2_down$基因名)

union_group = rbind(sub_group1_up[, 1:2], sub_group1_down[, 1:2], sub_group2_up[, 1:2], sub_group2_down[, 1:2])
up_intersect = intersect(sub_group1_up$基因名, sub_group2_up$基因名)
down_intersect = intersect(sub_group1_down$基因名, sub_group2_down$基因名)
sub_group1_up_uniq = setdiff(sub_group1_up$基因名, sub_group2_up$基因名)
sub_group1_down_uniq = setdiff(sub_group1_down$基因名, sub_group2_down$基因名)
sub_group2_up_uniq = setdiff(sub_group2_up$基因名, sub_group1_up$基因名)
sub_group2_down_uniq = setdiff(sub_group2_down$基因名, sub_group1_down$基因名)

write.table(data.frame(基因名=up_intersect), file=glue('{group}.UpGenes.Intersect.txt'), quote=F, row.names=F)
write.table(data.frame(基因名=down_intersect), file=glue('{group}.DownGenes.Intersect.txt'), quote=F, row.names=F)
write.table(data.frame(基因名=sub_group1_up_uniq), file=glue('R{sample_name1}_vs_R{control_name1}.UpGenes.Uniq.txt'), quote=F, row.names=F)
write.table(data.frame(基因名=sub_group2_up_uniq), file=glue('R{sample_name2}_vs_R{control_name2}.UpGenes.Uniq.txt'), quote=F, row.names=F)
# write.table(data.frame(基因名=sub_group3_up_uniq), file=glue('R{sample_name3}_vs_R{control_name3}.UpGenes.Uniq.txt'), quote=F, row.names=F)
write.table(data.frame(基因名=sub_group1_down_uniq), file=glue('R{sample_name1}_vs_R{control_name1}.DownGenes.Uniq.txt'), quote=F, row.names=F)
write.table(data.frame(基因名=sub_group2_down_uniq), file=glue('R{sample_name2}_vs_R{control_name2}.DownGenes.Uniq.txt'), quote=F, row.names=F)
# write.table(data.frame(基因名=sub_group3_down_uniq), file=glue('R{sample_name3}_vs_R{control_name3}.DownGenes.Uniq.txt'), quote=F, row.names=F)

# 富集分析
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
enrichment = function(enrichData, group, type) {
    go = enrichGO(gene=enrichData$基因id, OrgDb="org.Hs.eg.db", ont="ALL", readable=T, pvalueCutoff=1, qvalueCutoff=1)
    write.table(go@result,
                file=glue('{group}.enrichGO.{type}.txt'),
                sep="\t", quote=F, row.names=F)
    dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    ggsave(glue('{group}.enrichGO.{type}.pdf'), width=12, height=16)
    ggsave(glue('{group}.enrichGO.{type}.png'), width=12, height=16)

    options(clusterProfiler.download.method = "wget")
    kegg = enrichKEGG(gene=enrichData$基因id, keyType = 'ncbi-geneid', pvalueCutoff=1, qvalueCutoff=1)
    kegg = setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
    write.table(kegg@result,
                file=glue('{group}.enrichKEGG.{type}.txt'),
                sep="\t", quote=F, row.names=F)
    dotplot(kegg, showCategory=20)
    ggsave(glue('{group}.enrichKEGG.{type}.pdf'), width=12, height=12)
    ggsave(glue('{group}.enrichKEGG.{type}.png'), width=12, height=12)
}


enrichment(union_group[which(union_group$基因名 %in% up_intersect), ], group, 'up_intersect')
enrichment(union_group[which(union_group$基因名 %in% down_intersect), ], group, 'down_intersect')
enrichment(union_group[which(union_group$基因名 %in% sub_group1_up_uniq), ], group, glue('R{sample_name1}_vs_R{control_name1}.up_uniq'))
enrichment(union_group[which(union_group$基因名 %in% sub_group1_down_uniq), ], group, glue('R{sample_name1}_vs_R{control_name1}.down_uniq'))
enrichment(union_group[which(union_group$基因名 %in% sub_group2_up_uniq), ], group, glue('R{sample_name2}_vs_R{control_name2}.up_uniq'))
enrichment(union_group[which(union_group$基因名 %in% sub_group2_down_uniq), ], group, glue('R{sample_name2}_vs_R{control_name2}.down_uniq'))
# enrichment(union_group[which(union_group$基因名 %in% sub_group3_up_uniq), ], group, glue('R{sample_name3}_vs_R{control_name3}.up_uniq'))
# enrichment(union_group[which(union_group$基因名 %in% sub_group3_down_uniq), ], group, glue('R{sample_name3}_vs_R{control_name3}.down_uniq'))


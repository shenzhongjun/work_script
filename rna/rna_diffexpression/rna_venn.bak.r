#!/home/zhouyj/anaconda3/envs/R4/bin/Rscript

library(getopt)

spec <- matrix(c(
  'help', 'h', 0, 'logical', 'Show this help.',
  'groups', 'g', 1, 'character', 'The group file path.',
  'workdir', 'o', 1, 'character', 'The work dir.'
), byrow = T, ncol = 5)
opt = getopt(spec = spec)

# 输入参数检查
if (!is.null(opt$help) ||
  is.null(opt$groups) ||
  is.null(opt$workdir)) {
  cat(getopt(spec = spec, usage = T), '\n')
  quit()
}

# 参数设置
opt$outdir = sub("/$", "", opt$outdir)    # 去掉末尾的'/'防止报错
groups = normalizePath(opt$groups)    # 获取绝对路径
wd = normalizePath(opt$workdir)

# START
library(glue)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

setwd(wd)

groups = read.table(groups, check.names = F, col.names=c('group', 'sample', 'control'))
group = 'group7'
dir.create(glue('{wd}/venn_enrichment/{group}'), showWarnings=FALSE)
setwd(glue('{wd}/venn_enrichment/{group}'))

sub_group1 = groups[which(groups$group == group), ][1, ]
sub_group2 = groups[which(groups$group == group), ][2, ]
# sub_group3 = groups[which(groups$group == group), ][3, ]
sample_name1 = sub_group1$sample
control_name1 = sub_group1$control
sample_name2 = sub_group2$sample
control_name2 = sub_group2$control
# sample_name3 = sub_group3$sample
# control_name3 = sub_group3$control
sub_group1_up = read.table(
  glue('{wd}/{sample_name1}_vs_{control_name1}/DiffExpression/Sample_vs_Control/R{sample_name1}_vs_R{control_name1}.DEG.final.Up.txt'),
  header = T, check.names = F)
sub_group1_down = read.table(
  glue('{wd}/{sample_name1}_vs_{control_name1}/DiffExpression/Sample_vs_Control/R{sample_name1}_vs_R{control_name1}.DEG.final.Down.txt'),
  header = T, check.names = F)
sub_group2_up = read.table(
  glue('{wd}/{sample_name2}_vs_{control_name2}/DiffExpression/Sample_vs_Control/R{sample_name2}_vs_R{control_name2}.DEG.final.Up.txt'),
  header = T, check.names = F)
sub_group2_down = read.table(
  glue('{wd}/{sample_name2}_vs_{control_name2}/DiffExpression/Sample_vs_Control/R{sample_name2}_vs_R{control_name2}.DEG.final.Down.txt'),
  header = T, check.names = F)
# sub_group3_up = read.table(
#   glue('{wd}/{sample_name3}_vs_{control_name3}/DiffExpression/Sample_vs_Control/R{sample_name3}_vs_R{control_name3}.DEG.final.Up.txt'),
#   header = T, check.names = F)
# sub_group3_down = read.table(
#   glue('{wd}/{sample_name3}_vs_{control_name3}/DiffExpression/Sample_vs_Control/R{sample_name3}_vs_R{control_name3}.DEG.final.Down.txt'),
#   header = T, check.names = F)

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

# VENN图
# venn_list_up = list(sub_group1_up$基因名, sub_group2_up$基因名, sub_group3_up$基因名)
# title_up = glue('{group}_UpGenes')
# names(venn_list_up) = c(glue('{sample_name1}_vs_{control_name1}'), glue('{sample_name2}_vs_{control_name2}'), glue('{sample_name3}_vs_{control_name3}'))
# venn.diagram(
#   disable.logging = T,
#   main = title_up,
#   x = venn_list_up,
#   height = 4000,
#   width = 4000,
#   filename = glue('{group}.UpGenes.Venn.png'),
#   fill = c("red", "blue", "green"),
#   alpha = 0.5,
#   cat.default.pos = "text",
# )
#
# venn_list_down = list(sub_group1_down$基因名, sub_group2_down$基因名, sub_group3_down$基因名)
# title_down = glue('{group}_DownGenes')
# names(venn_list_down) = c(glue('{sample_name1}_vs_{control_name1}'), glue('{sample_name2}_vs_{control_name2}'), glue('{sample_name3}_vs_{control_name3}'))
# venn.diagram(
#   disable.logging = T,
#   main = title_down,
#   x = venn_list_down,
#   height = 4000,
#   width = 4000,
#   filename = glue('{group}.DownGenes.Venn.png'),
#   fill = c("red", "blue", "green"),
#   alpha = 0.5,
#   cat.default.pos = "text",
# )
venn_list_up = list(sub_group1_up$基因名, sub_group2_up$基因名)
title_up = glue('{group}_UpGenes')
names(venn_list_up) = c(glue('{sample_name1}_vs_{control_name1}'), glue('{sample_name2}_vs_{control_name2}'))
venn.diagram(
  disable.logging = T,
  main = title_up,
  x = venn_list_up,
  height = 4000,
  width = 4000,
  filename = glue('{group}.UpGenes.Venn.png'),
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.default.pos = "text",
)

venn_list_down = list(sub_group1_down$基因名, sub_group2_down$基因名)
title_down = glue('{group}_DownGenes')
names(venn_list_down) = c(glue('{sample_name1}_vs_{control_name1}'), glue('{sample_name2}_vs_{control_name2}'))
venn.diagram(
  disable.logging = T,
  main = title_down,
  x = venn_list_down,
  height = 4000,
  width = 4000,
  filename = glue('{group}.DownGenes.Venn.png'),
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.default.pos = "text",
)

# 富集分析
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


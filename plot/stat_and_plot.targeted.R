# 2022年3月9日：因生产流程0.7版old_pipline文件改动，修改name变量
library("ggplot2")
library("plyr")
library(argparser)
library(data.table)
library(ggsignif)
library(patchwork)

argv <- arg_parser('')
argv <- add_argument(argv,"--infile", help="", default="/mnt/share05/clinical_project/projects/blood_tumor/test/01.HRD_lixx/03.HRD_standard_product/03.purity_rank/plot/new.result.stat")
argv <- add_argument(argv,"--outdir", help="the output path", default="/mnt/share05/clinical_project/projects/blood_tumor/test/01.HRD_lixx/03.HRD_standard_product/03.purity_rank/plot")
argv <- parse_args(argv)

infile <- argv$infile
out <- argv$outdir

all <- read.table(infile, header=TRUE, sep="\t", check.names=F)
all <- as.data.frame(all)

name <- c("chr", "start", "end", "region", "gene", "dep/len", "average_depth", "cov/len", "coverage_rate")
result <- as.data.frame(matrix(nrow=0, ncol=6))
colnames(result) <- c("sample", "type", "value", "number", "total", "percent")
j <- 1
for (i in 1:nrow(all)){
	prefix <- all[i, "pre"]
	file <- all[i, "file"]
	data <- read.table(file, header=F, sep="\t", check.names=F)
	data <- as.data.frame(data)	
	colnames(data) <- name
	mean <- mean(data$average_depth)
	for (per in c(0.2, 0.5, 1)){
		mean_tmp <- mean * per
		mean_tmp_num <- nrow(subset(data, data$average_depth>mean_tmp))
		mean_tmp_percent <- round(mean_tmp_num/nrow(data) * 100, 1)
		result[j, ] <- c(prefix, paste(per, "mean", sep= "x "), mean_tmp, mean_tmp_num, nrow(data), mean_tmp_percent)
		j <- j+1 }
}

write.table(result, file=paste(out, "result.stat.txt", sep="/"), quote = FALSE, sep = "\t", row.names = FALSE)

result <- read.table(paste(out, "result.stat.txt", sep="/"), header=TRUE, sep="\t", check.names=F)
result <- as.data.frame(result)
type_factor <- factor(result$type, levels=c("0.2x mean", "0.5x mean", "1x mean"))
#barplot(rep(1,n),col=cm.colors(n), border=NA, main="cm.colors")
sample_list <- unique(result$sample)
print(sample_list)
n <- length(sample_list)
col_list <- cm.colors(n)
colors <- colorRampPalette(c("darkcyan", "lightcyan"))(n)
print(col_list)
ggplot(result, aes(x = type, y=percent, fill=sample)) + 
	geom_bar(stat="identity", width = 0.6, position = position_dodge(width=0.8)) + 
	geom_text(aes(label = percent),position=position_dodge(width = 0.8),size = 1.5,vjust = -0.25)+ 
	guides(fill = guide_legend(reverse = T)) +  
	scale_y_continuous(limits=c(0,110), breaks=c(0,20,40,60,80,100), expand = c(0, 0)) +
	scale_fill_manual(values=colors)+
	xlab('')+ylab('Targets covered(%)') +
	theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggsave(paste(out, 'bar.png', sep="/"), width = 12, height = 5)

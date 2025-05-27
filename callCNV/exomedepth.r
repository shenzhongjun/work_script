# 用法：/mnt/share02/zhouyj/software/miniconda3/envs/exomedepth/bin/Rscript /mnt/share02/zhouyj/script/callCNV/exomedepth.r tumorbam normalbam chip.bed 0.01
library(ExomeDepth)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Not enough arguments. Please provide at least 3 arguments.")
}
print(args)

transition.probability = as.numeric(args[4])
fasta = '/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta'
bed <- read.table(args[3], header = FALSE, stringsAsFactors = FALSE)
names(bed) <- c("chromosome", "start", "end", "x", 'name')
bed = bed[, c("chromosome", "start", 'end', "name")]
bed$chromosome <- gsub("chr", "", bed$chromosome)

args <- commandArgs(trailingOnly = TRUE)

bam_files = c(args[1], args[2])

bam_counts <- getBamCounts(bed.frame = bed, bam.files = bam_files, include.chr = TRUE, referenceFasta = fasta)

test_sample <- basename(args[1])
reference_sample <- basename(args[2])
reference_set <- apply(bam_counts[, reference_sample, drop = FALSE],MAR = 1,FUN = sum)
exome_model <- new("ExomeDepth", test = bam_counts[, test_sample], reference = reference_set, formula = "cbind(test, reference) ~ 1")
exome_model <- CallCNVs(x = exome_model, transition.probability = transition.probability,
    chromosome = bam_counts$chromosome, start = bam_counts$start, end = bam_counts$end, name = bam_counts$exon)
cnv_results <- exome_model@CNV.calls
write.table(cnv_results, file = "CNV_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)


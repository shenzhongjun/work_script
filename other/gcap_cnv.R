#!/usr/bin/env Rscript

library(optparse)
option_list = list(
  make_option(c("-f", "--cnvfile"),       type = "character", metavar="character", default=NULL, help="{sample}.genemetrics.addCN.xls file by cnvkit[required]."),
  make_option(c("-s", "--sample"),        type = "character", metavar="character", default=NULL, help="Tumor sample name[required]."),
  make_option(c("-o", "--outdir"),        type = "character", metavar="character", default=getwd(), help="Result output path, default is getwd()."),
  make_option(c("-p", "--purity"),        type = "numeric",   metavar="numeric", default=1, help="Tumor purity, default is 1."),
  make_option(c("-m", "--model"),         type = "character", metavar="character", default="XGB11", help="Trained model name, should be one of XGB11, XGB32, XGB56, default is XGB11."),
  make_option(c("-g", "--genome"),        type = "character", metavar="character", default="hg19", help="Genome build version, default is hg19."),
  make_option(c("-t", "--tightness"),     type = "integer",   metavar="integer", default=1L, help="Control the tightness to be a circular amplicon. If the value is larger, it is more likely a fCNA assigned to 'noncircular' instead of 'circular', default is 1."),
  make_option(c("-c", "--gapCN"),         type = "integer",   metavar="integer", default=3L, help="A gene with copy number above background (ploidy + gapCN in general) would be treated as focal amplicon. Smaller, more amplicons, default is 3."),
  make_option(c("-n", "--onlyOncogenes"), type = "logical",   metavar="logical", default=FALSE, help="Only known oncogenes are kept for circular prediction, default is FALSE.")
)
opt = parse_args(OptionParser(option_list=option_list))

if (is.null(opt$cnvfile) || is.null(opt$sample)) {
  print_help(parser)
  stop("Error: the required argument missing.")
}

cnvfile = opt$cnvfile
sample = opt$sample
outdir = opt$outdir
purity = opt$purity
model = opt$model
genome = opt$genome
tightness = opt$tightness
gapCN = opt$gapCN
onlyOncogenes = opt$onlyOncogenes

cnv = read.table(cnvfile, header=T)
cnv = cnv[, c(2,3,4,9)]
cnv$total_cn = cnv$cn
cnv$minor_cn = NA
cnv$sample=sample
cnv$purity=purity

library(gcap)
gcap.ASCNworkflow(cnv,
                  genome_build=genome,
                  model=model,
                  tightness=tightness,
                  gap_cn=gapCN,
                  only_oncogenes=onlyOncogenes,
                  outdir=outdir,
                  result_file_prefix=sample)

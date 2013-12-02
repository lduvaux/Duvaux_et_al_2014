#!/bin/Rscript
library('multicore')


source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./params.R")
source("./functions.R")

main <- function(argv){

    tab <- read.delim(TARG_INFO, stringsAsFactors = F)

	# 1) stats per gene
	counts <- Count_Genes_targets(tab)
	mat_targets <- counts$mat_targets
	mat_whole_gene <- counts$mat_whole_gene

	# 2) table of targets
	Total <- rowSums(mat_targets)
	mat_targets_f <- cbind(Category=rownames(mat_targets), Total, mat_targets)

	inds <- PrePro_findIndex(ord, rownames(mat_targets))
	mat_targets_f <- mat_targets_f[inds,]
	write.table(mat_targets_f, file=OUT_TARG, sep="\t", quote=F, row.names=F)

	# 3) table of genes
	Total <- rowSums(mat_whole_gene)
	Remaining <- mat_whole_gene[,1]+mat_whole_gene[,2]
	mat_whole_gene_f <- cbind(Category=rownames(mat_whole_gene), Total, Remaining, mat_whole_gene)

	mat_whole_gene_f <- mat_whole_gene_f[inds,]
	write.table(mat_whole_gene_f, file=OUT_GENES, sep="\t", quote=F, row.names=F)

    ######################
    dummy <- numeric()
    outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    save (dummy, file=outFileName)
}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










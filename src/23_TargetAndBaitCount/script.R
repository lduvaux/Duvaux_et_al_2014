#!/bin/Rscript
library('parallel')


source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./params.R")
source("./functions.R")

main <- function(argv){

    ######################
    # 0) set up the table of counts
	count_table <- data.frame(GENE_CATEG, matrix(0, nrow=length(GENE_CATEG), ncol=length(SETS)))
	colnames(count_table)[2:ncol(count_table)] <- SETS
	
    # 1) Initial number of targets and baits (total: 3610 targets, 486 genes (no PMT))
	Cleaning_sum_up <- fill_count_table(count_table, ALL_TARG_INFO, GENE_CATEG, PRO_SEGMENT, BAITS_BE4_CLEAN)
	count_table <- Cleaning_sum_up$Count_table
	NonTrimGenes <- Cleaning_sum_up$NonTrimGenes
    ######################
    write.table(count_table, file=OUTFIL, row.names=F, sep="\t", quote=F)
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    save(NonTrimGenes, file=outFileName)

}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










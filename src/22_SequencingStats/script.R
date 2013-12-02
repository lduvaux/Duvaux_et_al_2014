#!/bin/Rscript
library('multicore')


source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./params.R")
source("./functions.R")

main <- function(argv){

    ######################
    # 1) median of sequencing depth
    tseqdep <- get_SeqDepth(SEQ_DEPTH)
    totmed <- median(tseqdep)
    cat(paste("Median sequencing depth over targets and individuals:", round(totmed)), file=OUTFIL1)
	med_indiv <- round(apply(tseqdep, 2, median), 0)

	# 2) number of reads per individuals
	treadsNbers <- get_readsNber(FASTQC_STATS, BAD_FILE)

	# 3) compute other parameters and set up the final table
	ftab <- set_ftab(treadsNbers, BAD_FCELL, BAD_INDIV, med_indiv)
    
    ######################
    write.table(ftab, file=OUTFIL2, row.names=F, sep="\t", quote=F)
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    dummy <- numeric()
    save(dummy ,file=outFileName)

}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










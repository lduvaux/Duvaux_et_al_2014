#!/bin/Rscript
rm(list=ls())
library('multicore')

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./params.R")
source("./functions.R")

main <- function(argv){

    ######################
    cat("\n")
    print("# 22.1) median of sequencing depth")
    tseqdep <- get_SeqDepth(SEQ_DEPTH)
    totmed <- median(tseqdep)
    cat(paste("Median sequencing depth over targets and individuals:", round(totmed)), file=OUTFIL1)
	med_indiv <- round(apply(tseqdep, 2, median), 0)

    cat("\n")
	print("# 22.2) number of reads per individuals")
	treadsNbers <- get_readsNber(FASTQC_STATS, BAD_FILE)

    cat("\n")
	print("# 22.3) compute other parameters and set up the final table")
	ftab <- as.data.frame(set_ftab(treadsNbers, BAD_FCELL, BAD_INDIV, med_indiv))

    cat("\n")
    print("# 22.4) insert size")
    tab_ins_size <- batch_ins_size(INS_PATH, PATTERN)

    # 5) set up the last table
    ind <- match(ftab$clones, tab_ins_size$sample_name)
    tab_ins_size <- tab_ins_size[ind,]
    ftab <- data.frame(ftab[,1:5], InsertSize=tab_ins_size$insert_size, ftab[,6:ncol(ftab)])
    
    
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










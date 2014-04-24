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
    print(totmed <- median(tseqdep))
    cat(paste("Median sequencing depth over targets and individuals:", round(totmed)), file=OUTFIL1)
	med_indiv <- round(apply(tseqdep, 2, median), 0)

    cat("\n")
	print("# 22.2) number of reads per individuals")
	treadsNbers <- get_readsNber(FASTQC_STATS, BAD_FILE)

    cat("\n")
	print("# 22.3) compute other parameters and set up the final table")
	ftab <- set_ftab(treadsNbers, BAD_FCELL, BAD_INDIV, med_indiv)
	ftab <- ftab[,!colnames(ftab)%in%c("dir","CNVestimation")]
    colnames(ftab) <- c("Sample Name", "Number of reads", "Read length","MTSD", "Encoding", "File1", "File2", "Sequencing Date", "Genepool Flowcell", "Lane", "Genepool Tag")

    cat("\n")
	print("# 22.4) Get capture metrics")
    tab_capture_metrics <- batch_capture_metrics(PATH, PATT_CAPTURE)
    ind <- match(ftab$"Sample Name", rownames(tab_capture_metrics))
    tab_capture_metrics <- tab_capture_metrics[ind,]
    
    cat("\n")
    print("# 22.5) insert size")
    tab_ins_size <- batch_ins_size(PATH, PATT_INSERT)
    ind <- match(ftab$"Sample Name", tab_ins_size$sample_name)
    tab_ins_size <- tab_ins_size[ind,]

    cat("\n")
    print("# 22.6) set up the last table")
    ftab <- data.frame(ftab[,1:4],tab_capture_metrics, MIS=tab_ins_size$insert_size, ftab[,5:ncol(ftab)])
    
    cat("\n")
    print("# 22.7) Draw plot capture metrics")
    pdf(HIST_CAPT)
    layout(matrix(1:4, nrow=2, ncol=2, byrow=T))
    hist(ftab$MTSD, breaks=12, main="Median of target sequencing depth", xlab="Sequencing depth")
    hist(ftab$Enrichment, breaks=12, main="Target enrichment", xlab="Fold enrichment relatively to  \nbackground sequencing depth")
    hist(ftab$Efficiency, breaks=12, main="Capture efficiency", xlab="Proportion of total sequenced bases\nmapping to targets")
    hist(ftab$PTB30X, breaks=12, main="Proportion of bp sequenced at 30X", xlab="Proportion of base pairs with\na sequencing depth of at least 30X")
    dev.off()
    
    ######################
    print("save data")
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










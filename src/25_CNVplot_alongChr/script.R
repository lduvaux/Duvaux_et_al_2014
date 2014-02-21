#!/bin/Rscript
rm(list=ls())
source("../utils/functions.R")
source("../utils/globalCtes.R")
source("params.R")
source("./functions.R")

main <- function(argv){

    cat("\n")
    print(" #### 1) load data and colors")
        # 1.0) load subtargets coordinates
    subtarg <- read.table(SUBTARG, stringsAsFactors=F)
    refBaits <- sapply(subtarg$V5, Prepro_fixName)
    subtarg[,"V5"] <- refBaits
    colnames(subtarg) <- c("Contig", "Start", "End", "Strand", "Name")
    
        # 1.1) load reads count
    load(PREVIOUS_DATA)
        # 1.2) remove bad cytisus
    ind_bad <- PrePro_findIndex(BAD_CYTISUS, colnames(alpha_matrix))
    ind_bad2 <- PrePro_findIndex(c(BAD_CYTISUS, RACE_UNKNOWN), colnames(alpha_matrix))
    alpha_matrix <- alpha_matrix[,-ind_bad]
    print(dim(alpha_matrix))

        # 1.3) set up marker names
    subtargets <- rownames(alpha_matrix)
    targets <- unique(sapply(subtargets, collapse_elements))
    genes0 <- unique(sapply(targets, collapse_elements, what=1:2))
    bad <- grep("PMT", genes0, fixed=T)
    genes <- genes0[-bad]
#~    PMT <- genes0[bad]

    # 2 test with th efirst gene
    gn <- genes[1]
    ind <- grep (paste(gn, "_", sep=""), rownames(alpha_matrix))
    lcnv <- alpha_matrix[ind,]
    ind <- grep (paste(gn, "_", sep=""), subtarg[,"Name"])
    lsubtarg <- subtarg[ind, ]

}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);

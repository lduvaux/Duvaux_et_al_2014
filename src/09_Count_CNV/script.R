#!/bin/Rscript
rm(list=ls())
source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("../utils/randomForest_helperFuns.R")
source("params.R")
source("./functions.R")
require(VennDiagram)


main <- function(argv){
	load(PREVIOUS_DATA)
#~	print(ls())

    cat("\n")
    print(" #### 1) load data and colors")
        # 1.1) reads count
	raw_data_file <- RAW_DATA
	alpha_matrix <- PrePro_roundToZeroFive(alpha_matrix)
        # 1.2) remove bad cytisus
    ind_bad <- PrePro_findIndex(BAD_CYTISUS, colnames(alpha_matrix))
    ind_bad2 <- PrePro_findIndex(c(BAD_CYTISUS, RACE_UNKNOWN), colnames(alpha_matrix))
    alpha_matrix <- alpha_matrix[,-ind_bad]

    clones <- colnames(alpha_matrix)
    races <- unique(sapply(clones, function(x) unlist(strsplit(x, "_"))[1]))
    subtargets <- rownames(alpha_matrix)
    targets <- unique(sapply(subtargets, collapse_elements))
    genes0 <- unique(sapply(targets, collapse_elements, what=1:2))
    bad <- grep("PMT", genes0, fixed=T)
    genes <- genes0[-bad]
    PMT <- genes0[bad]


    cat("\n")
    print(" #### 2) count global proportion of deletion/insertions/both")
        # 2.1) subtargets
    CNV_count_bait <- t(apply(alpha_matrix, 1, caract_bait))
    colnames(CNV_count_bait) <- c("Deletion", "Insertion", "Both")
    toto <- apply(CNV_count_bait, 2, sum)
    list4venn_bait <- list(Deletion = 1:(toto[1]+toto[3]), Insertion = (toto[1]+1):sum(toto))

        # 2.2) targets
    CNV_count_targ <- t(sapply(targets, get_caract_targ, subtargets, alpha_matrix))
    colnames(CNV_count_targ) <- c("Deletion", "Insertion", "Both")
    toto_targ <- apply(CNV_count_targ, 2, sum)
    list4venn_targ <- list(Deletion = 1:(toto_targ[1]+toto_targ[3]), Insertion = (toto_targ[1]+1):sum(toto_targ))

        # 2.3) genes
    CNV_count_gn <- t(sapply(genes, get_caract_targ, subtargets, alpha_matrix))
    colnames(CNV_count_gn) <- c("Deletion", "Insertion", "Both")
    toto_gn <- apply(CNV_count_gn, 2, sum)
    list4venn_gn <- list(Deletion = 1:(toto_gn[1]+toto_gn[3]), Insertion = (toto_gn[1]+1):sum(toto_gn))

        # 2.2) promoters
    CNV_count_pmt <- t(sapply(PMT, get_caract_targ, subtargets, alpha_matrix))
    colnames(CNV_count_pmt) <- c("Deletion", "Insertion", "Both")
    toto_pmt <- apply(CNV_count_pmt, 2, sum)
    list4venn_pmt <- list(Deletion = 1:(toto_pmt[1]+toto_pmt[3]), Insertion = (toto_pmt[1]+1):sum(toto_pmt))


    cat("\n")
    print(" #### 3) count global proportion of deletion/insertions/both per race")
        # 3.1) subtargets
    CNV_count_race <- mclapply(1:nrow(alpha_matrix), function(x) sapply(races, caract_bait_race, clones, alpha_matrix[x,], x))
    names(CNV_count_race) <- subtargets

    # x) save results
	outFileName <- argv[1]
    ver(sprintf("Saving *DUMMY* data to %s",outFileName))
    dummy <- numeric()
    save(dummy,file=outFileName)
}





argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










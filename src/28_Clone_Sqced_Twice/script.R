#!/bin/Rscript
    # load global ressources
source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

    # load local ressources
source("./params.R")
source("./functions.R")

main <- function(argv){
	load(PREVIOUS_DATA)
	set.seed(0)

    cat("\n")
    print(" #### 1) load data and colors")
        # 1.1) reads count
	alpha_matrix <- PrePro_roundToZeroFive(alpha_matrix)
        # 1.2) remove bad cytisus
    ind_bad <- PrePro_findIndex(BAD_CYTISUS, colnames(alpha_matrix))
    ind_bad2 <- PrePro_findIndex(c(BAD_CYTISUS, RACE_UNKNOWN), colnames(alpha_matrix))
    alpha_matrix <- alpha_matrix[,-ind_bad]
    print(dim(alpha_matrix))

    clones <- colnames(alpha_matrix)
    races <- unique(sapply(clones, function(x) unlist(strsplit(x, "_"))[1]))
    subtargets <- rownames(alpha_matrix)
    targets <- unique(sapply(rownames(alpha_matrix), collapse_elements))
    
    cat("\n")
    print("##### 2) compare results for the clone sequenced twice")
    mat_twice <- alpha_matrix[,CLONE_TWICE]
    vec_diff <- mat_twice[,1]-mat_twice[,2]
    ind_diff <- which(vec_diff!=0)
    mat_diff <- mat_twice[ind_diff,]
    targets_diff <- unique(sapply(rownames(mat_diff), collapse_elements))

    cat("\n")
    print("##### 3) number of False Positive deletions and insertions")
    print("   ## 3.1) percentage of FP for subtargets")
    print(pc_subtarg <- nrow(mat_diff)*100/nrow(mat_twice))
    print("   ## 3.2) number of FP insertion and deletions for subtargets")
    sum_row_subtarg <- apply(mat_diff, 1, sum)
    test_indel_subt <-ifelse(sum_row_subtarg>2.5, "NA", ifelse(sum_row_subtarg<=1.5, "Del", "Ins"))
    print(FP_subtarg <- table(test_indel_subt))

    print("   ## 3.3) number of FP for targets")
    print(pc_targ <- length(targets_diff)*100/length(targets))
    print("   ## 3.2) number of FP insertion and deletions for subtargets")
    mat_diff2 <- mat_diff[TO_KEEP,]
    sum_row_targ <- apply(mat_diff2, 1, sum)
    test_indel_targ <- ifelse(sum_row_targ>2.5, "NA", ifelse(sum_row_targ<=1.5, "Del", "Ins"))
    print(FP_targ <- table(test_indel_targ))


    cat("\n")
    print(" #### save results")
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    dummy <- numeric()
    save.image(file=outFileName)

}

cat("dummy logfile\n")


argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










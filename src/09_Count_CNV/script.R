#!/bin/Rscript

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("../utils/randomForest_helperFuns.R")
source("params.R")
source("./functions.R")

main <- function(argv){
	load(PREVIOUS_DATA)
#~	print(ls())

    # 1) load data and colors
        # 1.1) reads count
	raw_data_file <- RAW_DATA
	alpha_matrix <- PrePro_roundToZeroFive(alpha_matrix)
        # 1.2) remove bad cytisus
    ind_bad <- PrePro_findIndex(BAD_CYTISUS, colnames(alpha_matrix))
    alpha_matrix <- alpha_matrix[,-ind_bad]

    # 2) count CNV
        # 2.1) globally
    ww <- as.vector(as.vector(alpha_matrix))
    table(ww)
    

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










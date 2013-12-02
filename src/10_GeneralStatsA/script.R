#!/bin/Rscript
source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")
source("../utils/randomForest_helperFuns.R")
source("./params.R")
source("./functions.R")

main <- function(argv){
	
	load(PREVIOUS_DATA)
	alpha_matrix_raw <- PrePro_roundToZeroFive(alpha_matrix)
	IndivRace_raw <- PrePro_fetchRacesPrior(PRIOR_ASSIGN)
	races_raw <- unique(IndivRace_raw)
	ind_order <- PrePro_findIndex(colnames(alpha_matrix_raw), names(IndivRace_raw))
	IndivRace_raw <- IndivRace_raw[ind_order]
	bad_indiv_ind <- which(is.na(IndivRace_raw)|IndivRace_raw=="hybrid")
	IndivRace <- IndivRace_raw[-bad_indiv_ind]
	alpha_matrix <- alpha_matrix_raw[,-bad_indiv_ind]
	IndivRace_raw <- IndivRace_raw[ind_order]
	ListRaces <- unique(IndivRace)
	indivperRace <- table(IndivRace)
	print(ls())
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    dummy <- numeric()
    save(dummy,file=outFileName)

	
}





argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










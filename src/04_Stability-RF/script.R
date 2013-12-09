#!/bin/Rscript

library(randomForest)
source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")
source("../utils/randomForest_helperFuns.R")
source("params.R")
source("./functions.R")

main <- function(argv){

    load(PREVIOUS_DATA)
	set.seed(0)
	# 0) set up the data set
		# 0.1) the matrix

    data_mat <- alpha_matrix
    rounded_ratio_mat <- PrePro_roundToZeroFive(data_mat)
    
	inval_rows <- PrePro_findMeaninglessRows(rounded_ratio_mat) # loci with no variation
	rounded_ratio_mat <- rounded_ratio_mat[!inval_rows,]
	tdata_mat <- t(rounded_ratio_mat)

		# 0.2) races, individuals and libraries
    races <- PrePro_fetchRaces(RAW_DATA, CTRL_GUYS , BAD_GUYS)
    races[races=="Medicago_ctrl"] <- "Medicago"
	races_uniq <- PrePro_fetchUniqRaces(races)
	
	indiv <- colnames(data_mat)
	lib_names <- PrePro_LibName(INDIV_DETAILS,indiv)
	
		# 0.3) removing bad cytisus from the matrix
	ind <- PrePro_findIndex(BAD_CYTISUS, indiv)
	tdata_mat <- tdata_mat[-ind,]
	lib_names <- lib_names[-ind]
	races <- races[-ind]
	indiv <- rownames(tdata_mat)
		# 0.4) set up xtest, ytest, x and y
	ind <- PrePro_findIndex(RACE_UNKNOWN, indiv)
	x <- tdata_mat[-ind,]
	xtest <- tdata_mat[ind,]
	y0 <- as.factor(races)
	
	y <- y0[-ind]
	ytest <- y0[ind]

		# 1) check results stability
			# 1.1) stability ~ number of trees (growing 100 RFs of 100 trees each)
	print(paste("Grow ", N_RF, " supervised RFs of ", NTREES, " trees", sep=""))
	print(paste("Matrix dimensions for supervised RFs:", paste(dim(tdata_mat), collapse=" ")))
	print(system.time(lrf_superv <- Imptce_Stability_superv(N_RF, NTREES, tdata_mat, races, SAMP_SIZE1)))
	
	print(paste("Grow ", N_RF, " unsupervised Rfs of ", NTREES, " trees", sep=""))
	print(paste("Matrix dimensions for unsupervised RFs:", paste(dim(x), collapse=" ")))
	print(system.time(lrf_unsuperv <- Imptce_Stability_unsuperv(N_RF, NTREES, x, y, xtest, ytest, SAMP_SIZE2)))
	
	
	print("Combine small supervised Rfs in Bigger supervised RFs")
	print(system.time(gradNtree_lrf_superv <- set_combis_nRF(lrf_superv, lg=100, by=10)))
	print("Combine small unsupervised Rfs in Bigger unsupervised RFs")
	print(system.time(gradNtree_lrf_unsuperv <- set_combis_nRF(lrf_unsuperv, lg=100, by=10)))


    plot_10_bestvar(gradNtree_lrf_superv, lrf_superv, ysup=0.02, pdf_name=PDF_NAME1)
    plot_10_bestvar(gradNtree_lrf_unsuperv, lrf_unsuperv, ysup=0.02, pdf_name=PDF_NAME2)
			# 1.2) stability ~ reference pools
	best20_superv <- sort_VarPerImp(gradNtree_lrf_superv)[1:20]
	best20_unsuperv <- sort_VarPerImp(gradNtree_lrf_unsuperv)[1:20]

	print(names(best20_superv)[names(best20_superv)%in%names(best20_unsuperv)])
	print(names(best20_superv)[!names(best20_superv)%in%names(best20_unsuperv)])
	print(names(best20_unsuperv)[!names(best20_unsuperv)%in%names(best20_superv)])
	
	
	
    ######################
    outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    print(head(y))
    save(y, indiv, x, xtest, ytest,SAMP_SIZE2 , gradNtree_lrf_unsuperv,file=outFileName)
#~     save (data_list, file=outFileName)
}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










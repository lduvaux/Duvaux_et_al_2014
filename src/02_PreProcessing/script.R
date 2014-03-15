#!/bin/Rscript
source("../utils/functions.R")
source("../utils/globalCtes.R")
source("../utils/getter_functions.R")

source("./transformDataByFit.R")
source("./params.R")
source("./functions.R")



##test
main <- function(argv){
    
    raw_data_file <- PREVIOUS_DATA

    set.seed(0)
    
    data_mat <- PrePro_fetchRawDataMat(raw_data_file, BAD_GUYS); dbgPrint(dim(data_mat))
    lib_name <- PrePro_LibName(indiv_details=INDIV_DETAILS,indiv=colnames(data_mat))

#~    races <- PrePro_fetchRaces(raw_data_file, CTRL_GUYS, BAD_GUYS)
#~    races_uniq <- unique(races)
    
    col_lib <- PrePro_LibColo(lib_name, col_lib_raw)
#~    col_races <- PrePro_fetchColours(races, races_uniq, race_colo_raw)
    
    control_vec <- PrePro_fetchRawCtrl(raw_data_file); dbgPrint(length(control_vec))
    quantiles_ctrl <- PrePro_quantiles_ctrl(control_vec, THRES)

    PLOT_MAINS <- list(mains = colnames(data_mat), col.mains=col_lib)
     
    data_list <- PrePro_prePrecessing(control_vec, data_mat, polyn_deg = POLY_DEG, use_intercept = USE_INTERCEPT, weight_error = TRUE, plot_mains=PLOT_MAINS, thresholds = quantiles_ctrl$quantiles)

    ######################
    outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    save (data_list, file=outFileName)
}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);

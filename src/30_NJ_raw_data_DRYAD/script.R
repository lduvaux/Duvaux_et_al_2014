#!/bin/Rscript
rm(list=ls())
library(ape)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("../utils/randomForest_helperFuns.R")
source("params.R")
source("./functions.R")

main <- function(argv){

    # 1load data
	load(PREVIOUS_DATA)
    set.seed(0)
#~	print(ls())
	inter <- ifelse (fit$use_intercept, "Intercept", "NoIntercept")
	raw_data_file <- RAW_DATA
    col_lanes_un <- col_lib_raw
    col_races_un <- race_colo_raw
    
	# 2) setup data and colors
	alpha_matrix <- PrePro_roundToZeroFive(alpha_matrix)
    
	indiv <- colnames(alpha_matrix)
	lanes <- PrePro_LibName(INDIV_DETAILS,indiv)   # lane info
    lanes_un <- unique(lanes)

    pool_info0 <- read.delim(POOL_INFO, stringsAsFactors = F)  # pool info
    pool_info <- PrePro_Clone2Pools(indiv, lanes, pool_info0)
    lty_pools <- PrePro_phylty (lanes_un, pool_info)    # probably slow! for loop behind
    lanepool_un <- PrePro_setLanePool(lanes_un, pool_info)

    races <- PrePro_fetchRaces(raw_data_file, CTRL_GUYS, BAD_GUYS)
    races_un <- unique(races)

    col_races <- PrePro_fetchColours(races, races_un, col_races_un)
    col_lanes <- PrePro_LibColo(lanes, col_lanes_un)
    
    subtarg <- rownames(alpha_matrix)
    subtarg_categ <- PrePro_BaitCateg(subtarg)
    gene_categ <- PreProNJ_GetGeneCateg(subtarg_categ)
    
	# 3) build figures
    myroot <- which(colnames(alpha_matrix)==ROOT)
	outpdf <- paste(PDF_NAME)
	pdf(file=outpdf, width=11, height=14)
	lapply(1:length(gene_categ), processing_NJ, gene_categ, alpha_matrix, METHOD, myroot,
        PLOT_LANES, col_races, col_lanes, lanes_un, races_un, col_lanes_un, col_races_un,
        PLOT_POOLS, lty_pools, lanepool_un,
        subtarg, LGD_POS)
	dev.off()

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










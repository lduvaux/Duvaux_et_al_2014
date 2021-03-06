#!/bin/Rscript
library(ape)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("../utils/randomForest_helperFuns.R")
source("params.R")
source("./functions.R")

main <- function(argv){
	load(PREVIOUS_DATA)
    set.seed(0)
	print(ls())
	inter <- ifelse (fit$use_intercept, "Intercept", "NoIntercept")
	raw_data_file <- RAW_DATA
	
	# 1) setup data and colors
	alpha_matrix <- PrePro_roundToZeroFive(alpha_matrix)
	
	indiv <- colnames(alpha_matrix)
	lib_name <- PrePro_LibName(INDIV_DETAILS,indiv)
    libraries <- unique(lib_name)
    
    races <- PrePro_fetchRaces(raw_data_file, CTRL_GUYS, BAD_GUYS)
    races_uniq <- unique(races)

    col_races <- PrePro_fetchColours(races, races_uniq, race_colo_raw)
    col_lib <- PrePro_LibColo(lib_name, col_lib_raw)
    
    baits <- rownames(alpha_matrix)
    All_bait_categ <- PrePro_BaitCateg(baits)
    gene_categ <- PreProNJ_GetGeneCateg(All_bait_categ)
    
	# 2) build figures
    root <- which(colnames(alpha_matrix)==ROOT)
	outpdf <- paste(PREFIX4PDF_NAME, fit$polyn_deg, "_", inter, "_", METHOD, ".pdf", sep="")
	pdf(file=outpdf, width=11, height=14)
	lapply(1:length(gene_categ), processing_NJ, gene_categ, alpha_matrix, METHOD, root, LANE, col_races, col_lib, libraries, races_uniq, col_lib_raw, race_colo_raw, baits)
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










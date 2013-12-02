#!/bin/Rscript
library('multicore')


source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./Process_glm.R")
source("./params.R")
source("./functions.R")

main <- function(argv){

    load(PREVIOUS_DATA)
	set.seed(0)
	
    ######################
    
    Pre_GLMtab <- setGLMtable(ListGenes, BaitsGeneNames, alpha_matrix, Genes_Info, CATEG_FOR_GLM)
    ######################
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










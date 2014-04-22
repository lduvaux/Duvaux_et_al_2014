#!/bin/Rscript
rm(list=ls())
source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./params.R")
source("./functions.R")

main <- function(argv){

    load(PREVIOUS_DATA)
	set.seed(0)

	print("##### 1) prepare required information for GLM")
	# 1) set up the matrix of relevant individuals
	alpha_matrix_raw <- PrePro_roundToZeroFive(alpha_matrix)
	alpha_matrix <- Trim_Alpha_mat(PRIOR_ASSIGN, alpha_matrix_raw) #Y (of course)

	# 2) fetch the table of gene info + add the field "GeneAlias"
	Genes_Info <- read.delim(INFO_TARGENE_FILE, stringsAsFactors=F)
	Genes_Info <- cbind(Genes_Info,GeneAlias=sapply(Genes_Info$NewTargetName, Pro_geneName))	# add gene name under alias

	# 3) prepare list baits, exons and genes names
	BaitsGeneNames <- sapply(rownames(alpha_matrix), Pro_geneName) #Y
	ListGenes <- unique(BaitsGeneNames) #Y

	# 4) set up list of complete genes
	load(COMPGENES_DATA)


    ######################
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    
    save(alpha_matrix, Genes_Info, BaitsGeneNames, ListGenes, NonTrimGenes, file=outFileName)

}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










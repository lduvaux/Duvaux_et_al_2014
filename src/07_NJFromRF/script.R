#!/bin/Rscript
rm(list=ls())
library(ape)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

#~ source("../utils/randomForest_helperFuns.R")
source("params.R")
source("./functions.R")


main <- function(argv){
	load(PREVIOUS_DATA)
	set.seed(0)
    mat_prox <- contig_rf$proximity
	dim(mat_prox)

    mat_prox_test <- contig_rf$test$proximity
	new <- colnames(mat_prox_test)%in%rownames(mat_prox_test)
	nindiv <- ncol(mat_prox_test)
	fullproxi <- matrix(NA, nrow=nindiv,ncol=nindiv)

	fin_new <- nrow(mat_prox_test)
	deb <- nrow(mat_prox_test)+1
	fin <- deb+ncol(mat_prox)-1

	fullproxi[1:fin_new,] <- mat_prox_test
	fullproxi[,1:fin_new] <- t(mat_prox_test)
	fullproxi[deb:fin, deb:fin] <- mat_prox
	colnames(fullproxi) <- colnames(mat_prox_test)
	rownames(fullproxi) <- colnames(mat_prox_test)

	distance_matrix <- 1-fullproxi
	grep("Lathyrus", colnames(mat_prox_test))
	grep("Lathyrus", colnames(mat_prox_test), value=T)

    races <- PrePro_fetchRaces(RAW_DATA, CTRL_GUYS, BAD_GUYS)
	races_uniq <- unique(races)	
	col_races <- PrePro_fetchColours(races, races_uniq, race_colo_raw)
	good <- names(col_races)%in%colnames(distance_matrix)
	col_races  <- col_races[good]
	ind <- PrePro_findIndex(colnames(distance_matrix), names(col_races))

	
	col_races <- col_races[ind]

    root <- which(colnames(distance_matrix)==ROOT)
	tree <- root(nj(distance_matrix), root)
	tree_edges <- get_tree_edges(tree)
	edge_colors <- apply_colors_2_edges(tree, tree_edges, lane=F, col_races)


	pdf(PDF_NAME, height=12.9, width=9)
    whole_fig(tree, col_races, edge_colors, races_uniq, race_colo_raw, contig_rf, leg_pos="bottomleft")
	dev.off()

	jpeg(JPG_NAME, height=688*2, width=480*2, quality=100, res=72*2)
    whole_fig(tree, col_races, edge_colors, races_uniq, race_colo_raw, contig_rf, leg_pos="bottomleft", cex_leaves=0.5)
	dev.off()

	outFileName <- argv[1]
    ver(sprintf("Saving *DUMMY* data to %s",outFileName))
    dummy <- numeric()
    save(dummy, file=outFileName)

}


argv <- commandArgs(TRUE)[1]
print(argv)
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);



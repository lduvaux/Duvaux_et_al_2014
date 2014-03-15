#!/bin/Rscript
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
	dim(contig_rf$proximity)
	new <- colnames(contig_rf$test$proximity)%in%rownames(contig_rf$test$proximity)

	nindiv <- ncol(contig_rf$test$proximity)
	fullproxi <- matrix(NA, nrow=nindiv,ncol=nindiv)

	fin_new <- nrow(contig_rf$test$proximity)
	deb <- nrow(contig_rf$test$proximity)+1
	fin <- deb+ncol(contig_rf$proximity)-1

	fullproxi[1:fin_new,] <- contig_rf$test$proximity
	fullproxi[,1:fin_new] <- t(contig_rf$test$proximity)
	fullproxi[deb:fin, deb:fin] <- contig_rf$proximity
	colnames(fullproxi) <- colnames(contig_rf$test$proximity)
	rownames(fullproxi) <- colnames(contig_rf$test$proximity)

	distance_matrix <- 1-fullproxi
	grep("Lathyrus", colnames(contig_rf$test$proximity))
	grep("Lathyrus", colnames(contig_rf$test$proximity), value=T)


	races <- PrePro_fetchRaces(RAW_DATA, CTRL_GUYS, BAD_GUYS)
	races_uniq <- unique(races)	

	col_races <- PrePro_fetchColours(races, races_uniq, race_colo_raw)
	good <- names(col_races)%in%colnames(distance_matrix)
	col_races  <- col_races[good]
	ind <- PrePro_findIndex(colnames(distance_matrix), names(col_races))

	
	col_races <- col_races[ind]
#identical(names(col_races),colnames(distance_matrix))

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



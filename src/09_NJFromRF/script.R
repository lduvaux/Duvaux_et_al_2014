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


	tree <- root(nj(distance_matrix), 31)	# Lathyrus
	#tree <- root(nj(distance_matrix), 21)	# cytisus
	tree_edges <- get_tree_edges(tree)
	edge_colors <- apply_colors_2_edges(tree, tree_edges, lane=F, col_races)


	pdf("./NJ_from_RF.pdf", height=12.9, width=9)
	par(mar=c(3, 2, 3, 0))
	plot(tree, show.tip.label = T, tip.color=col_races, edge.color=edge_colors, edge.width=1.5, root.edge=F, cex = 0.75, font=4)

	legend("bottomleft", legend=races_uniq, col = race_colo_raw, lwd=2, bg = 'gray92', cex=1)

	title(paste("NJ tree based on random forest proximity matrix\n(", nrow(contig_rf$importance), " informative genes)", sep=""), cex.main=1.5, font=4)
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










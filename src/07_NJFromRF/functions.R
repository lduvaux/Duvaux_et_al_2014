
PreProNJ_GetGeneCateg <-  function(All_bait_categ){
	cate <- unique(All_bait_categ)
	cate <- c("All", cate)
	return(cate)
}


define_alpha_matrix <- function(alpha_matrix, gene_categ, i_categ, baits)
{
	categ <- gene_categ[i_categ]
	if (i_categ==1)
		talpha <- alpha_matrix
	else {
		ind <- grep(paste("^", categ, "_", sep=""), baits)
		talpha <- alpha_matrix[ind,]
	}
	return(talpha)
}

compute_ludo_dtce <- function(talpha)
{
	nind <- ncol(talpha)
	nb_baits=nrow(talpha)
	
	vdis <- numeric(sum(1:(nind-1)))
	forbid <- NULL
	inc <- 1
	for (ite in 1:(nind-1))
	{
		forbid <- c(forbid, ite)
		vec <- talpha[,ite]
		tab_temp <- talpha[,-forbid]
		if (is.matrix(tab_temp)) 
			dis <- sapply(1:ncol(tab_temp), function(x) nb_baits-sum(vec==tab_temp[,x])) 
		else 
			dis <- nb_baits-sum(vec==tab_temp)
		vdis[inc:(inc+nind-length(forbid)-1)] <- dis
		inc <- inc+nind-length(forbid)
	}
	return(vdis)
}

put_vector_into_matrix <- function(nind, vdis, indiv_names)
{
	M <- matrix(0, nind, nind)
	M[lower.tri(M)] <- vdis
	M <- t(M)
	M[lower.tri(M)] <- vdis
	dimnames(M) <- list(indiv_names)
	return(M)
}

compute_dstce <- function(talpha, method)
{
	if (method=="ludo") {
		vdis <- compute_ludo_dtce(talpha)
		M <- put_vector_into_matrix(ncol(talpha), vdis, colnames(talpha))
	}
	else if (method=="manhattan")
		M <- as.matrix(dist(t(talpha), method="manhattan"))
	else if ("wrong distance method")
	return(M)
}

get_tree <- function(distance_matrix, racine)
{
	tree <- nj(distance_matrix)
	tree <- ladderize(root (tree, racine), right=T)
	return(tree)
}

get_tree_edges <- function(tree)
# select terminal branches
{
	edges <-  which(tree$edge[,2] <= Ntip(tree))
	edge_nb <- tree$edge[edges,2]
	tree_edges <- list(edges=edges, edge_nb=edge_nb)
	return(tree_edges)
}

apply_colors_2_edges <- function(tree, tree_edges, lane, col_races, col_lib)
{
	edge_colors <- rep("black", length(tree$edge[,1]))
	if (lane) 
		edge_colors[tree_edges$edges] <- col_lib[tree_edges$edge_nb] 
	else
		edge_colors[tree_edges$edges] <- col_races[tree_edges$edge_nb]
	return(edge_colors)
}

plot_tree <- function(tree, col_races, edge_colors, lane, libraries, races_uniq, col_lib_raw, race_colo_raw, categ, nbaits)
{
	plot(tree, show.tip.label = T, tip.color=col_races, edge.color=edge_colors, edge.width=1.5, root.edge=F, cex = 0.7, font=4)

	if (lane) 
		legend("topright", legend=c(libraries, races_uniq), col = col_lib_raw, lty = c(rep(1, length(libraries)), rep(-1, length(race_colo_raw))), lwd=2, bg = 'gray92', cex=0.8, text.col=c(rep("black", length(libraries)), race_colo_raw))
	else 
		legend("topright", legend=races_uniq, col = race_colo_raw, lwd=2, bg = 'gray92', cex=0.8)
	title(paste("NJ tree of CNV distribution across aphid individuals\n", categ, " genes (", nbaits, " baits)", sep=""), cex.main=1.5, font=4)
}

whole_fig <- function(tree, col_races, edge_colors, races_uniq, race_colo_raw, contig_rf, leg_pos="bottomleft", cex_leaves=0.75, leg_cex=1)
{
	par(mar=c(3, 2, 3, 0))
	plot(tree, show.tip.label = T, tip.color=col_races, edge.color=edge_colors, edge.width=1.5, root.edge=F, cex = cex_leaves, font=4)

	legend(leg_pos, legend=races_uniq, col = race_colo_raw, lwd=2, bg = 'gray92', cex=leg_cex, text.font=3)

	title(paste("NJ tree based on random forest proximity matrix\n(", nrow(contig_rf$importance), " informative genes)", sep=""), cex.main=1.5, font=4)
}




















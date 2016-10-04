#~ source("./02_ProcessScripts/transformDataByFit.R")
#~ source("./00_Accessory_Scripts/getter_functions.R")
PreProNJ_GetGeneCateg <-  function(All_bait_categ){
	cate <- unique(All_bait_categ)
	cate <- c("All", cate)
	return(cate)
}

PrePro_Clone2Pools <- function(clones, lane_info, tab_pool_info){
# collect pool information for clones
    matri <- sapply(clones, function(x) unlist(strsplit(x, "_"))[2])
    
    matri_EBI <- sapply(tab_pool_info[,"sample_alias"], function(x) unlist(strsplit(x, "_"))[3])
    bad <- which(!matri_EBI%in%matri)
    matri_EBI <- matri_EBI[-bad]
    tab_pool_info <- tab_pool_info[-bad,]
    
    ind <- match(matri, matri_EBI)
    tab_pool_info <- tab_pool_info[ind,]

    res <- data.frame(Indiv=clones, Lane=lane_info, Pool=tab_pool_info$library_name, file1=tab_pool_info$forward_file_name, file2=tab_pool_info$reverse_file_name, stringsAsFactors = F)
    return(res)
}

PrePro_lanelty <- function(lane, tab_pool_info){

    tab <- subset(tab_pool_info, Lane==lane)
    vpools <- tab$Pool
    uniq_pool <- unique(vpools)
#~    print(length(uniq_pool))
    if (length(uniq_pool)!=2)
        stop(paste("For lane ", lane, ", the number of pool is not 2! Check the data\n", sep=""))
    lty <- ifelse(vpools==uniq_pool[1], 1, 2)
    return(lty)
}

PrePro_phylty <- function(all_lanes, tab_pool_inf){
    all_lty <- numeric(nrow(tab_pool_inf))
    for (lan in all_lanes){
#~        print(lan)
        vlty <- PrePro_lanelty (lan, tab_pool_inf)
        ind <- which(tab_pool_inf$Lane==lan)
        all_lty [ind] <- vlty
    }
    return(all_lty)
}

PrePro_setLanePool <- function(lanes_uni, tab_pool_inf){
    v <- NULL
    for (i in lanes_uni){
        tab <- subset(tab_pool_inf, Lane==i)
        pools <- unique(tab$Pool)
        aa <- paste(i, pools, sep=".")
        v <- c(v, aa)
    }
    return(v)
}

define_alph_mat <- function(alph_mat, gn_cat, i_cat, baits)
{
	categ <- gn_cat[i_cat]
	if (i_cat==1)
		talpha <- alph_mat
	else {
		ind <- grep(paste("^", categ, "_", sep=""), baits)
		talpha <- alph_mat[ind,]
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

get_tree <- function(distance_matrix, my_root)
{
	tree <- nj(distance_matrix)
	tree <- ladderize(root (tree, my_root), right=T)
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

processing_NJ <- function(
    i_cat, gn_cat, alph_mat, method, root,
    plot_lanes, colo_races, colo_lanes, lanes_uni, races_uni, col_lanes_uni, col_races_uni,
    plot_pools, lity_pools, lanepool_uni,
    sub_targ, lgd_pos)
{
	categ <- gn_cat[i_cat]
    lgdx <- lgd_pos[i_cat]
	talpha <- define_alph_mat(alph_mat, gn_cat, i_cat, sub_targ)
	nsubtarg <- nrow(talpha)
	dtce_mat <-  compute_dstce(talpha, method)
	tree <- get_tree(dtce_mat, root)
	
	tree_edges <- get_tree_edges(tree)
	if (plot_lanes)
        edge_colors <- apply2edges("black", tree, tree_edges, colo_lanes)
    else
        edge_colors <- apply2edges("black", tree, tree_edges, colo_races)
    if (plot_pools)
        edge_ltys <- apply2edges(1, tree, tree_edges, lity_pools)
    else
        edge_ltys <-  1

    plot_tree(tr=tree, tip_col=colo_races, edge_col=edge_colors, edge_lty=edge_ltys, plt_lanes=plot_lanes, lanepool_un=lanepool_uni, races_un=races_uni, col_lane=col_lanes_uni, race_colo_raw=col_races_uni, categ=categ, nbaits=nsubtarg, legx=lgdx)
	
	print(categ)
	print(dim(talpha))

}

plot_tree <- function(tr, tip_col, edge_col, edge_lty, plt_lanes, lanepool_un, races_un, col_lanes_un, race_colo_raw, categ, nbaits, legx)
{
	plot(tr, tip.color=tip_col, edge.color=edge_col, edge.lty=edge_lty, edge.width=1.5, show.tip.label = T, root.edge=F, cex = 0.68, font=4)

	if (plt_lanes) 
		legend(legx, legend=c(lanepool_un, races_un), col = rep(col_lanes_un, each=2), lty = c(rep(1:2, length(col_lanes_un)), rep(-1, length(race_colo_raw))), lwd=2, bg = 'gray92', cex=0.8, text.col=c(rep("black", length(lanepool_un)), race_colo_raw))
	else 
		legend("topright", legend=races_un, col = race_colo_raw, lwd=2, bg = 'gray92', cex=0.8)
	title(paste("NJ tree of CNV distribution across aphid individuals\n", categ, " genes (", nbaits, " subtargets)", sep=""), cex.main=1.5, font=4)
}

apply2edges <- function(common, tree, tree_edges, feature_races)
{
    edge_feature <- rep(common, length(tree$edge[,1]))
    edge_feature[tree_edges$edges] <- feature_races[tree_edges$edge_nb]
    return(edge_feature)
}



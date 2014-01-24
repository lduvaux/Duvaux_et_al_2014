#~ source("./02_ProcessScripts/transformDataByFit.R")
#~ source("./00_Accessory_Scripts/getter_functions.R")
#~ library(parallel)

################ new functions for rank test ################
get_elements <- function(x, sep="_", what=c(1)){
    return(strsplit(x, sep)[[1]][what])
}

collapse_elements <- function(nom, sep="_", what=1:3, colla="_"){
    res <- paste(get_elements(nom, sep, what), collapse=colla)
    return(res)
}

count_categ <- function(){
    load(PREVIOUS_DATA2)
    noms <- rownames(alpha_matrix)
    print(system.time(count_categ <- table(sapply(noms, get_elements))))
    return(count_categ)
}

get_contig <- function(target, info_TargGene){
    return(info_TargGene$contigV2[info_TargGene$NewTargetName==target])
}

get_P_same_contig <- function(noms1, noms2, info_TargGene_fil, verbose=T){


    v1 <- sapply(noms1, collapse_elements, sep="_", what=1:3, colla="_") 
    v2 <- sapply(noms2, collapse_elements, sep="_", what=1:3, colla="_") 
    info_TargGene <- read.delim(info_TargGene_fil)
    vcontigs1 <- sapply(v1, get_contig, info_TargGene)
    vcontigs2 <- sapply(v2, get_contig, info_TargGene)

    res0 <- outer(vcontigs1, vcontigs2, "==")
    res <- table(res0)
    if (identical(vcontigs1, vcontigs2)) {
        if (verbose) print('Symetrical matrix')
        res[2] <- res[2]-length(vcontigs1) # remove diagonal
        res <- res/2    # discard one half of the matrix (since it is symetric)
    } else if (verbose) print('Asymetrical matrix')
    P <- res[2]/sum(res)

    return(P)
}

get_distr_rdom_P_same_contig <- function(ite, bait_names, info_TargGene_fil, Gn=T, PMT=T) {
    PMTs_ind <- grep("^PMT_", bait_names)

    # 1) randomly chose the remaining genes before contig selection
    if (Gn){
    Gns <- bait_names[-PMTs_ind]
    rdom_Gns_rk <- get_rdom_Gn_rk(Gns, N_cate_rfGn, bait_names)
    noms_a <- names(rdom_Gns_rk[rdom_Gns_rk!=0])}

    # 2) randomly chose the remaining PMTs before contig selection
    if (PMT){
    PMTs <- bait_names[PMTs_ind]
    rdom_PMTs_rk <- get_rdom_PMT_rk(PMTs, N_cate_rfGn, bait_names)
    noms_b <- names(rdom_PMTs_rk[rdom_PMTs_rk!=0])}

    # 3) compute P
    if (Gn & PMT)
        {noms1 <- noms_a ; noms2 <- noms_b }
    else if (!PMT)
        {noms1 <- noms_a ; noms2 <- noms_a }
    else if (!Gn)
        {noms1 <- noms_b ; noms2 <- noms_b }

    if (ite%%100==0) print(paste(ite, "iterations done"))
    return(get_P_same_contig(noms1, noms2, info_TargGene_fil, verbose=F))
}

get_rdom_Gn_rk <- function(gns, N_cate_rfGn, bait_names){
    # 1) set up global and informative set of genes
    N_gns <- length(gns)    # the number of informative genes from the 'best baits per gene' data set (step 2.3)
    N_inf_gns <- sum(N_cate_rfGn) - as.numeric(N_cate_rfGn["PMT"])

    # 2) set up importance vector
    imp <- sample(1:N_gns, N_gns)
    names(imp) <- gns
    
    imp[imp<(N_gns+1-N_inf_gns)] <- 0   # should stay ~120 genes with some importance
    imp <- sort(imp, decreasing=T)  # has to be sorted for the function 'get_BestBaitperContig'
    return(imp)
}

get_rdom_PMT_rk <- function(bait_names, pmts, N_cate_rfGn, bait_names){
    # 1) set up global and informative set of PMTs
    N_pmts <- length(pmts)    # the number of informative PMTs from the 'best baits per gene' data set
    N_inf_pmts <- as.numeric(N_cate_rfGn["PMT"])

    # 2) set up importance vector
    imp <- sample(1:N_pmts, N_pmts)
    names(imp) <- pmts
    
    imp[imp<(N_pmts+1-N_inf_pmts)] <- 0   # should stay ~120 genes with some importance
    imp <- sort(imp, decreasing=T)  # has to be sorted for the function 'get_BestBaitperContig'
    return(imp)
}

get_baits_per_pairs <- function(bait_names, P_PMT, P_Gn, P_Gn_PMT, info_TargGene_fil, gini_gene_rf, ds)
#draw pairs of baits randomly but by respecting the probabilities of being on the same contig
{
    # 1) set up the number of PMTs and Gns
        # informative
    info_TargGene <- read.delim(info_TargGene_fil, stringsAsFactors=F)
    ind_inf_PMT <- grep("^PMT_", names(gini_gene_rf))
    inf_PMT <- names(gini_gene_rf)[ind_inf_PMT]
    inf_PMT_exon <- sapply(inf_PMT, collapse_elements, sep="_", what=1:3, colla="_")
    Ninf_PMT <- length(ind_inf_PMT)
    N_inf_Gns <- length(gini_gene_rf) - Ninf_PMT

        # total
    ind_all_PMT <- grep("^PMT_", bait_names)
    all_PMT <- bait_names[ind_all_PMT]
    all_PMT_exon <- sapply(all_PMT, collapse_elements, sep="_", what=1:3, colla="_")
    N_all_PMT <- length(all_PMT)
    all_Gns <- bait_names[-ind_all_PMT]
    N_all_Gns <- length(all_Gns)

    # 2) draw some initial random PMTs
    n_init <- round(1/5*Ninf_PMT)
    PMTs <- all_PMT[sample(1:N_all_PMT, n_init)]

    # 3) check % same contigs
    while (length(PMTs) != Ninf_PMT) {
        PMT_exons <- sapply(PMTs, collapse_elements, sep="_", what=1:3, colla="_")
        vcontigs_PMT <- sapply(PMT_exons, get_contig, info_TargGene)
        m_pairtest_PMT <- outer(vcontigs_PMT, vcontigs_PMT, "==")
        v_pairtest_PMT <- m_pairtest_PMT[upper.tri(m_pairtest_PMT, diag = FALSE)]
        count <- table(v_pairtest_PMT)
        test <- as.numeric(count[2]/sum(count))
        print(round(test, 4))

        if (test < (P_PMT- ds*0.5)) {
            print(paste("# test=", round(test,4), ", so draw baits from common contigs"))
            # pick 5 baits from the same contigs the baits already selected
            contig_test <- unique(subset(info_TargGene, NewTargetName%in%PMT_exons)[,"contigV2"])
            tab <- info_TargGene[info_TargGene$contigV2%in%contig_test,]
            tab <- subset(tab, NewTargetName%in%all_PMT_exon)
            tab <- subset(tab, !NewTargetName%in%PMT_exons)

            if (nrow(tab)!=0){
                new_exons <- tab$NewTargetName[sample(1:nrow(tab), 5)]
                new_baits <- all_PMT[sapply(new_exons, function(x) grep(paste("^", x, "_", sep=""), all_PMT))]
            } else {
                print("#### no other bait available on common contigs, so draw baits at random")
                veve <- all_PMT[!all_PMT%in%PMTs]
                new_baits <- veve[sample(1:length(veve), 5)]
            }
        } else {
            print("# draw baits at random")
            veve <- all_PMT[!all_PMT%in%PMTs]
            new_baits <- veve[sample(1:length(veve), 5)]
        }
        if (length(PMTs)> Ninf_PMT-5){
            ins <- Ninf_PMT-length(PMTs)
            PMTs <- c(PMTs, new_baits[1:ins])
        } else {
            PMTs <- c(PMTs, new_baits)}
    }
}


nullHDraw <- function(ite, g, n_info, bait_names, info_TargGene, verbose=F){

    PMTs_ind <- grep("^PMT_", bait_names)
    # 1) setup the number of remaining genes before contig selection
    Gns <- bait_names[-PMTs_ind]
    rdom_Gns_rk <- get_rdom_Gn_rk(Gns, N_cate_rfGn, bait_names)

    # 2) setup the number of remaining PMTs before contig selection
    PMTs <- bait_names[PMTs_ind]
    rdom_PMTs_rk <- get_rdom_PMT_rk(PMTs, N_cate_rfGn, bait_names)

    # 2) keep only a rankable importance for the best bait of non discarded genes only
        # attribute an EQUAL low importance of 0 for all other baits
        # corresponds to step 2.1) in the main script
    imp[imp<(n+1-n_info)] <- 0   # should stay ~240 genes with some importance
    imp <- sort(imp, decreasing=T)  # has to be sorted for the function 'get_BestBaitperContig'

    # 3) select best bait per contig
    imp1 <- imp[imp!=0]
    if (verbose) cat(imp1, sep="\n", file="importance.txt")
    all_targ <- sapply(names(imp1), collapse_elements, sep="_", what=1:3, colla="_")
    
#~    info_TargGene <- info_TargGene[infoTarg$NewTargetName%in%all_targ,]
    info_TargGene <- get_info_AllTargGene(info_TargGene, all_targ)
    
    if (verbose) cat(names(imp1), sep="\n", file="genes.txt")
    bestBait_PerContig <- as.character(get_1baitPerContig(info_TargGene, imp1, verbose=F))
    if (verbose) print(bestBait_PerContig)
    ind <- PrePro_findIndex(bestBait_PerContig, names(imp))
    
    # 4) keep only a rankable importance for the best bait best bait per contig
        # attribute an EQUAL low importance of 0 for all other baits
        # corresponds to step 3.1) in the main script
    imp[-ind] <- 0

    # 5) compute sum of ranks per grp
    df_fc <- data.frame(rnk = rank(imp), grp = g)
    res <- aggregate(rnk ~ grp,df_fc,sum)
    v <- res$rnk
    names(v) <- res$grp
    
    if (ite%%100==0) print(paste(ite, "iterations done"))
    return(v)
}

addZeroImpGenes <- function(xx){

	import_genes <- sapply(names(xx),function(x){
		st <- strsplit(x,"_")[[1]][1:2]
		return(paste(st[1],st[2],sep="_"))
		})

	all_genes <- getAllGeneNames()
	zero_imp <- all_genes[!(all_genes %in% import_genes)]
	v <- rep(0,length(zero_imp))
	names(v) <- zero_imp
	out <- c(xx,v)
	return(out)
}

getAllGeneNames <- function(){
    load(PREVIOUS_DATA2)
    data_mat <- alpha_matrix
	out <- unique(sapply(rownames(data_mat),function(x){
		st <- strsplit(x,"_")[[1]][1:2]
		return(paste(st[1],st[2],sep="_"))
		}))
	return(out)
}

get_pval <- function(observ, dist_exp, two_sided=F)
# if two_sided=F, P is always P(observ>expected)
{
    pval <- round(1-length(subset(dist_exp, observ>dist_exp))/length(dist_exp), 3)
    if (two_sided) {
        pval2 <- abs(1-pval)
        pval <- ifelse(pval<=pval2, pval, pval2)*2}
    return(pval)
}

get_1rdom_bait_per_gn <- function (gene_names){
    
    load(PREVIOUS_DATA2)
    noms <- rownames(alpha_matrix)
    res <- sapply(gene_names, grep_first, noms, value=T)
    return(res)
}

grep_first <- function(patt, vec, value=T){
    return(grep(paste("^", patt, sep=""), vec, value=value)[1])
}

############################################################# 

getGenePerCategTable <- function(){
	load(PREVIOUS_DATA2)
	N_gene_categ <- as.vector(gene_per_categ[,2])
	names(N_gene_categ) <- as.vector(as.character(gene_per_categ[,1]))
	return(N_gene_categ)
}


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
#	plot(tree, show.tip.label = T, tip.color=col_races, edge.color=edge_colors, edge.width=1.5, root.edge=F, cex = 0.45, font=4)

	if (lane) 
		legend("topright", legend=c(libraries, races_uniq), col = col_lib_raw, lty = c(rep(1, length(libraries)), rep(-1, length(race_colo_raw))), lwd=2, bg = 'gray92', cex=0.8, text.col=c(rep("black", length(libraries)), race_colo_raw))
	else 
		legend("topright", legend=races_uniq, col = race_colo_raw, lwd=2, bg = 'gray92', cex=0.8)
	title(paste("NJ tree of CNV distribution across aphid individuals\n", categ, " genes (", nbaits, " baits)", sep=""), cex.main=1.5, font=4)
}

processing_NJ <- function(i_categ, gene_categ, alpha_matrix, method, racine, lane, col_races, col_lib, libraries, races_uniq, col_lib_raw, race_colo_raw, baits)
{

	categ <- gene_categ[i_categ]
	talpha <- define_alpha_matrix(alpha_matrix, gene_categ, i_categ, baits)
	nbaits <- nrow(talpha)
	distance_matrix <-  compute_dstce(talpha, method)
	tree <- get_tree(distance_matrix, racine)
	
	tree_edges <- get_tree_edges(tree)
	edge_colors <- apply_colors_2_edges(tree, tree_edges, lane, col_races, col_lib)
	plot_tree(tree, col_races, edge_colors, lane, libraries, races_uniq, col_lib_raw, race_colo_raw, categ, nbaits)
	
	print(categ)
	print(dim(talpha))

}

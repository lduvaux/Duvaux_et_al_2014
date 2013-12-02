Count_Genes_targets <- function(tabl)
{

	nb_exons_categ <- table(tabl$GeneCateg)
	genes_names <- unique(tabl$Gene)
	mat_whole_gene <- Set_mat_whole_gene(nb_exons_categ)
	mat_targets <- Set_mat_targets(nb_exons_categ)
	for (ite in seq(length(genes_names)))
	{
		gen <- genes_names[ite]
		subtab <- subset(tabl, Gene==gen)
		ntarg <- nrow(subtab)

		categ <- subtab$GeneCateg[1]
		ind <- grep(paste("^", categ, "$", sep=""), rownames(mat_whole_gene))
		if (length(ind)>1) 
			stop("several matches")
		mat_targets[ind, 1] <- mat_targets[ind, 1]+ntarg
		
		# test
		vtest <- is.na(subtab$contigV2)
		ntarg_removed <- length(grep(T, vtest))
		if (ntarg_removed>ntarg)
			stop("ntarg_removed>ntarg")
		mat_targets[ind, 2] <- mat_targets[ind, 2]+ntarg_removed
		test1 <- T%in%vtest
		test2 <- ifelse(test1 & length(grep(T, vtest)) == nrow(subtab), T, F)
		
		if (test2)
			mat_whole_gene[ind, 3] <- mat_whole_gene[ind, 3]+1
		else if (!test1)
			mat_whole_gene[ind, 1] <- mat_whole_gene[ind, 1]+1
		else
			mat_whole_gene[ind, 2] <- mat_whole_gene[ind, 2]+1
	}
	res  <- list(mat_whole_gene=mat_whole_gene, mat_targets=mat_targets)
	return(res)
}


Set_mat_whole_gene <- function(Nb_exons_categ)
{
	tab <- matrix(length(Nb_exons_categ), 3, data=0)
	rownames(tab) <- names(Nb_exons_categ)
	colnames(tab) <- c("Complete", "Partially_deleted", "Removed")
	return(tab)
}

Set_mat_targets <- function(Nb_exons_categ)
{
	tab <- matrix(length(Nb_exons_categ), 2, data=0)
	rownames(tab) <- names(Nb_exons_categ)
	colnames(tab) <- c("Targets", "Targets_Removed")
	return(tab)
}




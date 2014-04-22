fetchRacesPrior <- function(tabrace){
    tabrace <- subset(tabrace, tabrace[,"Included.in.CNV.RF.training.set.and.GLMM.analyses"]==1)
    vrace <- sapply(tabrace$Sample.Name, get_elements)
    return(vrace)
}

Trim_Alpha_mat <- function(prior_assign, raw_alpha_mat)
# remove individuals designated as irrelavant for statistics uses
{
	# 1) fetch info about individuals
	IndivRace_raw <- fetchRacesPrior(prior_assign)
    
	# 2) keep only good individuals
    good <- pmatch(names(IndivRace_raw), colnames(raw_alpha_mat))
	alpha_matrix <- raw_alpha_mat[,good]
	return(alpha_matrix)
}

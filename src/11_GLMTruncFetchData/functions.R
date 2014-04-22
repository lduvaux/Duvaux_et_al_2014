Trim_Alpha_mat <- function(PRIOR_ASSIGN, alpha_matrix_raw)
# remove individuals designated as irrelavant for statistics uses
{
	# 1) fetch info about individuals
	IndivRace_raw <- PrePro_fetchRacesPrior(PRIOR_ASSIGN)
    
	# 2) check correspondance between alpha matrix and the info matrix
	ind_order <- PrePro_findIndex(colnames(alpha_matrix_raw), names(IndivRace_raw)) # Y
	IndivRace_raw <- IndivRace_raw[ind_order] # Y
	if (!identical(names(IndivRace_raw), colnames(alpha_matrix_raw)))
		stop("Problem of correspondance between alpha_matrix & PRIOR_ASSIGN")
        
	# 3) remove bad individuals
	bad_indiv_ind <- which(is.na(IndivRace_raw)|IndivRace_raw=="hybrid")
	alpha_matrix <- alpha_matrix_raw[,-bad_indiv_ind]
	return(alpha_matrix)
}

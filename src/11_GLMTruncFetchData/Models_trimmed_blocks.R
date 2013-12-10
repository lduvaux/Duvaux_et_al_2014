# models
	# with trimmed
MOD_ALL1 <- Duplication ~ Race + LnGeneLength + LnExonLength + Family + trimmed + LnGeneLength * Family + LnExonLength * Family + Family * trimmed + Race * Family

FIXED_TERMS <- c("Race","trimmed", "Family")
M_MAX <- 6

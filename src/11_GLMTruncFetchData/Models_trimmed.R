# models
	# with trimmed
MOD_ALL1 <- Duplication ~ LnGeneLength + LnExonLength + Family + trimmed + LnGeneLength * Family + LnExonLength * Family + Family * trimmed

FIXED_TERMS <- c("trimmed", "Family")
M_MAX <- 5

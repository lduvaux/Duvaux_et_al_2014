# models
	# with trimmed

MOD_ALL1 <- fqcy ~ LnGeneLength + LnExonLength + Family + trimmed + CompDup + LnGeneLength * Family + LnExonLength * Family + Family * trimmed + Family * CompDup + (1|race)

FIXED_TERMS <- c("trimmed", "Family", "CompDup")
M_MAX <- 6

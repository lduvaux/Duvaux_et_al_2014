# models
	# with trimmed

#~MOD_ALL1 <- Fqcy_all ~ LnGeneLength + LnExonLength + Family + trimmed + CpDup + LnGeneLength * Family + LnExonLength * Family + Family * trimmed + Family * CpDup + (1|race)  # I have removed all the interactions as they don't appear significant
MOD_ALL1 <- Fqcy_all ~ LnGeneLength + LnExonLength + Family + trimmed + CpDup + (1|race)
MOD_ALL2 <- Fqcy_all ~ LnGeneLength + LnExonLength + Family + trimmed + CpDup + (1|race)
FAMILY <- "binomial"
FIXED_TERMS <- c("trimmed", "Family", "CpDup")
M_MAX <- 6


#~fm1 <- lme4(fqcy ~ LnGeneLength + LnExonLength + Family + trimmed + CpDup + LnGeneLength * Family + LnExonLength * Family + Family * trimmed + Family * CpDup + (1|race), data = GLMtab_all2, family=binomial)
#~test_trimmed <- dredge(fm1, fixed=c("trimmed", "Family", "CpDup"), m.max=7)

invlogit(x) # binomial

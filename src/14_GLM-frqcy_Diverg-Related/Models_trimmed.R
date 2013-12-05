# models
	# with trimmed

MOD_ALL1 <- Fqcy_all ~ LnGeneLength + LnExonLength + Family + CpDup + LnGeneLength * Family + LnExonLength * Family + Family * CpDup + (1|race) + (1|trimmed)

# interactions:
# LnGeneLength:Family => small interaction
# LnExonLength:Family => no interaction
# FamilyP450:trimmedYes => 0.055
# FamilyXx:CpDupTRUE => no interaction

MOD_ALL2<- Fqcy_all ~ LnGeneLength + LnExonLength + Family + CpDup + LnGeneLength * Family + (1|race) + (1|trimmed)

FAMILY <- "binomial"
FIXED_TERMS <- c("trimmed", "Family", "CpDup")
M_MAX <- 6

#~invlogit(x) # binomial

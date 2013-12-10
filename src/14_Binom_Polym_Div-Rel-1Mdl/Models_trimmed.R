
# 1) model binomial 1 (duplication events)
MOD_ALL1 <- Duplication ~ LnGeneLength + LnExonLength + Family + Phylog_lvl + Family * Phylog_lvl + LnGeneLength * LnExonLength + (1|Phylog_lvl/Race) + (1|trimmed) + (1|Family/Gene) # + (1|LnGeneLength/LnExonLength)
FAMILY <- "binomial"
DELTA1 <- 10
FIXED_TERMS1 <- NULL  # dredge is quite long for MOD_ALL1 and 'Familly' term strongly significant for MOD_ALL1

M_MAX <- 10

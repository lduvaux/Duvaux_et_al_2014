
# 1) model binomial step 2 (complete duplication events)
MOD_ALL2 <- Duplication ~ LnGeneLength + LnExonLength + Family + Phylog_lvl + Family * Phylog_lvl + LnGeneLength * Family + LnExonLength * Family + (1|Phylog_lvl/Race) + (1|trimmed) + (1|Family/Gene) # + (1|LnGeneLength/LnExonLength) # virer Phylog_lvl car jamais significatif
FAMILY <- "binomial"
DELTA2 <- 10
FIXED_TERMS2 <- NULL
M_MAX <- 10


# 1) model binomial 1 (duplication events)
MOD_ALL1 <- Polymorphic ~ LnGeneLength + LnExonLength + Family + LnGeneLength * Family + LnExonLength * Family + (1|Race) + (1|trimmed) + (1|Family/Gene) # + (1|LnGeneLength/LnExonLength)
FAMILY <- "binomial"
DELTA1 <- 10
FIXED_TERMS1 <- c("Family")  # dredge is quite long for MOD_ALL1 and 'Familly' term strongly significant for MOD_ALL1

# 2) model binomial 2 (complete duplication events)
MOD_ALL2 <- Duplication ~ LnGeneLength + LnExonLength + Family + LnGeneLength * Family + LnExonLength * Family + (1|Race) + (1|trimmed) + (1|Family/Gene) # + (1|LnGeneLength/LnExonLength)
DELTA2 <- c(10, 10)
FIXED_TERMS2 <- NULL

M_MAX <- 10

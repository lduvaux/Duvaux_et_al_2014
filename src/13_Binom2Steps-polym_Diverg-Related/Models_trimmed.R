
# 1) model binomial 1 (duplication events)
MOD_ALL1 <- Dup ~ LnGeneLength + LnExonLength + Family + LnGeneLength * Family + LnExonLength * Family + (1|Race) + (1|trimmed)


# 2) model binomial 2 (complete duplication events)
MOD_ALL2<- Polymorphism ~ LnGeneLength + LnExonLength + Family + LnGeneLength * Family + LnExonLength * Family + (1|Race) + (1|trimmed)

FAMILY <- "binomial"

FIXED_TERMS <- c("Family")
M_MAX <- 5

#~invlogit(x) # binomial

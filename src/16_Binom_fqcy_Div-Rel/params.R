#~ All explit parameters for the related script.R should be declared here as constantes ~#

# 1) fetch data
PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")
CLUSTERS <- list(divergent=c("Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Ononis"), related=c("Medicago", "Pisum", "Trifolium"))
SUBCLUSTERS <- list(Cytisus="Cytisus", Lathyrus="Lathyrus", L.corn.="L.corn.", L.ped.="L.ped.", Ononis="Ononis", Medicago="Medicago", Pisum="Pisum", Trifolium="Trifolium")

# 2) models for pair plots
MODP2 <- Fqcy_all ~ LnIntronLength + LnExonLength + Family + trimmed + Race + Phylog_lvl + CpDup

# 3) model binomial step 2 (complete duplication?)
MOD_ALL2 <- Fqcy_all ~ LnExonLength + Family + CpDup + Phylog_lvl + Family * LnExonLength + Family * CpDup + (1|Race) + (1|Gene)  # trimmed removed as never significant ; LnIntronLength removed as never significant ; FamilyGr:Phylog_lvl removed as never significant 
FAMILY <- "binomial"
DELTA2 <- 10
FIXED_TERMS2 <- NULL
M_MAX <- 10

# 3) outputs
PAIRS_ALL_EXON_LENGTH <- "./Res16_PairPlot_DupFqcy.pdf"

GLM_DUPFQCY_MAX <- "./GLM16-1_FqcyMax.txt"
GLM_DUPFQCY_DREDGE <- "./GLM16-2_DupFqcyDredge.txt"
GLM_DUPFQCY_BEST <- "./GLM16-3_CpDupFqcybestMdl.txt"
GLM_DUPFQCY_AVG <- "./GLM16-4_DupFqcyAvg.txt"
DUPFQCY_PDF <- "./Res16-5_PredictPrCpDupFqcy.pdf"

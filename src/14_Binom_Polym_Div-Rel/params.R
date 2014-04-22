#~ All explit parameters for the related script.R should be declared here as constantes ~#

# 0) hardware params
N_CORE <- 8

# 1) fetch data
PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")
CLUSTERS <- list(divergent=c("Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Ononis"), related=c("Medicago", "Pisum", "Trifolium"))
SUBCLUSTERS <- list(Cytisus="Cytisus", Lathyrus="Lathyrus", L.corn.="L.corn.", L.ped.="L.ped.", Ononis="Ononis", Medicago="Medicago", Pisum="Pisum", Trifolium="Trifolium")

# 2) models for pair plots
MODP2 <- Polymorphic ~ LnIntronLength + LnExonLength + Family + trimmed + Race + Phylog_lvl

# 3) model binomial 1 (polymorphic or not polymorphic?)
MOD_ALL1 <- Polymorphic ~ LnIntronLength + LnExonLength + Family + Phylog_lvl + trimmed + (1|Race) + (1|Gene) # 'lengthX:family' interactions removed as 'lengthX' never significant, idem for + LnIntronLength * LnExonLength, + Family * Phylog_lvl 
FAMILY <- "binomial"
DELTA1 <- 10
FIXED_TERMS1 <- NULL
M_MAX <- 8

# 4) outputs
PAIRS_ALL_EXON_LENGTH <- "./Res14-0_PairPlot_Polym.pdf"

GLM_POL_MAX <- "./GLM14-1_PolMax.txt"
GLM_POL_DREDGE <- "./GLM14-2_PolDredge.txt"
GLM_POL_BEST <- "./GLM14-3_PolbestMdl.txt"
GLM_POL_AVG <- "./GLM14-4_PolAvg.txt"
POL_PDF <- "./Res14-4_PredictPrPolymorph.pdf"

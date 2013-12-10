#~ All explit parameters for the related script.R should be declared here as constantes ~#

# 1) fetch data
PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")
CLUSTERS <- list(divergent=c("Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Ononis"), related=c("Medicago", "Pisum", "Trifolium"))
SUBCLUSTERS <- list(Cytisus="Cytisus", Lathyrus="Lathyrus", L.corn.="L.corn.", L.ped.="L.ped.", Ononis="Ononis", Medicago="Medicago", Pisum="Pisum", Trifolium="Trifolium")


# 2) models for pair plots
MODP2 <- Duplication ~ LnGeneLength + LnExonLength + Family + trimmed + Polymorphic + Race + Phylog_lvl

# 3) outputs
PAIRS_ALL_EXON_LENGTH <- "./Res13-0_PairPlot_Dup.pdf"

GLM_POL_MAX <- "./GLM13-1_PolMax.txt"
GLM_POL_DREDGE <- "./GLM13-2_PolDredge.txt"
GLM_POL_AVG <- "./GLM13-3_PolAvg.txt"
POL_PDF <- "./Res13-4_PredictPrPolymorph.pdf"

GLM_CPDUP_MAX <- "./GLM_13-5_CpDupMax.txt"
GLM_CPDUP_DREDGE <- "./GLM_13-6_CpDupDredge.txt"
GLM_CPDUP_AVG <- "./GLM_13-7_CpDupAvg.txt"
CPDUP_PDF <- "./Res13-8_PredictPrCpDup.pdf"

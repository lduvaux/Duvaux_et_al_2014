#~ All explit parameters for the related script.R should be declared here as constantes ~#

# 1) fetch data
PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")
CLUSTERS <- list(divergent=c("Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Ononis"), related=c("Medicago", "Pisum", "Trifolium"))
SUBCLUSTERS <- list(Cytisus="Cytisus", Lathyrus="Lathyrus", L.corn.="L.corn.", L.ped.="L.ped.", Ononis="Ononis", Medicago="Medicago", Pisum="Pisum", Trifolium="Trifolium")


# 2) models for pair plots
MODP2 <- Duplication ~ LnGeneLength + LnExonLength + Family + trimmed + Race + Phylog_lvl

# 3) outputs
PAIRS_ALL_EXON_LENGTH <- "./Res15_PairPlot_DupType.pdf"

GLM_DUP_MAX <- "./GLM15-1_DupMax.txt"
GLM_DUP_DREDGE <- "./GLM15-2_DupDredge.txt"
GLM_DUP_AVG <- "./GLM15-3_DupAvg.txt"
DUP_PDF <- "./Res15-4_PredictPrCpDup.pdf"

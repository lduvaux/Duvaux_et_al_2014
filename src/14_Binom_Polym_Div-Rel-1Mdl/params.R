#~ All explit parameters for the related script.R should be declared here as constantes ~#

# 1) fetch data
PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")
CLUSTERS <- list(divergent=c("Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Ononis"), related=c("Medicago", "Pisum", "Trifolium"))
SUBCLUSTERS <- list(Cytisus="Cytisus", Lathyrus="Lathyrus", L.corn.="L.corn.", L.ped.="L.ped.", Ononis="Ononis", Medicago="Medicago", Pisum="Pisum", Trifolium="Trifolium")


# 2) models for pair plots
MODP2 <- Polymorphic ~ LnGeneLength + LnExonLength + Family + trimmed + Race + Phylog_lvl

# 3) outputs
PAIRS_ALL_EXON_LENGTH <- "./Res_Pair_Polym_ExLgCovar_AllGenes_DRall.pdf"

GLM_POL_MAX <- "./GLM_PolMax_DRall.txt"
GLM_POL_DREDGE <- "./GLM_PolDredge_DRall.txt"
GLM_POL_AVG <- "./GLM_PolAvg_DRall.txt"
POL_PDF <- "./Res_PredictProbaOfPollication_DRall.pdf"

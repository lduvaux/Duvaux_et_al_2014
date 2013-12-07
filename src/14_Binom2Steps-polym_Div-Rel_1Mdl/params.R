#~ All explit parameters for the related script.R should be declared here as constantes ~#

# 1) fetch data
PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")
CLUSTERS <- list(divergent=c("Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Ononis"), related=c("Medicago", "Pisum", "Trifolium"))
SUBCLUSTERS <- list(Cytisus="Cytisus", Lathyrus="Lathyrus", L.corn.="L.corn.", L.ped.="L.ped.", Ononis="Ononis", Medicago="Medicago", Pisum="Pisum", Trifolium="Trifolium")


# 2) models for pair plots
MODP2 <- Polymorphism ~ LnGeneLength + LnExonLength + Family + trimmed + Dup + Race + Phylog_lvl

# 3) outputs
PAIRS_ALL_EXON_LENGTH <- "./Pair_Polym_ExLgCovar_AllGenes_DRall.pdf"

GLM_DUP_MAX <- "./GLM_DupMax_DRall.txt"
GLM_DUP_DREDGE <- "./GLM_DupDredge_DRall.txt"
GLM_DUP_AVG <- "./GLM_DupAvg_DRall.txt"
DUP_PDF <- "./PredictProbaOfDuplication_DRall.pdf"

GLM_CPDUP_MAX <- "./GLM_CpDupMax_DRall.txt"
GLM_CPDUP_DREDGE <- "./GLM_CpDupDredge_DRall.txt"
GLM_CPDUP_AVG <- "./GLM_CpDupAvg_DRall.txt"
CPDUP_PDF <- "./PredictProbaOfCpDuplication_DRall.pdf"

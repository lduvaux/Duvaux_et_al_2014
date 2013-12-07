#~ All explit parameters for the related script.R should be declared here as constantes ~#

# 1) fetch data
PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")
CLUSTERS <- list(divergent=c("Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Ononis"), related=c("Medicago", "Pisum", "Trifolium"))
SUBCLUSTERS <- list(Cytisus="Cytisus", Lathyrus="Lathyrus", L.corn.="L.corn.", L.ped.="L.ped.", Ononis="Ononis", Medicago="Medicago", Pisum="Pisum", Trifolium="Trifolium")


# 2) models for pair plots
MODP2 <- Polymorphism ~ LnGeneLength + LnExonLength + Family + trimmed + Dup + Race

# 3) outputs
PAIRS_ALL_EXON_LENGTH <- "./Pair_Polym_ExLgCovar_AllGenes.pdf"

GLM_DUP_MAX <- "./GLM_DupMax.txt"
GLM_DUP_DREDGE <- "./GLM_DupDredge.txt"
GLM_DUP_AVG <- "./GLM_DupAvg.txt"
DUP_PDF <- "./PredictProbaOfDuplication.pdf"

GLM_CPDUP_MAX <- "./GLM_CpDupMax.txt"
GLM_CPDUP_DREDGE <- "./GLM_CpDupDredge.txt"
GLM_CPDUP_AVG <- "./GLM_CpDupAvg.txt"
CPDUP_PDF <- "./PredictProbaOfCpDuplication.pdf"

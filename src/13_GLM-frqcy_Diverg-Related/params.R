#~ All explit parameters for the related script.R should be declared here as constantes ~#

PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")

CLUSTERS <- list(divergent=c("Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Ononis"), related=c("Medicago", "Pisum", "Trifolium"))
SUBCLUSTERS <- list(Cytisus="Cytisus", Lathyrus="Lathyrus", L.corn.="L.corn.", L.ped.="L.ped.", Ononis="Ononis", Medicago="Medicago", Pisum="Pisum", Trifolium="Trifolium")

# models
	# for pair plots
MODP2 <- Fqcy_all ~ LnGeneLength + LnExonLength + Family + trimmed + CpDup + race

	# without trimmed
# outputs
PAIRS_ALL_RATIO <- "./Pair_Polym_RatioCovar_AllGenes.pdf"
PAIRS_ALL_EXON_LENGTH <- "./Pair_Polym_ExLgCovar_AllGenes.pdf"
INTERACTION_HIST <- "./Hist_FactInteract_AllGenes.pdf"

GLM_RES_BEST_TRIMMED <- "./GLM_CNVpolym_All_Best_trimmed.txt"
GLM_RES_TRIMMED_NoINTERACT <- "./GLM_CNVpolym_All_trimmed_NoInteract.txt"

HEAT_MAP_COMP <- "ProbaCompleteDup_HeatMap_3levs.pdf"
HEAT_MAP_PART <- "ProbaPartialDup_HeatMap_3levs.pdf"

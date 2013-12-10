#~ All explit parameters for the related script.R should be declared here as constantes ~#

# 1) fetch data
PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")
CLUSTERS <- list(divergent=c("Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Ononis"), related=c("Medicago", "Pisum", "Trifolium"))


# 2) models for pair plots
MODP2 <- Duplication ~ LnGeneLength + LnExonLength + Family + trimmed

# 3) outputs
PAIRS_ALL_EXON_LENGTH <- "./Res_Pair_Polym_ExLgCovar_AllGenes.pdf"
INTERACTION_HIST <- "./Res_Hist_FactInteract_AllGenes.pdf"
BOXPLOTS_PDF <- "./Res_PolymDistribution.pdf"
GLM_RES_BEST_TRIMMED <- "./Res_GLM_CNVpolym_All_Best_trimmed.txt"

#~ All explit parameters for the related script.R should be declared here as constantes ~#

PREVIOUS_DATA <- "../11_GLMTruncFetchData.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
CATEG_FOR_GLM <- c("Control", "Gr", "Or","P450")

# outputs
PAIRS_COMPLETE_EXON_LENGTH <- "./Pair_Polym_ExLgCovar_CompleteGenesOnly.pdf"
INTERACTION_HIST <- "./Hist_FactInteract_AllGenes.pdf"
INTERACTION_HIST_COMPLETE <- "./Hist_FactInterac_CompleteGenesOnly.pdf"

GLM_RES_BEST_TRIMMED <- "./GLM_CNVpolym_All_Best_trimmed_2lev.txt"
GLM_RES_AVG_TRIMMED <- "./GLM_CNVpolym_All_Avg_trimmed_2lev.txt"
GLM_RES_TRIMMED_NoINTERACT <- "./GLM_CNVpolym_All_trimmed_NoInteract_2lev.txt"

HEAT_MAP_COMP <- "ProbaCompleteDup_HeatMap_2levs.pdf"
HEAT_MAP_PART <- "ProbaPartialDup_HeatMap_2levs.pdf"

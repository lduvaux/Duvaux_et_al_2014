#~ All explicit parameters for the related script.R should be declared here as constantes ~#

source("../05_Stability-RF/params.R")
PREVIOUS_DATA <- "../03_ProcessingSegmentation.Rdata"
#~RAW_DATA <- "../01_GetData.Rdata"

# outputs
TAB_CNV_DISTR <- "Res_CNVs_distrib_subtargets.txt"
P_CNV1 <- "Res_Proba_CNV_per_marker.csv"
BARPLOT1.1 <- "Res_Fig2a_Distrib_CN_GbalSample.jpg"
BARPLOT1.2 <- "Res_Fig2a_Distrib_CN_GbalSample.pdf"

VENN1.1 <- "Res_Fig2b_Venn_CNV_per_marker.jpg"
VENN1.2 <- "Res_Fig2b_Venn_CNV_per_marker.pdf"

#~ All explicit parameters for the related script.R should be declared here as constantes ~#

source("../05_Stability-RF/params.R")
PREVIOUS_DATA <- "../03_ProcessingSegmentation.Rdata"
#~RAW_DATA <- "../01_GetData.Rdata"

# outputs
TAB_CNV_DISTR <- "CNVs_distrib_subtargets.txt"
P_CNV1 <- "Proba_CNV_per_marker.txt"
VENN1.1 <- "Fig2_Venn_CNV_per_marker.jpg"
VENN1.2 <- "Fig2_Venn_CNV_per_marker.pdf"

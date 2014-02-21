#~ All explit parameters for the related script.R should be declared here as constantes ~#

PREVIOUS_DATA <- "../03_ProcessingSegmentation.Rdata"
RAW_DATA <- "../01_GetData.Rdata"

Rdata <- "result_prosegment.Rdata"
fil <- "Stability_RF.Rdata"
N_RF <- 100
NTREES <- 100

SAMP_SIZE1 <- c(5, rep(10, 5), 9, 5)
SAMP_SIZE2 <- c(5, rep(10, 5), 9, 5)

PDF_NAME1 <- "Res_10_bestVariab_superv.pdf"
PDF_NAME2 <- "Res_10_bestVariab_unsuperv.pdf"

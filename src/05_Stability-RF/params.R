#~ All explit parameters for the related script.R should be declared here as constantes ~#

PREVIOUS_DATA <- "../03_ProcessingSegmentation.Rdata"
RAW_DATA <- "../01_GetData.Rdata"

Rdata <- "result_prosegment.Rdata"
fil <- "Stability_RF.Rdata"
N_RF <- 100
NTREES <- 100

# sample size in RF
    # races are sorted as follow:
    # "Cytisus", "Lathyrus", "L.corn.", "L.ped.", "Medicago", "Ononis", "Pisum", "Trifolium"
SAMP_SIZE1 <- c(5, rep(10, 5), 8, 6)
SAMP_SIZE2 <- c(5, rep(10, 5), 8, 6)

PDF_NAME1 <- "Res_10_bestVariab_superv.pdf"
PDF_NAME2 <- "Res_10_bestVariab_unsuperv.pdf"

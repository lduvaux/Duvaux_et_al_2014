#~ All explit parameters for the related script.R should be declared here as constantes ~#

PREVIOUS_DATA <- "../03_ProcessingSegmentation.Rdata"
RAW_DATA <- "../01_GetData.Rdata"

Rdata <- "result_prosegment.Rdata"
fil <- "Stability_RF.Rdata"
N_RF <- 100
NTREES <- 100
RACE_UNKNOWN <- c("Medicago_204_T60","Pisum_248_T98","Pisum_193_T94","Lathyrus_167_T18","L.corn._243_T36","Trifolium_N319_T118","L.ped._N184_T54","Ononis_277_T85","Pisum_246_T96","Medicago_221_T69","Medicago_138_T57","Lathyrus_282_T25","Pisum_247_T97","Pisum_192_T93","Medicago_N143_T74","L.corn._239_T35","Trifolium_136_T110","Pisum_249_T99","Medicago_N341_T78", "Ononis_N123_T89")

SAMP_SIZE1 <- c(5, rep(12, 7))
SAMP_SIZE2 <- c(5, rep(10, 5), 9, 10)

PDF_NAME1 <- "Res_10_bestVariab_superv.pdf"
PDF_NAME2 <- "Res_10_bestVariab_unsuperv.pdf"

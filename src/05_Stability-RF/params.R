#~ All explit parameters for the related script.R should be declared here as constantes ~#

PREVIOUS_DATA <- "../03_ProcessingSegmentation.Rdata"
RAW_DATA <- "../01_GetData.Rdata"

Rdata <- "result_prosegment.Rdata"
fil <- "Stability_RF.Rdata"
N_RF <- 100
NTREES <- 100
RACE_UNKNOWN <- c("L.corn._239_T35", "L.corn._243_T36", "L.ped._N184_T54", "Lathyrus_167_T18", "Lathyrus_282_T25", "Ononis_277_T85", "Pisum_192_T93", "Pisum_193_T94", "Pisum_246_T96", "Pisum_247_T97", "Pisum_248_T98", "Pisum_249_T99", "Trifolium_106_T106", "Trifolium_136_T110", "Trifolium_198_T111", "Trifolium_47_T112", "Trifolium_48_T113", "Trifolium_57_T114", "Trifolium_N102_T116", "Trifolium_N185_T117", "Trifolium_N322_T119", "Trifolium_N324_T120")  # Med2 as been tested as unknown => very similar to Med1, so not included here anymore.

SAMP_SIZE1 <- c(5, rep(10, 5), 9, 5)
SAMP_SIZE2 <- c(5, rep(10, 5), 9, 5)

PDF_NAME1 <- "Res_10_bestVariab_superv.pdf"
PDF_NAME2 <- "Res_10_bestVariab_unsuperv.pdf"

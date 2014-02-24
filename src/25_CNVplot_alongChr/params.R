#~ All explit parameters for the related script.R should be declared here as constantes ~#
SUBTARG <- "../00_RawData/CaptureTargets.PropSeq.bait51.interval_list.noheader"
TARG <- "../00_RawData/MappingOnV2.NewExonName_20131009.txt"
PREVIOUS_DATA <- "../03_ProcessingSegmentation.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
#~GN_FIG <- c("Control_g130", "Control_g194", "Control_g12", "Control_g49", "Or_g54", "Or_g21", "Control_g86", "Gr_g42", "Gr_g21", "P450_g24", "P450_g47", "P450_g11", "P450_g13", "Control_g163", "PMT_313", "Control_g56")
GN_FIG2 <- c("Control_g194", "Control_g49", "Gr_g21", "Control_g56", "Or_g6", "P450_g13")
YLI <- c(1.26, 2.5, 4.7, 1.6, 1.5, 1.5)   # top limit for paper plot
TRANSP <- F
ALPHA <- 120

# outputs
PDF1 <- "Plot_CNV_along_Chr.pdf"
PDF2 <- "Fig3_CNV_along_Chr_paper.pdf"
FIL_LGTH <- "LengthDifference_V1-V2.csv"

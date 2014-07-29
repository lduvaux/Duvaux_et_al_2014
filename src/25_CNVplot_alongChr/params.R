#~ All explit parameters for the related script.R should be declared here as constantes ~#
SUBTARG <- "../00_RawData/02_baitsFiles/CaptureTargets.PropSeq.bait51.interval_list.noheader"
TARG <- "../00_RawData/MappingOnV2.NewExonName_20131009.txt"
PREVIOUS_DATA <- "../03_ProcessingSegmentation.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
#~GN_FIG <- c("Control_g130", "Control_g194", "Control_g12", "Control_g49", "Or_g54", "Or_g21", "Control_g86", "Gr_g42", "Gr_g21", "P450_g24", "P450_g47", "P450_g11", "P450_g13", "Control_g163", "PMT_313", "Control_g56")
GN_FIG2 <- c("Control_g194", "Control_g49", "Gr_g21", "P450_g13")
MAT_LAYOUT <- matrix(1:4, nrow=2, ncol=2, byrow=T)
YLI1 <- -0.25   # bottom limit for paper plot
YLI2 <- c(1.26, 2.5, 4.7, 1.6, 1.5, 1.5)   # top limit for paper plot
TRANSP <- F
ALPHA <- 120
SKIP_MAIN <- F

# outputs
PDF1 <- "Res_Plot_CNV_along_Chr.pdf"
PDF2 <- "Fig1_CNV_along_Chr_paper.pdf"
JPG <- "Fig1_CNV_along_Chr_paper.jpg"
FIL_LGTH <- "Res_LengthDifference_V1-V2.csv"

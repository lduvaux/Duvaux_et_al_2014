source("../utils/getter_functions.R")
DEBUG <- F
VERBOSE <- 1

CTRL_GUYS <- c("Medicago_204_T60", "Medicago_32_T72", "Medicago_207_T63", "Medicago_209_T64", "Medicago_212_T67", "Medicago_N334_T76", "Medicago_29_T71", "Medicago_36_T73", "Medicago_205_T61", "Medicago_214_T68")# the third first from lane "light green", the 3 following from the "purple" lane and the last 4 are from the first lane (try to smooth the lane effect amongst several lanes). They are all from Medicago2.

BAD_GUYS <- c("L.ped._82_T47", "L.ped._84_T49", "Medicago_210_T65", "Pisum_5_T100")

BAD_CYTISUS <- c("Cytisus_115_T1", "Cytisus_127_T2", "Cytisus_128_T3", "Cytisus_16_T7", "Cytisus_76_T10", "Cytisus_77_T11", "Cytisus_79_T12", "Cytisus_80_T13", "Cytisus_87_T14", "Cytisus_89_T15", "Cytisus_90_T16", "Cytisus_92_T17")

RACE_UNKNOWN <- c("L.corn._239_T35", "L.corn._243_T36", "L.ped._N184_T54", "Lathyrus_167_T18", "Lathyrus_282_T25", "Ononis_277_T85", "Pisum_192_T93", "Pisum_193_T94", "Pisum_246_T96", "Pisum_247_T97", "Pisum_248_T98", "Pisum_249_T99", "Trifolium_106_T106", "Trifolium_136_T110", "Trifolium_198_T111", "Trifolium_47_T112", "Trifolium_48_T113", "Trifolium_57_T114", "Trifolium_N185_T117", "Trifolium_N322_T119", "Trifolium_N324_T120")  # Med2 as been tested as unknown => very similar to Med1, so not included here anymore.

PRIOR_ASSIGN <- read.delim("../00_RawData/TableS3_SamplingAndAssignement.csv", stringsAsFactors =F)
INFO_TARGENE_FILE <- "../00_RawData/MappingOnV2.NewExonName_20131009.txt"

INDIV_DETAILS <- read.delim("../00_RawData/01_Sequencing_stats/02-1_SummarySequencingPerindividual.txt", stringsAsFactors =F)
BAIT_GOOD_Q10RATIO <- "../00_RawData/02_baitsFiles/02-1_Target.GoodRatio.baits51.Q10.txt"	# 



race_colo_raw <- c("black", "brown", "cyan", "orange", "blue", "darkmagenta", "coral", "green", "red")
race_colo_raw_noctrl <- c("black", "brown", "cyan", "orange", "darkmagenta", "coral", "green", "red")
col_lib_raw <- c("red", "blue", "plum1", "cyan", "coral", "green", "darkmagenta", "gray40","mediumseagreen", "yellow4")


ver <- function(str){
try(
	if(VERBOSE > 0)
		print(str)
	)
}

dbgPrint <- function(str){
    if(DEBUG)
        print(str)
}

source("../utils/getter_functions.R")
DEBUG <- F
VERBOSE <- 1

CTRL_GUYS <- c("Medicago_204_T60", "Medicago_32_T72", "Medicago_207_T63", "Medicago_209_T64", "Medicago_212_T67", "Medicago_N334_T76", "Medicago_29_T71", "Medicago_36_T73", "Medicago_205_T61", "Medicago_214_T68")# the third first from lane "light green", the 3 following from the "purple" lane and the last 4 are from the first lane (try to smooth the lane effect amongst several lanes). They are all from Medicago2.

BAD_GUYS <- c("L.ped._82_T47", "L.ped._84_T49", "Medicago_210_T65")


BAD_CYTISUS <- c("Cytisus_115_T1", "Cytisus_127_T2", "Cytisus_128_T3", "Cytisus_16_T7", "Cytisus_76_T10", "Cytisus_77_T11", "Cytisus_79_T12", "Cytisus_80_T13", "Cytisus_87_T14", "Cytisus_89_T15", "Cytisus_90_T16", "Cytisus_92_T17")


INDIV_DETAILS <- read.delim("../00_RawData/02-1_SummarySequencingPerindividual.txt", stringsAsFactors =F)
PRIOR_ASSIGN <- read.delim("../00_RawData/00_filtering_good_Baits/IndivAssign_STRUCTUREandRF_20130814.final.csv", stringsAsFactors =F)
BAIT_GOOD_Q10RATIO <- "../00_RawData/00_filtering_good_Baits/02-1_Target.GoodRatio.baits51.Q10.txt"	# 
RAW_DATA_FILE <- "../00_RawData/02-0.Data4PropoSeq.Rdata"
INFO_TARGENE_FILE <- "../00_RawData/MappingOnV2.NewExonName_20131009.txt"


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

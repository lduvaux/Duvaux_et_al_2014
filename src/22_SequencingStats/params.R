#~ All explit parameters for the related script.R should be declared here as constantes ~#
SEQ_DEPTH <- "../00_RawData/01_Sequencing_stats/CoveragePerTargetPerInd.txt"
FASTQC_STATS <- "../00_RawData/01_Sequencing_stats/02_SummarySequencing_FastQC_pair-end.txt"
BAD_FILE <- "111214_0267_D067GACXX_1_SA-PE-020"
BAD_FCELL <- "D0CM0ABXX"
BAD_INDIV <- c(82, 84, 210, 5)  # correspond to the Genepool directory (second number in the clone name)
    # insert sizes
INS_PATH <- "../00_RawData/01_Sequencing_stats/01_bamstats"
PATTERN <- "InsSizeMetrics$"

# outputs
OUTFIL1 <- "./Res_MedianSqcingDepth.csv"
OUTFIL2 <- "./TableS4_SequencinqStatistics.csv"

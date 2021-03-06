#~ All explit parameters for the related script.R should be declared here as constantes ~#
# params
GENE_CATEG <- c("Or", "Gr", "OBP", "IR", "CSP", "ApisSNMP", "P450", "SP", "PS", "Control", "Promoter")
SETS <- c("Targets", "Genes", "GMAP_Remaining_targets", "GMAP_Remaining_genes", "GMAP_NonTrimmed_genes", "Final_Remaining_targets", "Final_NonTrimmed_targets", "Final_Remaining_genes", "Final_NonTrimmed_genes")

# inputs
ALL_TARG_INFO <- "../00_RawData/MappingOnV2.NewExonName_20131009.txt"
TARG_UNIQ_NOVERLAP <- "../00_RawData/Targets_uniq_V2_NoOverlap.txt"
PRO_SEGMENT <- "../03_ProcessingSegmentation.Rdata"
BAITS_BE4_CLEAN <- "../00_RawData/02_baitsFiles/CaptureTargets.PropSeq.bait51.interval_list.noheader"

# outputs
OUTFIL <- "./TableS2_SummaryTargetLoci.csv"

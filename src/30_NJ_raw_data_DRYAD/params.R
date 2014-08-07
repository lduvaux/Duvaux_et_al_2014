#~ All explit parameters for the related script.R should be declared here as constantes ~#
POOL_INFO <- "../00_RawData/experiment_paired_fastq_spreadsheet_Duvaux-et-al-2014.tsv"
PREVIOUS_DATA <- "../03_ProcessingSegmentation.Rdata"
RAW_DATA <- "../01_GetData.Rdata"
METHOD <- "manhattan"	# "manhattan" or "ludo"
ROOT <- "Lathyrus_281_T24"
PLOT_LANES <- T
PLOT_POOLS <- T
LGD_POS <- c("topright", "topright", "bottomleft", "topright",
             "topright", "bottomleft", "topright", "bottomright",
             "bottomright", "topright", "bottomright", "bottomright")
PDF_NAME <- "Res_CNV-NJ_PoolsAndLanes.pdf"

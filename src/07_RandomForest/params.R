#~ All explit parameters for the related script.R should be declared here as constantes ~#

PREVIOUS_DATA <- "../04_Stability-RF.Rdata"
PREVIOUS_DATA2 <- "../03_ProcessingSegmentation.Rdata"
GENE_GER_CATEG_DATA <- "../07_GenesPerCateg.Rdata"
RAW_DATA <- "../01_GetData.Rdata"

NTREES <- 5000
TWOSIDED <- F

# outputs
TAB50EXON <- "./Res_Tab_best50_exons.csv"
TABIMPEXON <- "./Res_RF_Imptce_exons.csv"
TAB20EXON <- "./Res_Tab_best20_genes.csv"
TABIMPGn <- "./Res_RF_imptce_genes.csv"
TAB20CONTIG <- "./Res_Tab_best20_contigs.csv"
TABIMPCONTIG <- "./Res_RF_imptce_contigs.csv"
ASSIGN_TAB <- "./Res_Indiv_Assign_contig_rf.csv"
PLOTDIAGGn <- "./Res_lastFig_GeneDistinguishRaces.pdf"

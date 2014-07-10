#!/bin/Rscript
library('parallel')
    # load global ressources
rm(list=ls())    
source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

    # load local ressources
source("./params.R")
source("./functions.R")
load(RF_DATA)
load(ALP_DATA)
#~load(CDD_DATA)
source("./params.R")
source("./functions.R")


main <- function(argv){
	# load data
    tab_targ <- read.delim(TARG, stringsAsFactors=F)
    InfSubtarg <- set_infoSubtarg(tem_tab)

    # 1) find gene name and family
        # find family
    Family <- sapply(InfSubtarg$Subtarget, get_elements)
        # find real gene name
    TargetName <- sapply(InfSubtarg$Subtarget, collapse_elements)
    v_ind <- match(TargetName, tab_targ$NewTargetName)
    Gene <- tab_targ$Gene[v_ind]

    # 2) new columns
        # 2.1) MDR
    tab_inf <- InfSubtarg[,2:9]
    MDR_ind <- apply(tab_inf, 1, which.max)
    MDR <- colnames(tab_inf)[MDR_ind]
    RF_rank <- 1:length(MDR)

        # 2.2) CDD MDR
    Gene_nom <- sapply(InfSubtarg$Subtarget, collapse_elements, what=1:2)
    alpha_matrix_RF_full <- set_alpha_mat(alpha_matrix, new_x3)
    good <- which(Family!="PMT")
    CDD_MDR <- rep(NA, length(Gene_nom))
    CDD_MDR[good] <- sapply(good, function(x) get_CDD_MDR(race=MDR[x], gene=Gene_nom[x], alp_mat=alpha_matrix_RF_full))

        # 2.3) MCN
    MCN_tab <- sapply(races_uniq, get_MCN, new_x3)
    MCN_tab <- MCN_tab[match(InfSubtarg$Subtarget, rownames(MCN_tab)),]

        # 2.4) CN MDR
    CN_MDR <- sapply(1:nrow(MCN_tab), function(x) MCN_tab[x, match(MDR[x],colnames(MCN_tab))])

    # 3) final table
    temp_tab <- as.data.frame(matrix(ncol=ncol(tab_inf)+ncol(MCN_tab), nrow=nrow(InfSubtarg), data=NA))
    temp_tab[,seq(1, ncol(temp_tab), by=2)] <- MCN_tab
    colnames(temp_tab)[seq(1, ncol(temp_tab), by=2)] <- paste("MCN", colnames(MCN_tab))
    temp_tab[,seq(2, ncol(temp_tab), by=2)] <- tab_inf
    colnames(temp_tab)[seq(2, ncol(temp_tab), by=2)] <- paste("RF score", colnames(tab_inf))
    f_InfSubtarg <- cbind(Gene, Family, Subtarget=InfSubtarg$Subtarget, MDR, CDD_MDR, CN_MDR, RF_rank, MDG=InfSubtarg$MeanDecreaseGini, temp_tab, MDA=InfSubtarg$MeanDecreaseAccuracy)
    write.table(f_InfSubtarg, file=OUTPUT1, quote=F, row.names=F, sep="\t")


    # 4) count misc stats relative to discriminating loci
        # 4.1) detect putative pseudogenes
    gr40 <- length(grep("Gr", head(Family, 40), fixed=T))
    or40 <- length(grep("Or", head(Family, 40), fixed=T))
    pseudo_gr40 <- length(grep("^Gr[0-9]{1,2}.*[P,C,I]$|^Gr[0-9]{1,2}[P,C,I].*$", head(Gene, 40), value=T))
    print(paste(pseudo_gr40, "out of the", gr40, "Gr genes in the top 40 are pseudogenes"))
    pseudo_or40 <- length(grep("^Or[0-9]{1,2}.*[P,C,I]$|^Or[0-9]{1,2}[P,C,I].*$", head(Gene, 40), value=T))
    print(paste(pseudo_or40, "out of the", or40, "Or genes in the top 40 are pseudogenes"))
   
        # 4.2) number of deletions for the MDR in the top40
    del40 <- length(which(head(CN_MDR, 40)<1))
    print(paste("Among the variants characterizing the MDR in the top 40,", del40, "are deletions"))

        # 4.3) proportion of the duplication in the top 40
    dup <- length(which(CN_MDR>1))
    dup40 <- length(which(head(CN_MDR, 40)>1))
    print(paste(dup40, "out of the", dup, "cases where MDR are characterized by duplications are in the top 40"))

        # 4.4.) distribution of CDD genes
    n_CDD <- length(which(CDD_MDR))
    n_CDD40 <- length(which(head(CDD_MDR, 40)))
    print(paste(n_CDD40, "out of the", n_CDD, "cases where MDR are characterized by CDD variants are in the top 40"))

    cat("\n")
    print(" #### save results")
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    dummy <- numeric()
    save(f_InfSubtarg, alpha_matrix_RF_full, file=outFileName)

}

cat("logfile\n")

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










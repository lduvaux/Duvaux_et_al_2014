#!/bin/Rscript
rm(list=ls())
source("../utils/functions.R")
source("../utils/globalCtes.R")
source("params.R")
source("./functions.R")

main <- function(argv){

    cat("\n")
    print(" #### 1) load data and colors")
        # 1.0) load subtargets coordinates
    subtarg <- read.table(SUBTARG, stringsAsFactors=F)
    refBaits <- sapply(subtarg$V5, Prepro_fixName)
    subtarg[,"V5"] <- refBaits
    colnames(subtarg) <- c("Contig", "Start", "End", "Strand", "Name")
    t_targ <-  read.delim(TARG, stringsAsFactors=F)
    
        # 1.1) load reads count
    load(PREVIOUS_DATA)
        # 1.2) remove bad cytisus
    ind_bad <- PrePro_findIndex(BAD_CYTISUS, colnames(alpha_matrix))
    ind_bad2 <- PrePro_findIndex(c(BAD_CYTISUS, RACE_UNKNOWN), colnames(alpha_matrix))
    alpha_matrix <- alpha_matrix[,-ind_bad]
    print(dim(alpha_matrix))

        # 1.3) set up marker names
    subtargets <- rownames(alpha_matrix)
    targets <- unique(sapply(subtargets, collapse_elements))
    genes0 <- unique(sapply(targets, collapse_elements, what=1:2))
    bad <- grep("PMT", genes0, fixed=T)
    genes <- genes0[-bad]
#~    PMT <- genes0[bad]

    # 2) test with the first gene
    l_gn <- c(genes[1], "Control_g190", "Control_g204", "Gr_g36", "Control_g144", "Control_g181")
    pdf("Plot_CNV_along_Chr.pdf")
    layout(matrix(1:6, nrow=3, ncol=2, byrow=T))
    par(mar=c(4, 4, 2, 2), mgp=c(1.5,0.5,0))
    for (gn in l_gn)
    {
        ind <- grep (paste(gn, "_", sep=""), t_targ[,"NewTargetName"])
        lt_targ <- t_targ[ind, ]
        
        ind <- grep (paste(gn, "_", sep=""), subtarg[,"Name"])
        lt_star <- subtarg[ind, ]
        
        ind <- grep (paste(gn, "_", sep=""), rownames(alpha_matrix))
        lcnv <- alpha_matrix[ind,]
        n_ctig <- sort(table(lt_star$Contig), decreasing=T)

        # check contig
        if (length(n_ctig)>1) {
        print("WARNING: gene on several contigs")
        gd_ctig <- names(n_ctig)[1]
        ind <- lt_star$Contig%in%gd_ctig
        lt_star <- lt_star[ind,]
        gd_star <- lt_star[ind,"Name"]
        gd <- sapply(gd_star, function(x) which(x==rownames(lcnv)))
        lcnv <- lcnv[gd, ]

        ind <- lt_targ$contigV2%in%gd_ctig
        lt_targ <- lt_targ[ind,]
        }

        plot_CNV_chr(tab_cnv=lcnv, tab_star=lt_star, tab_tar=lt_targ, c_ex=0.5, l_wd=.5)
    }
    dev.off()

}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);

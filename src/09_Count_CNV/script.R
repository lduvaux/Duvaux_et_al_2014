#!/bin/Rscript
rm(list=ls())
source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("../utils/randomForest_helperFuns.R")
source("params.R")
source("./functions.R")
library(VennDiagram)
library(grid)
library(gridBase)
library(lattice)

main <- function(argv){
	load(PREVIOUS_DATA)
	set.seed(0)

    cat("\n")
    print(" #### 1) load data and colors")
        # 1.1) reads count
	alpha_matrix <- PrePro_roundToZeroFive(alpha_matrix)
        # 1.2) remove bad cytisus
    ind_bad <- PrePro_findIndex(BAD_CYTISUS, colnames(alpha_matrix))
    ind_bad2 <- PrePro_findIndex(c(BAD_CYTISUS, RACE_UNKNOWN), colnames(alpha_matrix))
    alpha_matrix <- alpha_matrix[,-ind_bad]
    print(dim(alpha_matrix))

    clones <- colnames(alpha_matrix)
    races <- unique(sapply(clones, function(x) unlist(strsplit(x, "_"))[1]))
    subtargets <- rownames(alpha_matrix)
    targets <- unique(sapply(subtargets, collapse_elements))
    genes0 <- unique(sapply(targets, collapse_elements, what=1:2))
    bad <- grep("PMT", genes0, fixed=T)
    genes <- genes0[-bad]
    PMT <- genes0[bad]
    categ <- unique(sapply(genes0, function(x) unlist(strsplit(x, "_"))[1]))

    cat("\n")
    print(" #### 2) count global proportion of deletion/insertions/both")
        print(" # 2.0) Compute CN distribution")
    CNV_distr0 <- table(as.vector(alpha_matrix))
    CNV_distr <- CNV_distr0[-3]
    CNV_distr_pro <- round(CNV_distr/sum(CNV_distr), 5)*100
    tab <- rbind(names(CNV_distr_pro), CNV_distr_pro)
    write.table(tab, file=TAB_CNV_DISTR, sep="\t", row.names=F, quote=F, col.names=F)
    
    
        print(" # 2.1) subtargets")
    CNV_count_bait <- t(apply(alpha_matrix, 1, caract_bait))
    colnames(CNV_count_bait) <- c("Deletion", "Insertion", "Both")
    toto <- apply(CNV_count_bait, 2, sum)
    N_bait <- nrow(CNV_count_bait)
    n_bait <- sum(toto[1:2])-toto[3]
    l4venn_bait <- list(Deletion = 1:(toto[1]), Insertion = (toto[1]-toto[3]+1):n_bait)
    p_bait <- n_bait/N_bait
        
    cat("\n")
    print(" # 2.2) targets")
    CNV_count_targ <- t(sapply(targets, get_caract_targ, subtargets, alpha_matrix))
    colnames(CNV_count_targ) <- c("Deletion", "Insertion", "Both")
    toto_targ <- apply(CNV_count_targ, 2, sum)
    n_targ <- (sum(toto_targ[1:2])-toto_targ[3])
    N_targ <- nrow(CNV_count_targ)
    l4venn_targ <- list(Deletion = 1:(toto_targ[1]), Insertion = (toto_targ[1]-toto_targ[3]+1):n_targ)
    p_targ <- n_targ/N_targ
    
    cat("\n")
    print(" # 2.3) genes")
    CNV_count_gn <- t(sapply(genes, get_caract_targ, subtargets, alpha_matrix))
    colnames(CNV_count_gn) <- c("Deletion", "Insertion", "Both")
    toto_gn <- apply(CNV_count_gn, 2, sum)
    n_gn <- (sum(toto_gn[1:2])-toto_gn[3])
    N_gn <- nrow(CNV_count_gn)
    l4venn_gn <- list(Deletion = 1:(toto_gn[1]), Insertion = (toto_gn[1]-toto_gn[3]+1):n_gn)
    p_gn <- n_gn/N_gn

    cat("\n")
    print(" # 2.4) promoters")
    CNV_count_pmt <- t(sapply(PMT, get_caract_targ, subtargets, alpha_matrix))
    colnames(CNV_count_pmt) <- c("Deletion", "Insertion", "Both")
    toto_pmt <- apply(CNV_count_pmt, 2, sum)
    n_pmt <- (sum(toto_pmt[1:2])-toto_pmt[3])
    N_pmt <- nrow(CNV_count_pmt)
    l4venn_pmt <- list(Deletion = 1:(toto_pmt[1]), Insertion = (toto_pmt[1]-toto_pmt[3]+1):n_pmt)

    p_pmt <- n_pmt/N_pmt


    cat("\n")
    print(" # 2.5) per gene family")
    CNV_count_gn_pmt <- rbind(CNV_count_gn, CNV_count_pmt)
    test_pol_gn <- sapply(1:nrow(CNV_count_gn_pmt), function(x) 1%in%CNV_count_gn_pmt[x,])
    tab_pol_gn <- CNV_count_gn_pmt[test_pol_gn,]
    mat_pol <- matrix(data=NA, nrow=length(categ), ncol=3, dimnames=list(categ, c("Nb of loci", "Nb of CNV", "Proportion")))
    for (i in seq(length(categ)))
    {
        ca <- categ[i]
        # total
        ind <- grep(paste(ca, "_", sep=""), rownames(CNV_count_gn_pmt))
        N_cat <- nrow(CNV_count_gn_pmt[ind,])

        # polym
        ind <- grep(paste(ca, "_", sep=""), rownames(tab_pol_gn))
        n_cat <- nrow(tab_pol_gn[ind,])
        if (is.null(n_cat)) n_cat <- 0

        mat_pol[i,] <- c(N_cat, n_cat, round(n_cat/N_cat,3))
    }
    mat_pol <- cbind(rownames(mat_pol), mat_pol)

    m_prob <- cbind(c("Subtargets", "Targets", "Genes", "Promoters"), c(N_bait, N_targ, N_gn, N_pmt), c(n_bait, n_targ, n_gn, n_pmt), round(c(p_bait, p_targ, p_gn, p_pmt),3))
    colnames(m_prob) <- c("Marker", "Nb of loci", "Nb of CNV", "Proportion")
    m_prob <- rbind(m_prob, mat_pol)
    write.table(m_prob, file=P_CNV1, sep="\t", quote=F, row.names=F)


    cat("\n")
    print(" # 2.6) draw venns")
    l_venns <- list(Subtargets=l4venn_bait, Targets=l4venn_targ, Genes=l4venn_gn, Promoters=l4venn_pmt)

    jpeg(VENN1.1)
    draw_venn(l_venns)
    dev.off()
    pdf(VENN1.2)
    draw_venn(l_venns)
    dev.off()


#~    cat("\n")
#~    print(" #### 3) count global proportion of deletion/insertions/both per race")
#~        # 3.1) subtargets
#~    CNV_count_race <- mclapply(1:nrow(alpha_matrix), function(x) sapply(races, caract_bait_race, clones, alpha_matrix[x,], x))
#~    names(CNV_count_race) <- subtargets

    #############################
    cat("\n")
    print(" #### save results")
	outFileName <- argv[1]
    ver(sprintf("Saving *DUMMY* data to %s",outFileName))
    dummy <- numeric()
    save(p_bait, p_targ, p_pmt, p_gn,file=outFileName)
}



argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);










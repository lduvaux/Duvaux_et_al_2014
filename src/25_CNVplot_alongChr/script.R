#!/bin/Rscript
rm(list=ls())
set.seed(0)
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
    v_races <- PrePro_fetchRaces(RAW_DATA, CTRL_GUYS, BAD_GUYS)
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
    PMT <- genes0[bad]

    # 2) Draw all genes
#~    l_gn <- c(genes[1], "Control_g190", "Control_g204", "Gr_g36", "Control_g144", "Control_g181")
    l_gn <- genes
    print(system.time({
        pdf(PDF1)
        layout(matrix(1:6, nrow=3, ncol=2, byrow=T))
        par(mar=c(4, 4, 2, 2), mgp=c(1.5,0.5,0))
        tab <- as.data.frame(matrix(data=NA, ncol=4, nrow=length(genes), dimnames=list(genes, c("Gene", "LengthV1", "LengthV2", "Fold"))))
        for (i in seq(l_gn))
        {
            # prepare data
            gen <- l_gn[i]
            l_data <- get_data4plot(gen, t_targ, subtarg, alpha_matrix, adapt=T)
            tab[i,] <- l_data$stats

            # draw plot
#~            pdf("test.pdf")
#~              tab_cnv <- l_data$lcnv; tab_star <- l_data$lt_star; tab_tar <- l_data$lt_targ; yli=l_data$rgg; centz=c(0.75,1.25); c_ex=.9; l_wd=.9; races=v_races
            plot_CNV_chr(tab_cnv=l_data$lcnv, tab_star=l_data$lt_star, tab_tar=l_data$lt_targ, c_ex=0.5, l_wd=.5, yli=l_data$rgg, races=v_races)
            cat ("###################")
            if (i%%20==0) {cat("\n") ; print(i)}
            cat ("###################\n")
#~            dev.off()
        }
        dev.off()
    }))
    ind <- order(as.numeric(tab[,"Fold"]))
    tabf <- tab[ind,]
    write.table(tabf, file=FIL_LGTH, sep="\t", row.names=F, quote=F)


    # 3) plots for the paper
    l_gn <- GN_FIG2
    m_rg <- cbind(rep(-0.5, 6), YLI)
    print(system.time({
        pdf(PDF2)
        layout(matrix(1:6, nrow=3, ncol=2, byrow=T))
        par(mar=c(4, 4, 2, 2), mgp=c(1.5,0.5,0))
        for (i in seq(l_gn))
        {
            # prepare data
            gen <- l_gn[i]
            l_data <- get_data4plot(gen, t_targ, subtarg, alpha_matrix, adapt=T)

            # draw plot
#~            pdf("test.pdf")
#~              tab_cnv <- l_data$lcnv; tab_star <- l_data$lt_star; tab_tar <- l_data$lt_targ; yli=l_data$rgg; centz=c(0.75,1.25); c_ex=.9; l_wd=.9; races=v_races
            colos_race <- plot_CNV_chr(tab_cnv=l_data$lcnv, tab_star=l_data$lt_star, tab_tar=l_data$lt_targ, c_ex=0.5, l_wd=.5, yli=m_rg[i,], races=v_races, transp=TRANSP, alphaa=ALPHA)
            cat ("###################")
            if (i%%20==0) {cat("\n") ; print(i)}
            cat ("###################\n")
#~            dev.off()
        }
        dev.off()
    }))


    #############################
    cat("\n")
    print(" #### save results")
	outFileName <- argv[1]
    ver(sprintf("Saving *DUMMY* data to %s",outFileName))
    dummy <- numeric()
    save (dummy, file=outFileName)
}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);

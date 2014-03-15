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
    t_targ0 <-  read.delim(TARG, stringsAsFactors=F)
    
        # 1.1) load reads count
    load(PREVIOUS_DATA)
    set.seed(0)
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

    # reduce t_targ
    v_new_gn_targ <- sapply(t_targ0$NewTargetName, collapse_elements, what=1:2)
    ind <- v_new_gn_targ%in%genes
    t_targ <- t_targ0[ind,]

    # reorder genes vector
        # reorder target table by V1 coordinates of contig and lower bp coordinate
    first_coord <- sapply(1:nrow(t_targ), function(x) min(as.numeric(t_targ[x,c("startV1", "stopV1")])))
    ord <- order(t_targ$contigV1, first_coord)
    t_targ2 <- t_targ[ord,]
        # reorder genes
    v_new_gn_targ2 <- sapply(t_targ2$NewTargetName, collapse_elements, what=1:2)
    ind_genes0 <- sapply(seq(genes), function(x) min(which(genes[x]==v_new_gn_targ2)))
    ind_genes <- order(ind_genes0)
    genes2 <- genes[ind_genes]
    

    cat("\n")
    print(" #### 2) Draw all genes")
    if (SKIP_MAIN==F){
        l_gn <- genes2
        print(system.time({
            pdf(PDF1)
            layout(matrix(1:6, nrow=3, ncol=2, byrow=T))
            par(mar=c(4, 4, 2, 2), mgp=c(1.5,0.5,0))
            tab <- as.data.frame(matrix(data=NA, ncol=4, nrow=length(l_gn), dimnames=list(l_gn, c("Gene", "LengthV1", "LengthV2", "Fold"))))
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
    }


    cat("\n")
    print(" #### 3) plots for the paper")
    l_gn <- GN_FIG2
    m_rg <- cbind(rep(YLI1, 6), YLI2)
    print(system.time({
        pdf(PDF2)
        layout(MAT_LAYOUT)
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

    jpeg(JPG, height=480*2, width=480*2, quality=100, res=72*2)
    {
        layout(MAT_LAYOUT)
        par(mar=c(4, 4, 2, 2), mgp=c(1.5,0.5,0))
        for (i in seq(l_gn))
        {
            # prepare data
            gen <- l_gn[i]
            l_data <- get_data4plot(gen, t_targ, subtarg, alpha_matrix, adapt=T)
            # draw plot
            colos_race <- plot_CNV_chr(tab_cnv=l_data$lcnv, tab_star=l_data$lt_star, tab_tar=l_data$lt_targ, c_ex=0.5, l_wd=.5, yli=m_rg[i,], races=v_races, transp=TRANSP, alphaa=ALPHA)
        }
        dev.off()
    }


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

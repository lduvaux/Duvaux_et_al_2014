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
#~    l_gn <- c(genes[1], "Control_g190", "Control_g204", "Gr_g36", "Control_g144", "Control_g181")
    l_gn <- genes
    print(system.time({
        pdf(PDF1)
        layout(matrix(1:6, nrow=3, ncol=2, byrow=T))
        par(mar=c(4, 4, 2, 2), mgp=c(1.5,0.5,0))
        tab <- as.data.frame(matrix(data=NA, ncol=4, nrow=length(genes), dimnames=list(genes, c("Gene", "LengthV1", "LengthV2", "Fold"))))
        for (i in seq(l_gn))
        {
            gn <- l_gn[i]
            ind <- grep (paste(gn, "_", sep=""), t_targ[,"NewTargetName"])
            lt_targ <- t_targ[ind, ]

            ind <- grep (paste(gn, "_", sep=""), subtarg[,"Name"])
            lt_star <- subtarg[ind, ]
            n_ctig <- sort(table(lt_star$Contig), decreasing=T)
            
            ind <- grep (paste(gn, "_", sep=""), rownames(alpha_matrix))
            lcnv <- alpha_matrix[ind,]
            if (is.vector(lcnv)) {
                lcnv <- t(as.matrix(lcnv))
                rownames(lcnv) <- rownames(alpha_matrix)[ind]
            }
            
            # check contig
            if (length(n_ctig)>1) {
                print("WARNING: gene on several contigs")
                
                # detect good subtarget in initial subtargets
                gd_ctig <- names(n_ctig)[1]
                ind <- lt_star$Contig%in%gd_ctig
                lt_star <- lt_star[ind,]    # subset initial table of subtargets
                
                # remove subtargets from bad contigs in lcnv
                    # fetch contigs of star in lt_cnv
#~                ind <- match(rownames(lcnv), subtarg[, "Name"])
                ind <- sapply(rownames(lcnv), function(x) which(x==subtarg[, "Name"]))
                cnv_contigs <- subtarg[ind,"Contig"]
                ind <- cnv_contigs%in%gd_ctig
                lcnv <- lcnv[ind, ]

                # same with lt_targ table
                ind <- lt_targ$contigV2%in%gd_ctig
                lt_targ <- lt_targ[ind,]
            }
            LengthV1 <- max(c(lt_targ$startV1,lt_targ$stopV1), na.rm=T) - min(c(lt_targ$startV1,lt_targ$stopV1), na.rm=T) + 1
            LengthV2 <- max(c(lt_targ$startV2,lt_targ$stopV2), na.rm=T) - min(c(lt_targ$startV2,lt_targ$stopV2), na.rm=T) + 1
            Fold <- round(LengthV1/LengthV2, 2)
            tab[i,] <- c(gn, LengthV1, LengthV2, Fold)

    #~        tab_cnv <- lcnv; tab_star <- lt_star; tab_tar <- lt_targ; yli=c(-0.5, 2.5); centz=c(0.75,1.25); c_ex=.9; l_wd=.9
            mmax <- max(lcnv)
            if (mmax< 2.5)
                rgg <- c(-0.5, 2.5)
            else
                rgg <- c(-0.5, ceiling(mmax))

#~            pdf("test.pdf")
            plot_CNV_chr(tab_cnv=lcnv, tab_star=lt_star, tab_tar=lt_targ, c_ex=0.5, l_wd=.5, yli=rgg)
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

# for the paper
    l_gn <- GN_FIG
    print(system.time({
        pdf(PDF2)
        layout(matrix(1:6, nrow=3, ncol=2, byrow=T))
        par(mar=c(4, 4, 2, 2), mgp=c(1.5,0.5,0))
        for (i in seq(l_gn))
        {
            gn <- l_gn[i]
            ind <- grep (paste(gn, "_", sep=""), t_targ[,"NewTargetName"])
            lt_targ <- t_targ[ind, ]
            
            ind <- grep (paste(gn, "_", sep=""), subtarg[,"Name"])
            lt_star <- subtarg[ind, ]
            n_ctig <- sort(table(lt_star$Contig), decreasing=T)
            
            ind <- grep (paste(gn, "_", sep=""), rownames(alpha_matrix))
            lcnv <- alpha_matrix[ind,]
            if (is.vector(lcnv)) {
                lcnv <- t(as.matrix(lcnv))
                rownames(lcnv) <- rownames(alpha_matrix)[ind]
            }
            
            # check contig
            if (length(n_ctig)>1) {
                print("WARNING: gene on several contigs")
                
                # detect good subtarget in initial subtargets
                gd_ctig <- names(n_ctig)[1]
                ind <- lt_star$Contig%in%gd_ctig
                lt_star <- lt_star[ind,]    # subset initial table of subtargets
                
                # remove subtargets from bad contigs in lcnv
                    # fetch contigs of star in lt_cnv
#~                ind <- match(rownames(lcnv), subtarg[, "Name"])
                ind <- sapply(rownames(lcnv), function(x) which(x==subtarg[, "Name"]))
                cnv_contigs <- subtarg[ind,"Contig"]
                ind <- cnv_contigs%in%gd_ctig
                lcnv <- lcnv[ind, ]

                # same with lt_targ table
                ind <- lt_targ$contigV2%in%gd_ctig
                lt_targ <- lt_targ[ind,]
            }

    #~        tab_cnv <- lcnv; tab_star <- lt_star; tab_tar <- lt_targ; yli=c(-0.5, 2.5); centz=c(0.75,1.25); c_ex=.9; l_wd=.9
            mmax <- max(lcnv)
            if (mmax< 2.5)
                rgg <- c(-0.5, 2.5)
            else
                rgg <- c(-0.5, ceiling(mmax))

#~            pdf("test.pdf")
            plot_CNV_chr(tab_cnv=lcnv, tab_star=lt_star, tab_tar=lt_targ, c_ex=0.5, l_wd=.5, yli=rgg)
            cat ("###################")
            if (i%%20==0) {cat("\n") ; print(i)}
            cat ("###################\n")
#~            dev.off()
        }
        dev.off()
    }))

    

}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);

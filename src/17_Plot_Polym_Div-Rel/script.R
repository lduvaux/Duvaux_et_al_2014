#!/bin/Rscript
rm(list=ls())
library(ggplot2)
library(lme4)
 
    # load global ressources
source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

    # load local ressources
source("./params.R")
source("./functions.R")

main <- function(argv){
	# load data
	load(PREVIOUS_DATA)
	load(PREVIOUS_DATA2)
    
    cat("\n")
    print("#### 1) proba per trimmed per gene family")
#~    print(samp_size_pol <- table(with(GLMtab_all1, interaction(Family, Polymorphic))))
    samp_size1 <- table(with(GLMtab_all1, interaction(Family, trimmed)))
    
    # 0) compute proba over-representation of multigene family genes
    genes <- unique(GLMtab_all1$Gene)
    Family <- as.factor(sapply(as.character(genes), get_elements, what=1))
    truncated <- !genes%in%NonTrimGenes
    tab <- data.frame(genes, truncated, Family)

        # 0.1)proportion of trimmed per category
    size_inter <- table(with(tab, interaction(Family, truncated)))
    size_fam <- table(tab$Family)
    p_trim <- size_inter/rep(size_fam,2)

        # 0.2) run GLM
    form <- truncated ~ Family
    glm2 <- glm(formula=form, data = tab, family="binomial")
    print(summary(glm2))


    # 1) prepare data for plot
    samp_size1_0 <- table(with(GLMtab_all1, interaction(Family, trimmed, Polymorphic)))
    obs_prob1.0 <- samp_size1_0/rep(samp_size1, 2)
    obs_prob1 <- obs_prob1.0[(length(obs_prob1.0)/2+1):length(obs_prob1.0)]

    colo <- c("Black", "Blue", "Purple", "DarkGreen")
    categ <- c("Control", "Gr", "Or", "P450")
    categ1 <- paste(categ, " (", samp_size1_0[(length(samp_size1_0)/2+1):length(samp_size1_0)], ")", sep="")

    jpeg(JPG, height=480*2, width=480*2, quality=100, res=72*2)
    par(mar=c(7, 6,4,2))
    draw_plot(obs_prob1, colo, categ1, yli=c(0, 0.55), yla="Proportion of observations\npolymorphic for CNV per family")
    text(x=c(2.5, 9), y=c(0.2, 0.5), labels=c("Non\ntruncated", "Truncated"), lheight=1.5)
    dev.off()
    pdf(PDF)
    par(mar=c(7, 6,4,2))
    draw_plot(obs_prob1, colo, categ1, yli=c(0, 0.55), yla="Proportion of observations\npolymorphic for CNV per family")
    text(x=c(2.5, 9), y=c(0.2, 0.5), labels=c("Non\ntruncated", "Truncated"), lheight=1.5)
    dev.off()


    cat("\n")
    print("#### 2) proba per trimmed per gene genetic proximity")
    samp_size2 <- table(with(GLMtab_all1, interaction(Phylog_lvl, trimmed)))
    samp_size2_0 <- table(with(GLMtab_all1, interaction(Phylog_lvl, trimmed, Polymorphic)))
    obs_prob2.0 <- samp_size2_0/rep(samp_size2, 2)
    obs_prob2 <- obs_prob2.0[(length(obs_prob2.0)/2+1):length(obs_prob2.0)]
    categ2 <- c("Divergent", "Related")
    categ2_0 <- paste(categ2, " (", samp_size2_0[(length(samp_size2_0)/2+1):length(samp_size2_0)], ")", sep="")

    jpeg(JPG2, height=480*2, width=480*2, quality=100, res=72*2)
    par(mar=c(7, 6,4,2))
    draw_plot(obs_prob2, colo=c("black", "white"), categ2_0, yli=c(0, 0.35), yla="Proportion of observations\npolymorphic for CNV per genetic proximity")
    text(x=c(1.5, 5), y=c(0.2, 0.32), labels=c("Non\ntruncated", "Truncated"), lheight=1.5)
    dev.off()

    pdf(PDF2)
    par(mar=c(7, 6,4,2))
    draw_plot(obs_prob2, colo=c("black", "white"), categ2_0, yli=c(0, 0.35), yla="Proportion of observations\npolymorphic for CNV per genetic proximity")
    text(x=c(1.5, 5), y=c(0.2, 0.32), labels=c("Non\ntruncated", "Truncated"), lheight=1.5)
    dev.off()


    cat("\n")
    print("#### 3) proba per trimmed per family per gene genetic proximity")
    samp_size3 <- table(with(GLMtab_all1, interaction(Family, Phylog_lvl, trimmed)))
    samp_size3_0 <- table(with(GLMtab_all1, interaction(Family, Phylog_lvl, trimmed, Polymorphic)))
    obs_prob3.0 <- samp_size3_0/rep(samp_size3, 2)
    obs_prob3 <- obs_prob3.0[(length(obs_prob3.0)/2+1):length(obs_prob3.0)]
    categ3 <- paste(categ, " (", samp_size3_0[(length(samp_size3_0)/2+1):length(samp_size3_0)], ")", sep="")

    jpeg(JPG3, height=480*2, width=480*2, quality=100, res=72*2)   
    par(mar=c(7, 6,4,2), lheight=1.2)
    draw_plot(obs_prob3, colo=colo, categ3, yli=c(0, 0.6), yla="Proportion of observations polymorphic\nfor CNV per family per genetic proximity", sep_ind=c(5, 9, 13))
    text(x=c(3, 9, 15, 21, 6, 18), y=c(.2, .2, .5, .5, .28, .55), labels=c(rep(c("Divergent", "Related"),2),"Non\ntruncated", "Truncated"))
    abline(v=11.5, lty=2)
    dev.off()

    pdf(PDF3)
    par(mar=c(7, 6,4,2), lheight=1.2)
    draw_plot(obs_prob3, colo=colo, categ3, yli=c(0, 0.6), yla="Proportion of observations polymorphic\nfor CNV per family per genetic proximity", sep_ind=c(5, 9, 13))
    text(x=c(3, 9, 15, 21, 6, 18), y=c(.2, .2, .5, .5, .28, .55), labels=c(rep(c("Divergent", "Related"),2),"Non\ntruncated", "Truncated"))
    abline(v=11.5, lty=2)
    dev.off()

    #########
    cat("\n")
    print(" #### save results")
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    dummy <- numeric()
    save(dummy,file=outFileName)

	
}

cat("dummy logfile\n")


argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);




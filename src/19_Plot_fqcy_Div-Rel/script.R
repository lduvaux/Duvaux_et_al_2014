#!/bin/Rscript
rm(list=ls())
library(ggplot2)

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
    set.seed(0)
    cat("\n")
    print("#### 1) Frequency of unusual CN in populations per gene families x trimmed categories")
    samp_size1 <- table(with(GLMtab_all2, interaction(Family, trimmed)))
    inter <- with(GLMtab_all2, interaction(Family, CpDup))
    samp_size <- table(inter)

    colo <- c("White", "Blue", "Purple", "DarkGreen")
    categ <- c("Control", "Gr", "Or", "P450")
    lab <- paste(rep(categ, 2), " (", samp_size, ")", sep="")

    jpeg(JPG, height=480*2, width=480*2, quality=100, res=72*2)
    p1 <- draw_plot(inter, GLMtab_all2$Fqcy_all, colo, categ1, yli=c(0, 1), yla="Frequency of unusual CN variants\nin populations", lab=lab, main="Family x CpDup", new_mar=c(9, 6,4,2), att=c(1:4, 6:9), adjust=.05)
    par(lheight=1.5)
    text(x=c(2.5, 7.5), y=c(0, 0), labels=c("Non CDD", "CDD"))
    abline(v=5, lty=2)
    dev.off()

    pdf(PDF)
    p1 <- draw_plot(inter, GLMtab_all2$Fqcy_all, colo, categ1, yli=c(0, 1), yla="Frequency of unusual CN variants\nin populations", lab=lab, main="Family x CpDup", new_mar=c(9, 6,4,2), att=c(1:4, 6:9), adjust=.05)
    par(lheight=1.5)
    text(x=c(2.5, 7.5), y=c(0, 0), labels=c("Non CDD", "CDD"))
    abline(v=5, lty=2)
    dev.off()



    cat("\n")
    print("#### 2) Frequency of unusual CN in populations per phylogenetic levels")
    jpeg(JPG2, height=480*2, width=480*2, quality=100, res=72*2)
    p1 <- draw_plot(GLMtab_all2$Phylog_lvl, GLMtab_all2$Fqcy_all, col=c("white", "grey"), categ1, yli=c(0, 1), yla="Frequency of unusual CN variants\nin populations", lab=c("Divergent", "Related"), main="Genetic relatedness", att=1:2)
    dev.off()
    
    pdf(PDF2)
    p1 <- draw_plot(GLMtab_all2$Phylog_lvl, GLMtab_all2$Fqcy_all, col=c("white", "grey"), categ1, yli=c(0, 1), yla="Frequency of unusual CN variants\nin populations", lab=c("Divergent", "Related"), main="Genetic relatedness", att=1:2)
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




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


    coord_txt_x <- c(2.5, 7.5)
    coord_txt_y <- c(0, 0)
    x_line <- 4.5
    jpeg(JPG, height=480*2, width=480*2, quality=100, res=72*2)
    p1 <- draw_plot(inter, GLMtab_all2$Fqcy_all, colo, categ1, yli=c(0, 1), lab=lab, main="", xla=XLAB1, yla=YLAB1, new_mar=c(9, 6,4,2), att=c(1:4, 6:9), adjust=.05, xline=x_line)
    par(lheight=1.5)
    text(x=coord_txt_x, y=coord_txt_y, labels=c("Non CDD", "CDD"))
    abline(v=5, lty=2)
    dev.off()

    pdf(PDF)
    p1 <- draw_plot(inter, GLMtab_all2$Fqcy_all, colo, categ1, yli=c(0, 1), lab=lab, main="", xla=XLAB1, yla=YLAB1, new_mar=c(9, 6,4,2), att=c(1:4, 6:9), adjust=.05, xline=x_line)
    par(lheight=1.5)
    text(x=coord_txt_x, y=coord_txt_y, labels=c("Non CDD", "CDD"))
    abline(v=5, lty=2)
    dev.off()



    cat("\n")
    print("#### 2) Frequency of unusual CN in populations per phylogenetic levels")
    x_line <- 4.2
    jpeg(JPG2, height=480*2, width=480*2, quality=100, res=72*2)
    p1 <- draw_plot(GLMtab_all2$Phylog_lvl, GLMtab_all2$Fqcy_all, col=c("white", "grey"), categ1, yli=c(0, 1), xla=XLAB2, yla=YLAB2, lab=c("Divergent", "Related"), main="", att=1:2, xline=x_line)
    dev.off()
    
    pdf(PDF2)
    p1 <- draw_plot(GLMtab_all2$Phylog_lvl, GLMtab_all2$Fqcy_all, col=c("white", "grey"), categ1, yli=c(0, 1), xla=XLAB2, yla=YLAB2, lab=c("Divergent", "Related"), main="", att=1:2, xline=x_line)
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




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

    # 1) generate a dataset from the estimated model
    samp_size <- table(with(GLMtab_all1, interaction(Family, trimmed)))
    obs_prob0 <- table(with(GLMtab_all1, interaction(Family, trimmed, Polymorphic)))/rep(samp_size, 2)
    obs_prob <- obs_prob0[(length(obs_prob0)/2+1):length(obs_prob0)]

    colo <- c("Black", "Blue", "Purple", "DarkGreen")
    categ <- rep(c("Control", "Gr", "Or", "P450"), 2)

    jpeg(JPG, height=480*2, width=480*2, quality=100, res=72*2)
    draw_plot(obs_prob, colo, categ, yli=c(0, 0.6))
    dev.off()
    pdf(PDF)
    draw_plot(obs_prob, colo, categ, yli=c(0, 0.6))
    dev.off()

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






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
    
    cat("\n")
    print("#### 1) proba of complete duplication/deletion per truncated/non-truncated gene")
    samp_size1 <- table(with(GLMtab_all2, interaction(Family, trimmed)))
    inter <- with(GLMtab_all2, interaction(Family, CpDup))
    samp_size <- table(inter)
    
    colo <- c("White", "Blue", "Purple", "DarkGreen")
    categ <- c("Control", "Gr", "Or", "P450")

    draw_plot(GLMtab_all2$Family, GLMtab_all2$Fqcy_all, colo, categ1, yli=c(0, 1), yla="Frequency of unusual CN variants\nin populations", lab=categ)
    x11()
    draw_plot(inter, GLMtab_all2$Fqcy_all, colo, categ1, yli=c(0, 1), yla="Frequency of unusual CN variants\nin populations", lab=categ)
    x11()
    draw_plot(GLMtab_all2$CpDup, GLMtab_all2$Fqcy_all, colo, categ1, yli=c(0, 1), yla="Frequency of unusual CN variants\nin populations", lab=unique(GLMtab_all2$CpDup))


    jpeg(JPG, height=480*2, width=480*2, quality=100, res=72*2)
    draw_plot(inter, GLMtab_all2$Fqcy_all, colo, categ1, yli=c(0, 1), yla="Frequency of unusual CN variants\nin populations", lab=categ)
    text(x=c(2.5, 9), y=c(0.4, 0.9), labels=c("Non\ntruncated", "Truncated"), lheight=1.5)
    abline(v=5.5, lty=2)
    dev.off()
    
    pdf(PDF)
    draw_plot(obs_prob1, colo, categ1, yli=c(0, 1), yla="Proportion of observations completely\nduplicated/deleted for CNV per family")
    text(x=c(2.5, 9), y=c(0.4, 0.9), labels=c("Non\ntruncated", "Truncated"), lheight=1.5)
    abline(v=5.5, lty=2)
    dev.off()

    
    cat("\n")
    print("#### 2) Correalation between exon length and proba of complete duplication/deletion ")
    Complete_CNV <- ifelse(GLMtab_all2$Duplication=="3_CpDup", 1, 0)
    GLMtab_all2[,"Complete_CNV"] <- Complete_CNV
    gplot <- ggplot(GLMtab_all2, aes(x = LnExonLength, y = Complete_CNV, color = trimmed)) + geom_point(alpha=0.7) + labs(colour = "Truncated") + coord_fixed(ratio=2)
    gplot <- gplot + stat_smooth(method = 'glm', family = 'binomial')

    ggsave(PDF2)
    ggsave(JPG2)
    
    #########
    cat("\n")
    print(" #### save results")
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    dummy <- numeric()
    save(dummy,file=outFileName)

	file.remove("Rplots.pdf")
}

cat("dummy logfile\n")


argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);




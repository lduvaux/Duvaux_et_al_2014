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
	load(PREVIOUS_DATA2)
    Genes_Info <- read.delim(INFO_TARGENE_FILE, stringsAsFactors=F)

    # 0) get number of exons per genes
    v_targets <- sapply(rownames(alpha_matrix), collapse_elements)
    # finir de choper le ombre de targets restantes par gene
    
    genes <- unique(GLMtab_all2$Gene)
    N_exon <- sapply(genes, get_exdon_number, Genes_Info)
    LnLengthExon <- sapply(genes, get_length, GLMtab_all2)
    tab <- data.frame(genes, N_exon, LnLengthExon)
    print(cor.test(tab$N_exon, tab$LnLengthExon, method="s"))

    # prepare data for plots
    cat("\n")
    print("#### 1) proba of complete duplication/deletion per truncated/non-truncated gene")
    samp_size1 <- table(with(GLMtab_all2, interaction(Family, trimmed)))
    samp_size1_0 <- table(with(GLMtab_all2, interaction(Family, trimmed, Duplication)))
    obs_prob1.0 <- samp_size1_0/rep(samp_size1, 2)
    obs_prob1 <- obs_prob1.0[(length(obs_prob1.0)/2+1):length(obs_prob1.0)]

    colo <- c("Black", "Blue", "Purple", "DarkGreen")
    categ <- c("Control", "Gr", "Or", "P450")
    categ1 <- paste(categ, " (", samp_size1_0[(length(samp_size1_0)/2+1):length(samp_size1_0)], ")", sep="")

    jpeg(JPG, height=480*2, width=480*2, quality=100, res=72*2)
    draw_plot(obs_prob1, colo, categ1, yli=c(0, 1), yla="Proportion of observations completely\nduplicated/deleted for CNV per family")
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




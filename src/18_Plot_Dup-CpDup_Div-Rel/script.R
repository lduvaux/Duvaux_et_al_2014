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
	load(PREVIOUS_DATA3)
    set.seed(0)
    Genes_Info <- read.delim(INFO_TARGENE_FILE, stringsAsFactors=F)

    # 0) get number of exons per genes
    subtarg <- rownames(alpha_matrix)
    v_targets <- sapply(subtarg, collapse_elements)
    genes <- unique(GLMtab_all2$Gene)
    categ <- as.factor(sapply(as.character(genes), get_elements, what=1))
    N_exon_init <- sapply(genes, get_exdon_number, Genes_Info$NewTargetName)
    N_exon_clean <- sapply(genes, get_exdon_number, unique(v_targets))
    truncated <- !genes%in%NonTrimGenes
    LnLengthExon <- sapply(genes, get_length, GLMtab_all2)
    tab <- data.frame(Genes=genes, LgCoding=exp(LnLengthExon), N_exon=N_exon_init, N_exon_trunc=N_exon_clean, Truncated=truncated, Categ=categ)

#~    boxplot(N_exon_init~truncated)
    inter <- interaction(truncated, categ)
#~    boxplot(N_exon_init~inter)
    
    print(cor.test(tab$N_exon, tab$LgCoding, method="s"))
#~    plot(tab$N_exon, tab$LgCoding, col=as.numeric(tab$Categ))
    p <- ggplot(tab, aes(x = N_exon, y = LgCoding, color = Categ)) + geom_point() + facet_wrap( ~ Categ, ncol=2)
#~    x11()
#~    plot(p)
    bad_Or <- c("Or_g64", "Or_g69", "Or_g55")
    bad_P450 <- c("P450_g39", "P450_g42", "P450_g8", "P450_g10")
    gn_categ <- unique(categ)
    for (i in seq(gn_categ)){
        cate <- gn_categ[i]
        print(cate)
        ind <- categ%in%cate
#~        if (cate =="P450") {
#~            bb <- pmatch(bad_P450, tab$Genes)
#~            ind[bb] <- F
#~        }
        print(table(ind)[2])
#~        plot(tab$N_exon[ind], tab$LgCoding[ind])
        print(cor.test(tab$N_exon[ind], tab$LgCoding[ind], method="s"))
    }

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


    coord_txt_x <- c(2.5, 9)
    coord_txt_y <- c(0.4, 0.9)
    jpeg(JPG, height=480*2, width=480*2, quality=100, res=72*2)
    draw_plot(obs_prob1, colo, categ1, yli=c(0, 1), xla=XLAB1, yla=YLAB1, m_gp=c(5, 1, 0))
    text(x=coord_txt_x, y=coord_txt_y, labels=c("Non\ntruncated", "Truncated"), lheight=1.5)
    abline(v=5.5, lty=2)
    dev.off()
    
    pdf(PDF)
    draw_plot(obs_prob1, colo, categ1, yli=c(0, 1), xla=XLAB1, yla=YLAB1, m_gp=c(5, 1, 0))
    text(x=coord_txt_x, y=coord_txt_y, labels=c("Non\ntruncated", "Truncated"), lheight=1.5)
    abline(v=5.5, lty=2)
    dev.off()

    
    cat("\n")
    print("#### 2) Correalation between exon length and proba of complete duplication/deletion ")
    Complete_CNV <- ifelse(GLMtab_all2$Duplication=="3_CpDup", 1, 0)
    GLMtab_all2[,"Complete_CNV"] <- Complete_CNV
    gplot <- ggplot(GLMtab_all2, aes(x = LnExonLength, y = Complete_CNV, color = trimmed)) + geom_point(alpha=0.7) + labs(colour = "Truncated") + coord_fixed(ratio=2) + labs(x = XLAB2, y=YLAB2)
    gplot <- gplot + stat_smooth(method = 'glm', family = 'binomial') + theme(axis.title=element_text(size=13,face="bold"))

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

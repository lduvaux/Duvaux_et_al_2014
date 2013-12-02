#!/bin/Rscript
library('multicore')
library(nnet)
# library(ggplot2)
# library(reshape2)
#QQ~ library(fields)
library(MuMIn)
library("epicalc")
library(methods)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./params.R")
source("./functions.R")
source("./Models_trimmed_blocks.R")

argv <- commandArgs(TRUE)[1]

# main <- function(argv){

    load(PREVIOUS_DATA)
    set.seed(0)

    ######################
    races <- unique(sapply(colnames(alpha_matrix), function(x) unlist(strsplit(x, "_"))[1]))
    groups <- races ; clusters <- as.list(races) ; names(clusters) <- groups
    print(groups)
    Pre_GLMtab <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
    , Genes_Info, CATEG_FOR_GLM, CompGenes, test_blocks=T, list_groups=clusters)
    Pre_GLMtab$Race <- relevel(Pre_GLMtab$Race, ref="Medicago")
    
    # 1) Check data
    GLMtab_all2 <- set_GLMtab(Pre_GLMtab, NonTrim_only=F, covar="LnExonLength")
    GLMtab_all2$Race <- relevel(GLMtab_all2$Race, ref="Medicago")
    print(sample_size3 <- table(GLMtab_all2$Polymorphism, interaction(GLMtab_all2$Family, GLMtab_all2$trimmed, GLMtab_all2$Race)))
    Draw_pdf(pairs_glm(MODP2, data=GLMtab_all2), PAIRS_ALL_EXON_LENGTH)

		# 1.c) draw histograms per class of interactions
# 	vec_int <- with(Pre_GLMtab, interaction(trimmed, Family))
# 	tab_draw <- cbind(Pre_GLMtab[,c("Polymorphism", "Family", "Race", "GeneLength")], LnGeneLength=log(Pre_GLMtab[,"GeneLength"]), TotExonLength=Pre_GLMtab[,"TotExonLength"], LnExonLength=log(Pre_GLMtab[,"TotExonLength"]))
# 	Draw_pdf(Draw_distrib(tab_draw, Vec_int1=vec_int, Vec_int2=tab_draw$trimmed, numeric_var=4:7, nrow=4), INTERACTION_HIST)

    # 2) Investigate the effect of trimmed on CNV (multinom (nnet))
    print("fit the main model")
    fm1 <- multinom(formula=MOD_ALL1, data = GLMtab_all2)
    print("start all tests")
    test_trimmed <- dredge(fm1, fixed=FIXED_TERMS, m.max=M_MAX)
    
        # 2.1) best
    best_mod <- get_best(test_trimmed, 20)
    test_best <-  get.models(test_trimmed, best_mod)[[1]]
    summ_best <- summary(test_best)
    Output_glm_res(summ_best, GLM_RES_BEST_TRIMMED)
    
    # 3) 


    ######################
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
#     dummy <- numeric()
    save(test_trimmed, sample_size3, file=outFileName)
	print(sort(test_trimmed$AIC))
# }

# if(DEBUG)
# 	traceback(main(argv));
# if(!DEBUG)
# 	main(argv);










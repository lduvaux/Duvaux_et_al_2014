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
source("./Models_trimmed.R")

argv <- commandArgs(TRUE)[1]

# main <- function(argv){

#     # load sources and data
#     if (is.na(argv))
#         stop("Missing argument")

    load(PREVIOUS_DATA)
    set.seed(0)

    models <- Polymorphism ~ LnGeneLength + LnExonLength + Family + trimmed + LnGeneLength * Family + LnExonLength * Family + Family * trimmed
    
    ######################
    
    Pre_GLMtab <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
    , Genes_Info, CATEG_FOR_GLM, CompGenes)
    
    # 1) Check data
        # 1.a) ratio as covariate
    GLMtab_all1 <- set_GLMtab(Pre_GLMtab, NonTrim_only=F, covar="ratioLength")
    Draw_pdf(pairs_glm(MODP1, data=GLMtab_all1), PAIRS_ALL_RATIO)
    
        # 1.b) LnExonLength as covariate
    GLMtab_all2 <- set_GLMtab(Pre_GLMtab, NonTrim_only=F, covar="LnExonLength")
    print(sample_size1 <- table(GLMtab_all2$Polymorphism, GLMtab_all2$Family))
    print(sample_size2 <- table(GLMtab_all2$Family, GLMtab_all2$trimmed))
    print(sample_size3 <- table(GLMtab_all2$Polymorphism, interaction(GLMtab_all2$Family, GLMtab_all2$trimmed)))
    Draw_pdf(pairs_glm(MODP2, data=GLMtab_all2), PAIRS_ALL_EXON_LENGTH)
    
        # 1.c) draw histograms per class of interactions
    vec_int <- with(Pre_GLMtab, interaction(trimmed, Family))
    tab_draw <- cbind(Pre_GLMtab[,c(1:2, 6, 3)], LnGeneLength=log(Pre_GLMtab[,3]), TotExonLength=Pre_GLMtab[,4], LnExonLength=log(Pre_GLMtab[,4]))
    Draw_pdf(Draw_distrib(tab_draw, Vec_int1=vec_int, Vec_int2=tab_draw$trimmed, numeric_var=4:7, nrow=4), INTERACTION_HIST)
    
    # 2) Investigate the effect of trimmed on CNV
    print("fit the main model")
    fm1 <- multinom(formula=models, data = GLMtab_all2)
    print("start all tests")
    test_trimmed <- dredge(fm1, fixed=FIXED_TERMS, m.max=M_MAX)
    # 	model.avg(test_trimmed, subset = delta < 4)
    
        # 2.1) best
    best_mod <- get_best(test_trimmed)
    test_best <-  get.models(test_trimmed, best_mod)[[1]]
    summ_best <- summary(test_best)
    Output_glm_res(summ_best, GLM_RES_BEST_TRIMMED)
    
        # 2.2) without interactions (3rd best)
    test <- multinom(formula=Polymorphism ~ LnGeneLength + LnExonLength + Family + trimmed, data=GLMtab_all2)
    summ_best <- summary(test)
    Output_glm_res(summ_best, GLM_RES_TRIMMED_NoINTERACT)
    
    # 3) draw the results?
    # 	pp_dpolym <- get_prediction(test_best, GLMtab_all2)
    # 	lpp <- melt(pp_dpolym$probaTab, id.vars = c("Family", "trimmed", "LnGeneLength", "LnExonLength"), value.name = "probability")
    # 	lpp$AllInter <- interaction(lpp$Family, lpp$trimmed)
    # 	ggplot(lpp, aes(x = LnGeneLength, y = probability, colour = AllInter)) + geom_line() + facet_grid(variable ~ ., scales = "free")
    #QQ~	Draw_pdf(plotmod_predic(pp_dpolym$heat_map_Complete, GLMtab_all2, "Yes"), HEAT_MAP_COMP)
    #QQ~	Draw_pdf(plotmod_predic(pp_dpolym$heat_map_Partial, GLMtab_all2, "No"), HEAT_MAP_PART)
    
    ######################
    outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    #     dummy <- numeric()
    save(test_trimmed, summ_best, sample_size1, sample_size2, sample_size3, file=outFileName)
    print(sort(test_trimmed$AIC))
# }

# if(DEBUG)
#     traceback(main(argv));
# if(!DEBUG)
#     main(argv);










#!/bin/Rscript
rm(list=ls())
library('multicore')
library(lme4)
library(MuMIn)
library("epicalc")
library(methods)
library(ggplot2)
library(multcomp)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./params.R")
source("./functions.R")
source("./Models_trimmed.R")

#~argv <- commandArgs(TRUE)[1]

# main <- function(argv){

    load(PREVIOUS_DATA)
	set.seed(0)

######################
groups <- names(CLUSTERS)
print(groups)

Pre_GLMtab0 <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
, Genes_Info, CATEG_FOR_GLM, NonTrimGenes, test_blocks=T, list_groups=SUBCLUSTERS, fqcy=F)

Pre_GLMtab0 <- add_Phylog_lvl(Pre_GLMtab0, CLUSTERS)
Pre_GLMtab0 <- add_DupField(Pre_GLMtab0)
#~Pre_GLMtab0 <- add_raceField(Pre_GLMtab0)

res_all_groups <- list()
v_samples <- paste("sample_size_", groups, sep=""); list_samples <- list()


for (ite in seq(length(groups)))
{
    gp <- groups[[ite]]
    Pre_GLMtab <- subset(Pre_GLMtab0, Phylog_lvl==gp)
    
    # 1) Check data
    GLMtab_all1 <- set_GLMtab(Pre_GLMtab, NonTrim_only=F, covar="LnExonLength")
    
    print(sample_size1 <- table(GLMtab_all1$Family, GLMtab_all1$trimmed))
    print(sample_size2 <- table(GLMtab_all1$Family, GLMtab_all1$Dup))
    print(sample_size3 <- table(GLMtab_all1$Family, with(GLMtab_all1, interaction(Dup, trimmed))))
    assign(v_samples[ite], list(sample_size1=sample_size1, sample_size2=sample_size2, sample_size3=sample_size3))
    list_samples[[ite]] <- get(v_samples[ite])
    nompdf <- sub(".pdf", paste("_", gp, ".pdf", sep=""), PAIRS_ALL_EXON_LENGTH)
    Draw_pdf(pairs_glm(MODP2, data=GLMtab_all1), nompdf)

    # 2) Effect of variables on duplication events
    print("Fit model 1")
    fm1 <- glmer(formula=MOD_ALL1, data = GLMtab_all1, family=FAMILY)
    summary(fm1)

    print("Test all terms")
    test_trimmed <- dredge(fm1, fixed=FIXED_TERMS, m.max=M_MAX)
#~    mdl_vag <- model.avg(test_trimmed, subset = delta < 2)
#~    summary(mdl_vag)

        # 2.1) best
    best_mod <- get_best(test_trimmed, Delta=5)
    test_best <-  get.models(test_trimmed, best_mod)[[1]]
    summ_best <- summary(test_best)
#~    nomfil <- sub(".txt", paste("_", gp, ".txt", sep=""), GLM_RES_BEST_TRIMMED)
#~    Output_glm_res(summ_best, nomfil)

    assign(gp, list(test_trimmed=test_trimmed, summ_best=summ_best))
    res_all_groups[[ite]] <- get(gp)

        # 2.2) draw the results?
    nompdf <- sub(".pdf", paste("_", gp, ".pdf", sep=""), BOXPLOTS_PDF)
    Draw_pdf(drw_pred(test_best, GLMtab_all1), nompdf)

    # 3) Effect of variables on complete duplication event
        # 3.1) data
    GLMtab_all2 <- GLMtab_all1[GLMtab_all1$Dup==T,]

        # 3.2) test th emodel
    print("Fit model 2")
    fm2 <- glmer(formula=MOD_ALL2, data = GLMtab_all2, family=FAMILY)
#~    summary(fm2)

    print("Test all terms")
    test_trimmed2 <- dredge(fm2, fixed=FIXED_TERMS, m.max=M_MAX)
#~    mdl_vag <- model.avg(test_trimmed, subset = delta < 2)
#~    summary(mdl_vag)

        # 2.1) best
    best_mod2 <- get_best(test_trimmed2, Delta=5)
    test_best2 <-  get.models(test_trimmed2, best_mod2)[[1]]
    test_best2.1 <-  get.models(test_trimmed2, 2)[[1]]
    summ_best2 <- summary(test_best2)




    
    fac <- with(GLMtab_all1, interaction(trimmed, CpDup, Family))
    boxplot(fitted(fm1)~fac)
    boxplot(GLMtab_all1$Fqcy_all~fac)
    
}
names(res_all_groups) <- groups
names(list_samples) <- groups
######################
outFileName <- argv[1]
ver(sprintf("Saving data to %s",outFileName))
#     dummy <- numeric()
save(res_all_groups, list_samples, file=outFileName)
print(sort(test_trimmed$AIC))
# }

# if(DEBUG)
# 	traceback(main(argv));
# if(!DEBUG)
# 	main(argv);










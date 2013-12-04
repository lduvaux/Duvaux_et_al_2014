#!/bin/Rscript
rm(list=ls())
library('multicore')
library(lme4)
library(MuMIn)
library("epicalc")
library(methods)

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
, Genes_Info, CATEG_FOR_GLM, NonTrimGenes, test_blocks=T, list_groups=SUBCLUSTERS, fqcy=T)

Pre_GLMtab0 <- add_Phylog_lvl(Pre_GLMtab0, CLUSTERS)
Pre_GLMtab0 <- add_CpDupField(Pre_GLMtab0)
Pre_GLMtab0 <- add_raceField(Pre_GLMtab0)

res_all_groups <- list()
v_samples <- paste("sample_size_", groups, sep=""); list_samples <- list()


for (ite in seq(length(groups)))
{
    gp <- groups[[ite]]
    Pre_GLMtab <- subset(Pre_GLMtab0, Phylog_lvl==gp)
    
    # 1) Check data
    GLMtab_all2 <- set_GLMtab(Pre_GLMtab, NonTrim_only=F, covar="LnExonLength")

    ####### => start again from here
    
    print(sample_size1 <- table(GLMtab_all2$Family, GLMtab_all2$trimmed))
    print(sample_size2 <- table(GLMtab_all2$Family, GLMtab_all2$CpDup))
    print(sample_size3 <- table(GLMtab_all2$Family, interaction(GLMtab_all2$trimmed, GLMtab_all2$CpDup)))
    assign(v_samples[ite], list(sample_size1=sample_size1, sample_size2=sample_size2, sample_size3=sample_size3))
    list_samples[[ite]] <- get(v_samples[ite])
    nompdf <- sub(".pdf", paste("_", gp, ".pdf", sep=""), PAIRS_ALL_EXON_LENGTH)
    Draw_pdf(pairs_glm(MODP2, data=GLMtab_all2), nompdf)

        # 1.c) draw histograms per class of interactions


    # 2) Investigate the effect of trimmed on CNV
    print("Fit the main model")
    fm1 <- glmer(formula=MOD_ALL1, data = GLMtab_all2, family=FAMILY)

    fac <- with(GLMtab_all2, interaction(trimmed, CpDup, Family))
    boxplot(fitted(fm1)~fac)
    boxplot(GLMtab_all2$Fqcy_all~fac)
    
    print("Test all terms")
    test_trimmed <- dredge(fm1, fixed=FIXED_TERMS, m.max=M_MAX)
    model.avg(test_trimmed, subset = delta < 4)









#### start back here

        # 2.1) best
    best_mod <- get_best(test_trimmed)
    test_best <-  get.models(test_trimmed, best_mod)[[1]]
    summ_best <- summary(test_best)
    nomfil <- sub(".txt", paste("_", gp, ".txt", sep=""), GLM_RES_BEST_TRIMMED)
    Output_glm_res(summ_best, nomfil)

    assign(gp, list(test_trimmed=test_trimmed, summ_best=summ_best))
    res_all_groups[[ite]] <- get(gp)
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










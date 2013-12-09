#!/bin/Rscript
library('multicore')
library(nnet)
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

load(PREVIOUS_DATA)
set.seed(0)

######################
groups <- names(CLUSTERS)
print(groups)
Pre_GLMtab0 <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
, Genes_Info, CATEG_FOR_GLM, NonTrimGenes, test_blocks=T, list_groups=CLUSTERS)
res_all_groups <- list()
v_samples <- paste("sample_size_", groups, sep=""); list_samples <- list()
for (ite in seq(length(groups)))
{
    gp <- groups[[ite]]
    Pre_GLMtab <- subset(Pre_GLMtab0, Race==gp)
    
    # 1) Check data
    GLMtab_all2 <- set_GLMtab(Pre_GLMtab, NonTrim_only=F, covar="LnExonLength")
    print(sample_size1 <- table(GLMtab_all2$Polymorphism, GLMtab_all2$Family))
    print(sample_size2 <- table(GLMtab_all2$Family, GLMtab_all2$trimmed))
    print(sample_size3 <- table(GLMtab_all2$Polymorphism, interaction(GLMtab_all2$Family, GLMtab_all2$trimmed)))
    assign(v_samples[ite], list(sample_size1=sample_size1, sample_size2=sample_size2, sample_size3=sample_size3))
    list_samples[[ite]] <- get(v_samples[ite])
    nompdf <- sub(".pdf", paste("_", gp, ".pdf", sep=""), PAIRS_ALL_EXON_LENGTH)
    Draw_pdf(pairs_glm(MODP2, data=GLMtab_all2), nompdf)

        # 1.c) draw histograms per class of interactions
    vec_int <- with(Pre_GLMtab, interaction(trimmed, Family))
    tab_draw <- cbind(Pre_GLMtab[,c("Polymorphism", "Family", "Race", "GeneLength")], LnGeneLength=log(Pre_GLMtab[,"GeneLength"]), TotExonLength=Pre_GLMtab[,"TotExonLength"], LnExonLength=log(Pre_GLMtab[,"TotExonLength"]))
    nompdf <- sub(".pdf", paste("_", gp, ".pdf", sep=""), INTERACTION_HIST)
    Draw_pdf(Draw_distrib(tab_draw, Vec_int1=vec_int, Vec_int2=tab_draw$trimmed, numeric_var=4:7, nrow=4), nompdf)

    # 2) Investigate the effect of trimmed on CNV
    print("Fit the main model")
    fm1 <- multinom(formula=MOD_ALL1, data = GLMtab_all2)
    MOD_ALL2 <- Polymorphism ~ LnGeneLength + LnExonLength + Family + trimmed + Family * trimmed

    print("Test all terms")
    test_trimmed <- dredge(fm1, fixed=FIXED_TERMS, m.max=M_MAX)
# 	model.avg(test_trimmed, subset = delta < 4)

        # 2.1) best
    best_mod <- get_best_nnet_multinom(test_trimmed)
    test_best <-  get.models(test_trimmed, best_mod)[[1]]
    summ_best <- summary(test_best)
    nomfil <- sub(".txt", paste("_", gp, ".txt", sep=""), GLM_RES_BEST_TRIMMED)
    Output_glm_res(summ_best, nomfil)
    
    # 3) draw the results?
    nompdf <- sub(".pdf", paste("_", gp, ".pdf", sep=""), BOXPLOTS_PDF)
    Draw_pdf(drw_pred(test_best, GLMtab_all2), nompdf)
    
    
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










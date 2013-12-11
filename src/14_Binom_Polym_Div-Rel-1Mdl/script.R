#!/bin/Rscript
rm(list=ls())
#~library('multicore')
library(lme4)
library(MuMIn)
library(methods)
library(ggplot2)

# 0) set up parallel dredge
library(parallel)
    # Set up the cluster
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))

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
cat("
##########################################
14_Binom_Polym_Div-Rel-1Mdl: test the effect of gene features on Pr(polymorphic)
##########################################\n")

# 14.1) Set up GLM data table
cat("\n\n# 14.1) Set up GLM data table\n")
Pre_GLMtab0 <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
, Genes_Info, CATEG_FOR_GLM, NonTrimGenes, test_blocks=T, list_groups=SUBCLUSTERS, fqcy=F)

Pre_GLMtab0 <- add_Phylog_lvl(Pre_GLMtab0, CLUSTERS)
Pre_GLMtab0 <- add_PolymField(Pre_GLMtab0)
Pre_GLMtab0 <- add_GeneField(Pre_GLMtab0)

print(head(Pre_GLMtab0))

##### main script
# 14.2) Check data
cat("\n\n     # 14.2) Check data distrisbution\n")
GLMtab_all1 <- set_GLMtab(Pre_GLMtab0, NonTrim_only=F, covar="LnExonLength")

cat("\n             # Dimension of GLMtab_all1\n")
print(dim(GLMtab_all1))

print(sample_size1 <- table(GLMtab_all1$Family, GLMtab_all1$trimmed))
print(sample_size2 <- table(GLMtab_all1$Family, GLMtab_all1$Polymorphic))
print(sample_size3 <- table(GLMtab_all1$Family, with(GLMtab_all1, interaction(Polymorphic, trimmed))))

list_samples <- list(sample_size1=sample_size1, sample_size2=sample_size2, sample_size3=sample_size3)
Draw_pdf(pairs_glm(MODP2, data=GLMtab_all1), PAIRS_ALL_EXON_LENGTH)

# 14.3) Effect of variables on the proba to be polymorphic
    # 14.3.1) fit the models
cat("\n\n     # 14.3) Effect of variables on Pr(Polymorphic)\n")
cat("\n         # 14.3.1) Fit the maximal model\n")
fm1 <- glmer(formula=MOD_ALL1, data = GLMtab_all1, family=FAMILY, control=glmerControl(optCtrl=list(maxfun=15000)))
sum_fm1 <- summary(fm1)
print(sum_fm1, corr=F)
nomfil <- GLM_POL_MAX
output_glm(sum_fm1, nomfil)

    # 14.3.2) Fit al other models & model averaging
cat("\n         # 14.3.2) Fit al other models & model averaging\n")
cat("\n             # Test all terms\n")
clusterExport(clust, c("GLMtab_all1", "FAMILY"))
clusterEvalQ(clust, library(lme4))
print(system.time(test_trimmed <- pdredge(fm1, cluster=clust, fixed=FIXED_TERMS1, m.max=M_MAX)))

print(test_trimmed)

nomfil <-  GLM_POL_DREDGE
cat(capture.output(test_trimmed), file=nomfil, sep="\n")

cat("\n\n             # get deviance of the best model\n")
best_fma <- formula(attributes(test_trimmed)$calls[[1]])
best_fm <- glmer(formula=best_fma, data = GLMtab_all1, family=FAMILY, control=glmerControl(optCtrl=list(maxfun=15000)))
sum_best_fm <- summary(best_fm)
print(sum_best_fm, corr=F)
output_glm(sum_best_fm, GLM_POL_BEST)

cat("\n             # model averaging\n", sep="")
mdl_avg <- model.avg(test_trimmed, subset = delta < DELTA1)
sum_avg1 <- summary(mdl_avg)
print(sum_avg1) # save in a file
nomfil <-  GLM_POL_AVG
cat(capture.output(sum_avg1, nomfil), file=nomfil, sep="\n")

    # 14.3.3) draw results
#~cat("\n         # 14.3.3) Draw predicted Pr(Polymorphic)\n")
#~nompdf <- POL_PDF
#~Draw_pdf(drw_pred(test_trimmed, GLMtab_all1, DELTA1), nompdf)

# 14.5) record results
GLMPolymorphic <- list(max_mdl_Pol=fm1, all_mdl_Pol=test_trimmed, mdl_avg_Pol=mdl_avg)

#### end of the script
outFileName <- argv[1]
ver(sprintf("Saving data to %s",outFileName))
#     dummy <- numeric()
save(GLMPolymorphic, list_samples, file=outFileName)
# }

# if(DEBUG)
# 	traceback(main(argv));
# if(!DEBUG)
# 	main(argv);










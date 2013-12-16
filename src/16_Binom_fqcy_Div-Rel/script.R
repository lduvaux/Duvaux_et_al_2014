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

argv <- commandArgs(TRUE)[1]

# main <- function(argv){

load(PREVIOUS_DATA)
set.seed(0)

######################
cat("
##########################################
16_Binom_fqcy_Div-Rel: test the effect of gene features on Pr(Dup_fqcy)
##########################################\n")

# 16.1) Set up GLM data table
cat("\n\n# 16.1) Set up GLM data table\n")
Pre_GLMtab0 <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
, Genes_Info, CATEG_FOR_GLM, NonTrimGenes, test_blocks=T, list_groups=SUBCLUSTERS, fqcy=T)

Pre_GLMtab0 <- add_Phylog_lvl(Pre_GLMtab0, CLUSTERS)
Pre_GLMtab0 <- add_PolymField_fqcy(Pre_GLMtab0)
Pre_GLMtab0 <- add_GeneField(Pre_GLMtab0)
Pre_GLMtab0 <- add_CpDupField_fqcy(Pre_GLMtab0)

print(head(Pre_GLMtab0))

##### main script
# 16.2) Check data
cat("\n\n     # 16.2) Check data distrisbution\n")
GLMtab_all1 <- set_GLMtab_fqcy(Pre_GLMtab0, NonTrim_only=F)

    # 16.2.1) data
cat("\n         # 16.2.1) Remove non polymorphic genes\n")
GLMtab_all2 <- GLMtab_all1[GLMtab_all1$Polymorphic==T,]

print(sample_size1 <- table(GLMtab_all2$Family, GLMtab_all2$trimmed))
print(sample_size2 <- table(GLMtab_all2$Family, GLMtab_all2$CpDup))
print(sample_size3 <- table(GLMtab_all2$Family, with(GLMtab_all2, interaction(CpDup, trimmed))))

list_samples <- list(sample_size1=sample_size1, sample_size2=sample_size2, sample_size3=sample_size3)
Draw_pdf(pairs_glm(MODP2, data=GLMtab_all2), PAIRS_ALL_EXON_LENGTH)

# 16.3) Effect of variables on complete duplication events
cat("\n\n     # 16.3) Effect of variables on complete duplication events\n")

    # 16.3.1)  Fit the maximal model
cat("\n         # 16.3.1) Fit the maximal model\n")
cat("\n             # Dimension of GLMtab_all2\n")
print(dim(GLMtab_all2))
cat("\n")

fm2 <- glmer(formula=MOD_ALL2, data = GLMtab_all2, family=FAMILY, control=glmerControl(optCtrl=list(maxfun=15000)))
sum_fm2 <- summary(fm2)
print(sum_fm2, corr=F)
nomfil <- GLM_DUPFQCY_MAX
output_glm(sum_fm2, nomfil)

    # 16.3.2) Fit all other models & model averaging
cat("\n         # 16.3.2) Fit al other models & model averaging\n")
cat("\n             # Test all terms\n")
clusterExport(clust, c("GLMtab_all2", "FAMILY"))
clusterEvalQ(clust, library(lme4))
test_trimmed2 <- pdredge(fm2, cluster=clust, fixed=FIXED_TERMS2, m.max=M_MAX)

print(test_trimmed2)

nomfil <- GLM_DUPFQCY_DREDGE
cat(capture.output(test_trimmed2), file=nomfil, sep="\n")

cat("\n\n             # get deviance of the best model\n")
best_fma2 <- formula(attributes(test_trimmed2)$calls[[1]])
best_fm2 <- glmer(formula=best_fma2, data = GLMtab_all2, family=FAMILY, control=glmerControl(optCtrl=list(maxfun=15000)))
sum_best_fm2 <- summary(best_fm2)
print(sum_best_fm2, corr=F)
output_glm(sum_best_fm2, GLM_DUPFQCY_BEST)

cat("\n             # model averaging\n", sep="")
mdl_avg2 <- model.avg(test_trimmed2, subset = delta < DELTA2)
sum_avg2 <- summary(mdl_avg2)
print(sum_avg2)
nomfil <- GLM_DUPFQCY_AVG
cat(capture.output(sum_avg2, nomfil), file=nomfil, sep="\n")

    # 16.3.3) draw the results
#~cat("\n         # 16.3.4) Draw predicted probability of complete duplication\n")
#~nompdf <-  DUPFQCY_PDF
#~Draw_pdf(drw_pred(test_trimmed2, GLMtab_all2, DELTA2[ite], draw_CpDup=T), nompdf)

# 16.4) record results
Res_CpDup <- list(max_mdl_CpDup=fm2, all_mdl_CpDup=test_trimmed2, mdl_avg_CpDup=mdl_avg2)


#### end of the script
outFileName <- argv[1]
ver(sprintf("Saving data to %s",outFileName))
#     dummy <- numeric()
save(Res_CpDup, list_samples, file=outFileName)
# }

# if(DEBUG)
# 	traceback(main(argv));
# if(!DEBUG)
# 	main(argv);










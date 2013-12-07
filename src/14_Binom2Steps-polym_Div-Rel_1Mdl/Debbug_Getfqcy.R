library('multicore')
library(lme4)
library(MuMIn)
library(methods)
library(ggplot2)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./params.R")
source("./functions.R")
source("./Models_trimmed.R")

load(PREVIOUS_DATA)
set.seed(0)


Pre_GLMtab0 <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
, Genes_Info, CATEG_FOR_GLM, NonTrimGenes, test_blocks=T, list_groups=SUBCLUSTERS, fqcy=F)
Pre_GLMtab0 <- add_Phylog_lvl(Pre_GLMtab0, CLUSTERS)
Pre_GLMtab0 <- add_DupField(Pre_GLMtab0)

vec_polym <- Pre_GLMtab0$Polymorphism
table(vec_polym)
sum(table(vec_polym))

GLMtab_all1 <- set_GLMtab(Pre_GLMtab0, NonTrim_only=F, covar="LnExonLength")
table(GLMtab_all1$Polymorphism)


Pre_GLMtab0_2 <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
, Genes_Info, CATEG_FOR_GLM, NonTrimGenes, test_blocks=T, list_groups=SUBCLUSTERS, fqcy=T)
vec_fqcy <- Pre_GLMtab0_2$Fqcy_all
length(which(vec_fqcy!=0))


vec_fqcy <- Pre_GLMtab0_2$Fqcy_all

v_ind <- !((vec_polym=="1_NoDup" & vec_fqcy==0)|((vec_polym=="2_PtDup"|vec_polym=="3_CpDup") & vec_fqcy!=0))
table(v_ind)

Pre_GLMtab0[v_ind,]
Pre_GLMtab0_2[v_ind,]

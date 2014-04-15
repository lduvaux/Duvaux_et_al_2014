#!/bin/Rscript
rm(list=ls())
library(ape)
library(randomForest)
library(parallel)

source("params.R")
load(PREVIOUS_DATA)
ll <- ls()
ll <- ll[-which(ll=="x")]
rm(list=ll)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")
source("../utils/randomForest_helperFuns.R")

source("params.R")  # needed as the previous will load former params!
source("./functions.R")

set.seed(0)
races <- sapply(rownames(x), get_elements)
vec_vst <- apply(x, 2, get_Vst, races)
hist(vec_vst, 50)

length(vec_vst[vec_vst==1])/length(vec_vst)

gd <- which(vec_vest==1)
ftab <- x[,gd]
write.table(cbind(rownames(ftab),ftab), col.names=T, row.names=F, sep="\t", quote=F, file="tab_Vst1.csv")

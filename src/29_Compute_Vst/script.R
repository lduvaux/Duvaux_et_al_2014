#!/bin/Rscript
rm(list=ls())
library(ape)
library(randomForest)
library(parallel)


cat("\n")
print(" #### 29.1) load data, sources, parameters and functions")
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
tab_targ <- read.delim(TARG, stringsAsFactors=F)
tab_cnv <- t(x)
inval_rows <- PrePro_findMeaninglessRows(tab_cnv)
tab_cnv <- tab_cnv[!inval_rows,]

cat("\n")
print(" #### 29.2) compute Vst")
    # 29.2.1) set races and targ names
races <- sapply(colnames(tab_cnv), get_elements)
targets <- sapply(rownames(tab_cnv), collapse_elements)

    # 29.2.2) compute Vst
vec_vst <- apply(tab_cnv, 1, get_Vst, races)
        # write a subsample in a table
gd <- which(vec_vst==1)
ftab <- tab_cnv[gd,]
write.table(cbind(rownames(ftab),ftab), col.names=T, row.names=F, sep="\t", quote=F, file=SUBSAMP1)

vec_vst_log2 <- apply(tab_cnv, 1, get_Vst, races, log_2=T)
#~big <- which(abs(vec_vst-vec_vst_log2)>0.5)

cat("\n")
print(" #### 29.2) compute Vst")
vec_fam <- sapply(rownames(tab_cnv), get_elements)
families <- unique(vec_fam)
list_ind <- lapply(families, function(x) which(vec_fam==x))
names(list_ind) <- families

cat("\n")
print(" #### 29.3) plot histograms")
pdf(PLOT_HIST_VST)
layout(matrix(1:4, nrow=2, ncol=2, byrow=T))
famm <- c("Gr", "Or", "P450")
colos <- c("blue", "purple", "green")
sapply(1:length(famm), function(x) plot_dble_hist(famm[x], list_ind, color=colos[x], alpha=70))
dev.off()

#~cat("\n")
#~print(" #### 29.4) keep the highest Vst per contig")
#~vec_cont <- tab_targ$contigV2[sapply(targets, function(x) which(x==tab_targ$NewTargetName))]
#~contigs <- unique(vec_cont)

#~vec_maxVst <- get_all_maxVst(vec_vst, vec_cont)
#~vec_maxVst_log2 <- get_all_maxVst(vec_vst_log2, vec_cont)





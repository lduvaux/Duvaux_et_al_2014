#!/bin/Rscript
rm(list=ls())
library(tseries)

source("./params.R")
load(PREVIOUS_DATA)
load(PREVIOUS_DATA1)
source("./params.R")
load(PREVIOUS_DATA2)
load(PREVIOUS_DATA3)
ll <- ls()
gd <- match(c("contig_rf", "GLMtab_all1", "alpha_matrix", "distr_rk_rdom_LD_mat"), ll)
ll <- ll[-gd]
rm(list=ll)
source("./params.R")
source("./functions.R")

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

set.seed(0)

cat("\n")
print("#####  1) write list of most informative genes (1 per contig)")
tab_best_obs_raw <- contig_rf$importance
tab_best_obs <- tab_best_obs_raw[order(tab_best_obs_raw[,"MeanDecreaseGini"], decreasing=T),]
Global_rank <- 1:nrow(tab_best_obs)
Ranking_score <- nrow(tab_best_obs):1
tab_best_obs <- data.frame(Subtarget=rownames(tab_best_obs), Global_rank, Ranking_score, tab_best_obs)
write.table(tab_best_obs_raw, file=FIL, sep="\t", quote=F, row.names=F)


cat("\n")
print("#####  2) Detect what Gr discriminate best")
ind_Gr <- grep("Gr_", rownames(tab_best_obs_raw))
tab_best_obs_Gr <- tab_best_obs_raw[ind_Gr,]
v_races <- colnames(tab_best_obs_Gr)[1:8]
best_discrim <- apply(tab_best_obs_Gr[,1:8], 1, which.max)
race_discrim <- v_races[best_discrim]
count_race_discrim <- table(race_discrim)


cat("\n")
print("#####  3) Estimate the number and ranks of CDD genes in the best genes")
gn_best_obs <- sapply(rownames(tab_best_obs), collapse_elements, what=1:2)
categ <- sapply(rownames(tab_best_obs), get_elements, what=1)
data_CDD_obs <- get_info_rking(tab_best_obs, GLMtab_all1)


cat("\n")
print("#####  4) Homemade test to account for sampling scheme (best disc race only)")
N_CDD <- length(which(v_CDD_best))
sum_rk_obs <- data_CDD_obs$sum_rank_CDD
simul_sum_rk <- replicate(5000, sum(sample(111:1, N_CDD)))
P_val_rk_obs <- get_pval(sum_rk_obs, simul_sum_rk, two_sided=F)

hist(simul_sum_rk)
abline(v=sum_rk_obs, col="red")


cat("\n")
print("##### save data")
dummy <- NULL
argv <- commandArgs(TRUE)[1]
outFileName <- argv[1]
ver(sprintf("Saving data to %s",outFileName))
save.image(file=outFileName)


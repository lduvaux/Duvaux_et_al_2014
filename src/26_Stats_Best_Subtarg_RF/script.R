#!/bin/Rscript
rm(list=ls())

source("./params.R")
load(PREVIOUS_DATA)
load(PREVIOUS_DATA1)
source("./params.R")
load(PREVIOUS_DATA2)
load(PREVIOUS_DATA3)
ll <- ls()
gd <- match(c("contig_rf", "GLMtab_all1", "alpha_matrix", "distr_rk_rdom_LD_mat", "infoTargGene"), ll)
ll <- ll[-gd]
rm(list=ll)
source("./params.R")
source("./functions.R")

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

set.seed(0)


cat("\n")
print("#####  0) stats on truncated genes (Gene length, intron length...) for duplicated genes only")
ind <- seq(1, nrow(GLMtab_all1), by=8)
GLMtab_1Pgen0 <- GLMtab_all1[ind,]
GLMtab_1Pgen <- subset(GLMtab_1Pgen0, Polymorphic==T)
table(GLMtab_1Pgen$trimmed)

pdf(PDF0)
boxplot(GLMtab_1Pgen$LnIntronLength ~ GLMtab_1Pgen$trimmed)
print(t.test(GLMtab_1Pgen$LnIntronLength ~ GLMtab_1Pgen$trimmed))
print(kruskal.test(GLMtab_1Pgen$LnIntronLength ~ GLMtab_1Pgen$trimmed))

boxplot(GLMtab_1Pgen$LnGeneLength ~ GLMtab_1Pgen$trimmed)
print(t.test(GLMtab_1Pgen$LnGeneLength ~ GLMtab_1Pgen$trimmed))
print(kruskal.test(GLMtab_1Pgen$LnGeneLength ~ GLMtab_1Pgen$trimmed))
dev.off()

cat("\n")
print("#####  1) write list of most informative genes (1 per contig)")
tab_best_obs_raw <- round(contig_rf$importance, 5)
ncoll <- ncol(tab_best_obs_raw)
tab_best_obs_raw <- tab_best_obs_raw[,c(ncoll,1:(ncoll-1))]
tab_best_obs <- tab_best_obs_raw[order(tab_best_obs_raw[,"MeanDecreaseGini"], decreasing=T),]
Global_rank <- 1:nrow(tab_best_obs)
Ranking_score <- nrow(tab_best_obs):1
tab_best_obs <- data.frame(Subtarget=rownames(tab_best_obs), Global_rank, Ranking_score, tab_best_obs)

nom_targ <- sapply(as.character(tab_best_obs$Subtarget), collapse_elements)
Gene <- sapply(nom_targ, get_real_name, infoTargGene)


cat("\n")
print("#####  2) Detect what Gr discriminate best")
ind_Gr <- grep("Gr_", rownames(tab_best_obs_raw))
tab_best_obs_Gr <- tab_best_obs_raw[ind_Gr,]
indiv <- colnames(alpha_matrix)
v_races <- unique(sapply(indiv, get_elements))
best_discrim <- apply(tab_best_obs_Gr[,v_races], 1, which.max)
race_discrim <- v_races[best_discrim]
count_race_discrim <- table(race_discrim)


cat("\n")
print("#####  3) Estimate the number and ranks of CDD genes in the best genes")
    # 3.1) gather information
gn_best_obs <- sapply(rownames(tab_best_obs), collapse_elements, what=1:2)
categ <- sapply(rownames(tab_best_obs), get_elements, what=1)
data_CDD_obs <- get_info_rking(tab_best_obs, GLMtab_all1, v_races, varia=VARIA)
tab_best_obs <- data.frame(Gene, data_CDD_obs$info_CDD, tab_best_obs)

    # 3.2) add CN info to table
        # 3.2.1) remove bad individuals and round alpha mat
ind <- PrePro_findIndex(c(BAD_CYTISUS, RACE_UNKNOWN), colnames(alpha_matrix))
alpha_matrix2 <- alpha_matrix[,-ind]
alpha_matrix2 <- PrePro_roundToZeroFive(alpha_matrix2)
        # 3.2.2) define list of races
indiv <- colnames(alpha_matrix2)
races <- unique(sapply(indiv, get_elements))
l_races <- lapply(paste(races, "_", sep=""), grep, indiv)
names(l_races) <- races
        # 3.2.3) compute median
ind_best <- pmatch(tab_best_obs$Subtarget, rownames(alpha_matrix2))
alpha_mat_best <- alpha_matrix2[ind_best,]
tab_med <- t(apply(alpha_mat_best, 1, get_CNmedian_race, l_races))

ind_tab_best_obs <- colnames(tab_best_obs)%in%races
nonom <- colnames(tab_best_obs)[ind_tab_best_obs]
coll_nam <- paste(c("Med-CN", "RF-score"), rep(colnames(tab_med), each=2), sep="_")
tab_temp <- as.data.frame(matrix(data=NA, ncol=16, nrow=nrow(tab_med), dimnames=list(rownames(tab_med), c(coll_nam))))
tab_temp[seq(1, 16, by=2)] <- tab_med
tab_temp[seq(2, 16, by=2)] <- tab_best_obs[ind_tab_best_obs]

Family <- sapply(as.character(tab_best_obs$Subtarget), get_elements)
Most_discriminated <- races[apply(tab_best_obs[,races], 1, which.max)]
print(table(Most_discriminated))
tab_f <- data.frame(Gene=tab_best_obs[,1], Family, tab_best_obs[,2:7], Most_discriminated, tab_temp, tab_best_obs[,ncol(tab_best_obs)])

jpeg(JPG1, height=480*2, width=480*2, quality=100, res=72*2)
par(mar=c(5,5.5,4,2))
plot(tab_f$MeanDecreaseGini, xlab="Subtarget rank", ylab="Subtarget information\n(Mean decrease Gini)")
abline(v=c(VARIA), lty=3)
dev.off()

pdf(PDF1)
par(mar=c(5,5.5,4,2))
plot(tab_f$MeanDecreaseGini, xlab="Subtarget rank", ylab="Subtarget information\n(Mean decrease Gini)")
abline(v=c(VARIA), lty=3)
dev.off()


write.table(tab_f, file=FIL, sep="\t", quote=F, row.names=F)


cat("\n")
print("#####  4) Over-representation of CDD genes in disciminating targets: homemade test to account for sampling scheme (best disc race only)")
    # 4.1) test the sum of the rank
N_CDD <- data_CDD_obs$N_CDD
sum_rk_obs <- data_CDD_obs$sum_rank_CDD
simul_sum_rk <- replicate(5000, sum(sample(nrow(tab_best_obs):1, N_CDD)))

print("## results of the test")
print(P_val_rk_obs <- get_pval(sum_rk_obs, simul_sum_rk, two_sided=F))

hist(simul_sum_rk)
abline(v=sum_rk_obs, col="red")

    # 4.2) test the number of CDD in the top varia
N_CDD_varia <- data_CDD_obs$N_CDD_varia
simul_N_in_varia <- replicate(5000, length(
    which(sample(1:nrow(tab_best_obs), N_CDD)<=VARIA)
    ))
P_val_N_in_varia <- get_pval(N_CDD_varia, simul_N_in_varia, two_sided=F)

hist(simul_N_in_varia)
#~abline(v=sum_rk_obs_varia, col="red")

    
cat("\n")
print("##### save data")
dummy <- NULL
argv <- commandArgs(TRUE)[1]
outFileName <- argv[1]
ver(sprintf("Saving data to %s",outFileName))
save.image(file=outFileName)

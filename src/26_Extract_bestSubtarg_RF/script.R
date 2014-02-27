#!/bin/Rscript
rm(list=ls())

source("./params.R")
load(PREVIOUS_DATA)
load(PREVIOUS_DATA0)
ll <- ls()
gd <- match(c("contig_rf", "GLMtab_all2"), ll)
ll <- ll[-gd]
rm(list=ll)
source("./params.R")

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

set.seed(0)

cat("\n")
print("#####  1) write list of most informative genes (1 per contig)")
tab <- contig_rf$importance
tab0 <- tab[order(tab[,"MeanDecreaseGini"], decreasing=T),]
Global_rank <- 1:nrow(tab)
Ranking_score <- nrow(tab):1
tab <- data.frame(Subtarget=rownames(tab0), Global_rank, Ranking_score, tab0)

write.table(tab, file=FIL, sep="\t", quote=F, row.names=F)


cat("\n")
print("#####  2) Detect what Gr discriminate best")
ind_Gr <- grep("Gr_", rownames(tab0))
tab1 <- tab0[ind_Gr,]
v_races <- colnames(tab1)[1:8]
best_discrim <- apply(tab1[,1:8], 1, which.max)
race_discrim <- v_races[best_discrim]
count_race_discrim <- table(race_discrim)


cat("\n")
print("#####  3) Check if best genes are complete")
gn_best <- sapply(rownames(tab0), collapse_elements, what=1:2)
categ <- sapply(rownames(tab0), get_elements, what=1)

v_genes <- NULL
v_test <-NULL
v_rank <- NULL
v_sc_rank <- NULL
v_categ <- NULL
for (i in 1:length(gn_best)){

    gn <- gn_best[i]
    subtarg <- rownames(tab0)[i]
    ind <- which(gn==GLMtab_all2$Gene)

    if (length(ind)>0) {
    v_genes <- c(v_genes, gn)
    v_categ <- c(v_categ, categ[i])

    
    GLM3 <- GLMtab_all2[ind,]
    test <- "3_CpDup"%in%GLM3$Duplication
    v_test <- c(v_test, test)

    sc_rk <- tab[i,"Ranking_score"]
    v_sc_rank <- c(v_sc_rank, sc_rk)
    
    rk <- tab[i,"Global_rank"]
    v_rank <- c(v_rank, rk)
    
    } else {
        print(paste("No information for marker", gn))
    }
    if (i%%100==0) print(i)
}
mat <- data.frame(Gene=v_genes, Family=v_categ, Rank=v_rank, Ranking_score=v_sc_rank, CDD=v_test)

table(mat$CDD)
boxplot(mat$Ranking_score~mat$CDD)
boxplot(mat$Ranking_score~mat$Family)

kd <- interaction(mat$Family, mat$CDD)
print(tt <- table(kd))
kd0 <- sapply(names(tt), get_elements, what=1, sep="\\.")
boxplot(mat$Ranking_score~interaction(mat$Family, mat$CDD), names=kd0)


cat("\n")
print("##### save data")
dummy <- NULL
argv <- commandArgs(TRUE)[1]
outFileName <- argv[1]
ver(sprintf("Saving data to %s",outFileName))
save.image(file=outFileName)


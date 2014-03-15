#!/bin/Rscript
rm(list=ls())
library(ape)
library(randomForest)
library(parallel)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("../utils/randomForest_helperFuns.R")
source("./params.R")
source("./functions.R")

load(PREVIOUS_DATA)
set.seed(0)

races <- PrePro_fetchRaces(RAW_DATA, CTRL_GUYS , BAD_GUYS)
races[races=="Medicago_ctrl"] <- "Medicago"
races_uniq <- PrePro_fetchUniqRaces(races)


pure_indiv <- names(y)
ind_pure_indiv <- PrePro_findIndex(pure_indiv, indiv)


    cat("\n")
print("###### 1) run a RF by removing all noisy loci ######")
        # 1.1) extract locus importance from the biggest RF
temp_tab <- sapply(gradNtree_lrf_unsuperv, function(x) x$importance[,"MeanDecreaseGini"])
gini_NmaxTrees <- temp_tab[, ncol(temp_tab)]
mat_imp_NmaxTrees <- gradNtree_lrf_unsuperv[[length(gradNtree_lrf_unsuperv)]]$importance

    # 1.2) remove uninformative variables (i.e. baits)
bad <- which(gini_NmaxTrees==0)

gini_NmaxTrees1 <- gini_NmaxTrees[-bad]
new_x <- x[,-bad]

    # 1.3) keep only the better bait of each informative exon
data_PerExons <- Keep_MaxBait(gini_NmaxTrees1, new_x)

    # 1.4) run the RF
new_x1 <- data_PerExons$new_xf
new_xtest1 <- xtest[,PrePro_findIndex(colnames(new_x1), colnames(xtest))]

print(sprintf("Growing random forest of %i trees...", NTREES))
exon_rf <- randomForest(x=new_x1, y=y, xtest=new_xtest1, ytest=ytest, ntree=NTREES, proximity = TRUE, importance=T, sampsize=SAMP_SIZE2)
print(exon_rf)

    # 1.5) extract best baits per exon
gini_exon_rf0 <- exon_rf$importance[,"MeanDecreaseGini"]
N_cate_rfEx <- table(sapply(names(gini_exon_rf0), get_elements))

gini_exon_rf <- sort(gini_exon_rf0, decreasing=T)
gini_exon_rf_ind <- order(gini_exon_rf0, decreasing=T)
best50_ExonPromot <- names(gini_exon_rf[1:50])
tab_best50 <- (exon_rf$importance[gini_exon_rf_ind,])[1:50,]
geneOfbest50 <- unique(sapply(best50_ExonPromot, get_elements, sep="\\."))

tem_tab <- cbind(rownames(tab_best50), round(tab_best50,4))
write.table(tem_tab, TAB50EXON, quote=T, sep="\t", row.names=F)
tem_tab <- exon_rf$importance[gini_exon_rf_ind,]
tem_tab <- cbind(rownames(tem_tab), round(tem_tab,4))
write.table(tem_tab, TABIMPEXON, quote=T, sep="\t", row.names=F)


    cat("\n")
print("###### 2) run a RF with one bait per gene ######")
        # 2.1) rm non interesting baits
data_PerGene <- Keep_MaxBait_PerGene(gini_exon_rf0, new_x1)

    # 2.2) run the RF
new_x2 <- data_PerGene$mat_x
new_xtest2 <- xtest[,PrePro_findIndex(colnames(new_x2), colnames(xtest))]

print(sprintf("Growing random forest of %i trees...", NTREES))
gene_rf <- randomForest(x=new_x2, y=y, xtest=new_xtest2, ytest=ytest, ntree=NTREES, proximity = TRUE, importance=T, sampsize=SAMP_SIZE2)
print(gene_rf)

    # 2.3) extract best baits per gene
gini_gene_rf <- sort(gene_rf$importance[,10], decreasing=T)
N_cate_rfGn <- table(sapply(names(gini_gene_rf), get_elements)) # give the total number of gene baits and promoter baits to be used in the bootstrap analysis

gini_gene_rf_ind <- order(gene_rf$importance[,10], decreasing=T)
best20_GenePromot <- names(gini_gene_rf[1:20])
tab_best20_genes <- (gene_rf$importance[gini_gene_rf_ind,])[1:20,]

tem_tab <- cbind(rownames(tab_best20_genes), round(tab_best20_genes,4))
write.table(tem_tab, TAB20EXON, quote=T, sep="\t", row.names=F)
tem_tab <- gene_rf$importance[gini_gene_rf_ind,]
tem_tab <- cbind(rownames(tem_tab), round(tem_tab,4))
write.table(tem_tab, TABIMPGn, quote=T, sep="\t", row.names=F)


    cat("\n")
print("###### 3) run a RF with one bait per contig (either PMT either exon) ######")
    # 3.1) rm non interesting baits
all_targ <- sapply(names(gini_gene_rf), Pro_ExonName)
infoTargGene <- get_info_AllTargGene(INFO_TARGENE_FILE, all_targ)
bestBait_PerContig <- as.character(get_1baitPerContig(infoTargGene, gini_gene_rf))

    # 3.2) run the RF
ind <- PrePro_findIndex(bestBait_PerContig, colnames(new_x2))
new_x3 <- new_x2[,ind]
new_xtest3 <- xtest[,PrePro_findIndex(colnames(new_x3), colnames(xtest))]
print(sprintf("Growing random forest of %i trees...", NTREES))
contig_rf <- randomForest(x=new_x3, y=y, xtest=new_xtest3, ytest=ytest, ntree=NTREES, proximity = TRUE, importance=T, sampsize=SAMP_SIZE2)
    # export RF contig results
tab <- cbind(rownames(contig_rf$votes), contig_rf$votes)
write.table(tab, file=TRAIN_SET_VOTES, sep="\t", quote=F, row.names=F)
tab <- cbind(rownames(contig_rf$test$votes), contig_rf$test$votes)
write.table(tab, file=TEST_SET_VOTES, sep="\t", quote=F, row.names=F)

print(contig_rf)

    # 3.3) extract best baits per contig
gini_contig_rf <- sort(contig_rf$importance[,10], decreasing=T)
gini_contig_rf_ind <- order(contig_rf$importance[,10], decreasing=T)
best20_contigPromot <- names(gini_contig_rf[1:20])
tab_best20_contig <- (contig_rf$importance[gini_contig_rf_ind,])[1:20,]

tem_tab <- cbind(rownames(tab_best20_contig), round(tab_best20_contig,4))
write.table(tem_tab, TAB20CONTIG, quote=F, sep="\t", row.names=F)
tem_tab <- contig_rf$importance[gini_contig_rf_ind,]
tem_tab <- cbind(rownames(tem_tab), round(tem_tab,4))
write.table(tem_tab, TABIMPCONTIG, quote=F, sep="\t", row.names=F)

indiv_votes <- round(contig_rf$votes,2)
test_indiv_votes <- round(contig_rf$test$votes,2)
tvotes <- rbind(cbind(rownames(indiv_votes),indiv_votes), cbind(rownames(test_indiv_votes), test_indiv_votes))
indiv_predic <- as.character(contig_rf$predicted)
test_indiv_predic <- as.character(contig_rf$test$predicted)
tab <-  cbind(tvotes, c(indiv_predic, test_indiv_predic))
write.table(tab, ASSIGN_TAB, quote=F, sep="\t", row.names=F)


print("###### 4) check CN for best baits and their contigous baits ######")
    # 4.1) chose relevant individuals
good_indiv <- which(as.character(y)==contig_rf$predicted)
x_prim <- x[good_indiv,]
y_prim <- y[good_indiv]

    # 4.2) index for races
ind_races <- mclapply(races_uniq, grep, rownames(x_prim))
names(ind_races) <- races_uniq
    # 4.3) check the best 20 baits
bests20 <- rownames(tab_best20_contig)

system("mkdir -p Res_RF_bestGenes/")
sapply(1:length(bests20), function(x) check_baits4CN(bests20[x], x_prim=x_prim, y_prim=y_prim, contig_rf=contig_rf, mat_imp_NmaxTrees=mat_imp_NmaxTrees, info_TarGene_file=INFO_TARGENE_FILE, ind_races=ind_races, outdir="./Res_RF_bestGenes/", rang=x,races_uniq = races_uniq))

    # 4.4) the 4 best loci per race
contig_rf_imp <- contig_rf$importance[,-(9:10)]
check_baits4CN_perRace(contig_rf_imp, x_prim, y_prim=y_prim, contig_rf, mat_imp_NmaxTrees, INFO_TARGENE_FILE, ind_races, outdir="./Res_RF_bestGenes/", n_best=4,races_uniq)

argv <- commandArgs(TRUE)[1]
outFileName <- argv[1]
ver(sprintf("Saving data to %s",outFileName))
rm(argv)
save.image(file=outFileName)


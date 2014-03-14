#!/bin/Rscript
rm(list=ls())
library(ape)
library(randomForest)
library(parallel)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")
source("../utils/randomForest_helperFuns.R")

source("params.R")
load(PREVIOUS_DATA)
source("params.R")  # needed as the previous will load former params!
source("./functions.R")

set.seed(1)

cat("\n")
print("########################################################
Perform test to detect gene category with significant effect to distinguish races
########################################################")


cat("\n")
print("##### 1) observed sum of ranks")
genes <- addZeroImpGenes(gini_contig_rf)
group <- sapply(names(genes), get_elements)
categ <- sort(unique(group))
final_nber <- length(gini_contig_rf)
sum_ranks <- get_rk_sum(names(genes), n_imp=final_nber, kept=length(genes), categ)
sum_ranks <- data.frame(grp=categ, rnk=sum_ranks)
rownames(sum_ranks) <- categ


cat("\n")
print("##### 2) P being same contig")
nn <- paste(sapply(names(genes), collapse_elements, sep="_", what=1:2, colla="_"), "_", sep="")
bait_nam <- get_1rdom_bait_per_gn(nn)

data4plot_same_cont <- get_data4plot_same_contig(gini_gene_rf, INFO_TARGENE_FILE, bait_nam, n_sim=N_SIM_OLD, n_cores = N_CORES)

P_Gn_obs <- data4plot_same_cont$P_Gn_obs
P_PMT_obs <- data4plot_same_cont$P_PMT_obs
P_Gn_PMT_obs <- data4plot_same_cont$P_Gn_PMT_obs
ds_Gn <- data4plot_same_cont$ds_Gn
ds_PMT <- data4plot_same_cont$ds_PMT
ds_Gn_PMT <- data4plot_same_cont$ds_Gn_PMT
distr_rdom_P_Gn_cont_Gn <- data4plot_same_cont$P_Gn_sim
distr_rdom_P_PMT_cont_PMT <- data4plot_same_cont$P_PMT_sim
distr_rdom_P_Gn_cont_PMT <- data4plot_same_cont$P_Gn_PMT_sim


cat("\n")
print("##### 3) P being same contig new random algo")
system.time(
    sims2 <- mclapply(1:N_SIM_NEW, function (x)
        get_baits_per_pairs(x, bait_names=bait_nam, P_PMT=P_PMT_obs, P_Gn=P_Gn_obs, P_Gn_PMT=P_Gn_PMT_obs, info_TargGene_fil=INFO_TARGENE_FILE, gini_gene_rf=gini_gene_rf, ds_pmt=ds_PMT, ds_gns=ds_Gn, ds_gns_pmt=ds_Gn_PMT, inc=10, inc2=5, inc3=5, verbose=0)
    , mc.cores=N_CORES)
)

P_PMTs <- sapply(sims2, function(sol) sol$P_PMT)
P_Gns <- sapply(sims2, function(sol) sol$P_Gn)
P_Gns_PMTs <- sapply(sims2, function(sol) sol$PGn_PMT)

draw_P_same_contig(P_Gn_obs,
    distr_rdom_P_Gn_cont_Gn, P_Gns, P_PMT_obs, distr_rdom_P_PMT_cont_PMT, P_PMTs, P_Gn_PMT_obs, distr_rdom_P_Gn_cont_PMT, P_Gns_PMTs)


cat("\n")
print("##### 4) compute the expected distribution of ranks per gene category")
        # 4.4.1) random drawing (Gns and PMTs distinguished)
#~    N_bait_alpMat <- count_categ()
distr_rk_rdom_mat <- get_rdom_rk_mat(bait_nam, N_cate_rfGn, gini_gene_rf, INFO_TARGENE_FILE, gini_contig_rf, n_sim=N_SIM_OLD, n_cores = N_CORES)

        # 4.4.2) random_LD drawing (Gns and PMTs distinguished)
distr_rk_rdom_LD_mat <- get_rdom_LD_rk_mat(sims2, bait_nam, gini_gene_rf, INFO_TARGENE_FILE, gini_contig_rf, n_cores = N_CORES)

        # 4.4.3) plots
            # 4.4.3.1) length(gini_contig_rf) genes with importance and all genes used for the sum of ranks
distr_rdom_rk <- apply(distr_rk_rdom_mat, 2, get_rk_sum, n_imp=final_nber, kept=length(genes), categ)
distr_rdom_LD_rk <- apply(distr_rk_rdom_LD_mat, 2, get_rk_sum, n_imp=final_nber, kept=length(genes), categ)
draw_rk_distrib(sum_ranks, distr_rdom_rk, distr_rdom_LD_rk, final_nber, length(genes))

            # 4.4.3.1) for a subset of genes
dat_obs <- names(genes)
for (i in c(N_IN_TEST, final_nber))
{
    print(paste("Sum of ranks over ", i, " baits (", i, " baits ranked from ", i, " (best) to 1 (less good))", sep=""))
    # data obs
    sum_ranks0 <- get_rk_sum(dat_obs, n_imp=i, kept=i, categ)
    sum_ranks0 <- data.frame(grp=categ, rnk=sum_ranks0)
    rownames(sum_ranks0) <- categ

    # simuls
    distr_rdom_rk0 <- apply(distr_rk_rdom_mat, 2, get_rk_sum,
        n_imp=i, kept=i, categ)
    distr_rdom_LD_rk0 <- apply(distr_rk_rdom_LD_mat, 2, get_rk_sum,
        n_imp=i, kept=i, categ)
    # plots
    res <- draw_rk_distrib(sum_ranks0, distr_rdom_rk0, distr_rdom_LD_rk0, i, i)
    nfil <- sub(".txt", paste(i, ".txt", sep=""), SIGNIF)
    write.table(res, file=nfil, row.names=F, quote=F, sep="\t")
    print("Done")
}


cat("\n")
print("##### 5) draw the expected number of gene per categ in the top x")
for (i in c(30, 50, final_nber))
{
    print(paste("Cunt baits par category over the top ", i, " baits", sep=""))
    categ <- sum_ranks[,1]
    obs_count <- get_count_top(names(gini_contig_rf), top=i, categ)
    rdom_count <- apply(distr_rk_rdom_mat, 2, get_count_top,
        top=i, categ)
    rdom_LD_count <- apply(distr_rk_rdom_LD_mat, 2, get_count_top,
        top=i, categ)
    draw_count_distrib(top=i, categ, obs_count, rdom_count, rdom_LD_count)
}


cat("\n")
print("##### 6) save data")
argv <- commandArgs(TRUE)[1]
outFileName <- argv[1]
rm(argv)
ver(sprintf("Saving data to %s",outFileName))
save.image(file=outFileName)











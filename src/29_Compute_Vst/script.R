#!/bin/Rscript
rm(list=ls())
library(ape)
library(randomForest)
library(parallel)

cat("\n")
print(" #### 29.1) load data, sources, parameters and functions")
source("params.R")
load(PREVIOUS_DATA)
ll <- ls() ; ll <- ll[-which(ll=="x")]
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
tab_cnv <- tab_cnv[!PrePro_findMeaninglessRows(tab_cnv),]


cat("\n")
print(" #### 29.2) compute Vst per polymorphic subtarget")
races <- sapply(colnames(tab_cnv), get_elements)
    #29.2.1) raw data
vec_vst <- apply(tab_cnv, 1, get_Vst, races)
        # write a subsample in a table
gd <- which(vec_vst==1)
ftab <- tab_cnv[gd,]
write.table(cbind(rownames(ftab),ftab), col.names=T, row.names=F, sep="\t", quote=F, file=SUBSAMP1)
    #29.2.2) log2 transformed
vec_vst_log2 <- apply(tab_cnv, 1, get_Vst, races, log_2=T)

    # 29.2.3) plot histograms
l_ind_fam <- sort_name_per_categ(rownames(tab_cnv)) # 29.3.1) prepare data per families
plot_dble_hist(vst_val=vec_vst, list_fam=l_ind_fam, nam_plot=PLOT_HIST_VST)  # 29.3.2) plot Vst per family


cat("\n")
print(" #### 29.3) Compute average Vst per gene/PMT - 1 value per segment")
print("# 29.3.0) load data")
load(RAW_ALPHA)
alpha_matrix0 <- alpha_matrix
alpha_matrix <- PrePro_roundToZeroFive(alpha_matrix)
gd <- colnames(alpha_matrix)%in%colnames(tab_cnv)
alpha_matrix <- alpha_matrix[,gd]
alpha_matrix0 <- alpha_matrix0[,gd]

print("# 29.3.1) for each gene, get one value per segment")
m_alpha_seg <- get_alpha_seg(alpha_matrix)
m_alpha_seg0 <- get_alpha_seg(alpha_matrix0)

print("# 29.3.2) compute Vst per gene (average of segments per genes)")
    # rounded
gene_Vst <- compute_gene_Vst(m_alpha_seg, races)
bad <- is.na(gene_Vst)
gene_Vst <- gene_Vst[!bad]

    # raw estimates
gene_Vst0 <- compute_gene_Vst(m_alpha_seg0, races)
bad <- is.na(gene_Vst0)
gene_Vst0 <- gene_Vst0[!bad]


print("# 29.3.3) prepare data per families")
l_ind_fam2 <- sort_name_per_categ(names(gene_Vst))
l_ind_fam02 <- sort_name_per_categ(names(gene_Vst0))

print("# 29.3.4) plot Vst per family")
    print("# 29.3.4.1) rounded alpha")
plot_dble_hist(vst_val=gene_Vst, list_fam=l_ind_fam2, nam_plot=PLOT_HIST_VST2, brks=seq(-1.8, 1, by=0.1), y_lim=c(0,3))
plot_dble_hist(vst_val=gene_Vst, list_fam=l_ind_fam2, nam_plot=PLOT_HIST_VST3, brks=seq(-1.8, 1, by=0.2), y_lim=c(0,2))
plot_dble_hist(vst_val=gene_Vst, list_fam=l_ind_fam2, nam_plot=PLOT_HIST_VST4, brks=seq(-2, 1, by=0.3), y_lim=c(0,2))

    print("# 29.3.4.1) raw alpha")
plot_dble_hist(vst_val=gene_Vst0, list_fam=l_ind_fam02, nam_plot=PLOT_HIST_VST02, brks=seq(-1.8, 1, by=0.1), y_lim=c(0,3))
plot_dble_hist(vst_val=gene_Vst0, list_fam=l_ind_fam02, nam_plot=PLOT_HIST_VST03, brks=seq(-1.8, 1, by=0.2), y_lim=c(0,2))
plot_dble_hist(vst_val=gene_Vst0, list_fam=l_ind_fam02, nam_plot=PLOT_HIST_VST04, brks=seq(-2, 1, by=0.3), y_lim=c(0,2))

















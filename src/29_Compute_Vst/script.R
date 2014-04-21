#!/bin/Rscript
rm(list=ls())
library(ape)
library(randomForest)
library(parallel)

cat("\n")
print("#### 29.1) load data, sources, parameters and functions")
print("# 29.1.1) load data1 and sources")
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

print("# 29.1.2) load data2 and prepare matrices")
print("# 29.3.0) load data")
tab_targ <- read.delim(TARG, stringsAsFactors=F)
tab_cnv <- t(x)
tab_cnv <- tab_cnv[!PrePro_findMeaninglessRows(tab_cnv),]

load(RAW_ALPHA)
alpha_matrix0 <- alpha_matrix
alpha_matrix <- PrePro_roundToZeroFive(alpha_matrix)
gd <- colnames(alpha_matrix)%in%colnames(tab_cnv)
alpha_matrix <- alpha_matrix[,gd]
alpha_matrix0 <- alpha_matrix0[,gd]
bad <- PrePro_findMeaninglessRows(alpha_matrix0)
alpha_matrix0 <- alpha_matrix0[!bad,]

set.seed(0)

cat("\n")
print("#### 29.2) compute Vst per polymorphic subtarget")
races <- sapply(colnames(tab_cnv), get_elements)
print("## 29.2.1) non tranformed alpha")
    print("# 29.2.1.1) compute Vst for rounded and raw alpha")
vec_vst <- apply(tab_cnv, 1, get_Vst, races)
vec_vst0 <- apply(alpha_matrix0, 1, get_Vst, races)
        # write a subsample in a table
gd <- which(vec_vst==1)
ftab <- tab_cnv[gd,]
write.table(cbind(rownames(ftab),ftab), col.names=T, row.names=F, sep="\t", quote=F, file=SUBSAMP1)

    print("# 29.2.1.2) plot histograms")
l_ind_fam <- sort_name_per_categ(rownames(tab_cnv))
l_ind_fam0 <- sort_name_per_categ(rownames(alpha_matrix0))
plot_dble_hist(vst_val=vec_vst, list_fam=l_ind_fam, nam_plot=HIST_VST, y_lim=c(0,2.5))
plot_dble_hist(vst_val=vec_vst0, list_fam=l_ind_fam0, nam_plot=HIST_VST0, y_lim=c(0,2.5))

print("## 29.2.2) log2 transformed alpha")
    print("# 29.2.2.1) compute Vst for log2 rounded and raw alpha")
vec_vst_log2 <- apply(tab_cnv, 1, get_Vst, races, log_2=T)
vec_vst0_log2 <- apply(alpha_matrix0, 1, get_Vst, races, log_2=T)

    print("# 29.2.2.2) plot histograms log2")
plot_dble_hist(vst_val=vec_vst_log2, list_fam=l_ind_fam, nam_plot=HIST_VST_L2, y_lim=c(0,2.5))
plot_dble_hist(vst_val=vec_vst0_log2, list_fam=l_ind_fam0, nam_plot=HIST_VST0_L2, y_lim=c(0,2.5))


cat("\n")
print("#### 29.3) Compute average Vst per gene/PMT - 1 value per segment")
print("## 29.3.1) for each gene, get one value per segment")
m_alpha_seg <- get_alpha_seg(alpha_matrix)
m_alpha_seg0 <- get_alpha_seg(alpha_matrix0)

print("## 29.3.2) compute Vst per gene (average of segments per genes)")
    # rounded
gene_Vst <- compute_gene_Vst(m_alpha_seg, races)
bad <- is.na(gene_Vst)
gene_Vst <- gene_Vst[!bad]

    # raw estimates
gene_Vst0 <- compute_gene_Vst(m_alpha_seg0, races)
bad <- is.na(gene_Vst0)
gene_Vst0 <- gene_Vst0[!bad]

    # rounded log2
gene_Vst_log2 <- compute_gene_Vst(m_alpha_seg, races, log_2=T)
bad <- is.na(gene_Vst_log2)
gene_Vst_log2 <- gene_Vst_log2[!bad]

    # raw estimates log2
gene_Vst0_log2 <- compute_gene_Vst(m_alpha_seg0, races, log_2=T)
bad <- is.na(gene_Vst0_log2)
gene_Vst0_log2 <- gene_Vst0_log2[!bad]

print("## 29.3.3) prepare data per families")
l_ind_fam2 <- sort_name_per_categ(names(gene_Vst))
l_ind_fam02 <- sort_name_per_categ(names(gene_Vst0))

print("## 29.3.4) plot Vst per family")
    print("# 29.3.4.1) rounded alpha")
plot_dble_hist(vst_val=gene_Vst, list_fam=l_ind_fam2, nam_plot=HIST_VST2, brks=seq(-1.8, 1, by=0.1), y_lim=c(0,3))
plot_dble_hist(vst_val=gene_Vst, list_fam=l_ind_fam2, nam_plot=HIST_VST3, brks=seq(-1.8, 1, by=0.2), y_lim=c(0,2))
plot_dble_hist(vst_val=gene_Vst, list_fam=l_ind_fam2, nam_plot=HIST_VST4, brks=seq(-2, 1, by=0.3), y_lim=c(0,2))

    print("# 29.3.4.1) raw alpha")
plot_dble_hist(vst_val=gene_Vst0, list_fam=l_ind_fam02, nam_plot=HIST_VST02, brks=seq(-1.8, 1, by=0.1), y_lim=c(0,3))
plot_dble_hist(vst_val=gene_Vst0, list_fam=l_ind_fam02, nam_plot=HIST_VST03, brks=seq(-1.8, 1, by=0.2), y_lim=c(0,3))
plot_dble_hist(vst_val=gene_Vst0, list_fam=l_ind_fam02, nam_plot=HIST_VST04, brks=seq(-2, 1, by=0.3), y_lim=c(0,3))

    print("# 29.3.4.3) rounded alpha log2")
plot_dble_hist(vst_val=gene_Vst_log2, list_fam=l_ind_fam2, nam_plot=HIST_VST2_L2, brks=seq(-1.8, 1, by=0.1), y_lim=c(0,3))
plot_dble_hist(vst_val=gene_Vst_log2, list_fam=l_ind_fam2, nam_plot=HIST_VST3_L2, brks=seq(-1.8, 1, by=0.2), y_lim=c(0,2))
plot_dble_hist(vst_val=gene_Vst_log2, list_fam=l_ind_fam2, nam_plot=HIST_VST4_L2, brks=seq(-2, 1, by=0.3), y_lim=c(0,2))

    print("# 29.3.4.4) raw alpha log2")
plot_dble_hist(vst_val=gene_Vst0_log2, list_fam=l_ind_fam02, nam_plot=HIST_VST02_L2, brks=seq(-1.8, 1, by=0.1), y_lim=c(0,3))
plot_dble_hist(vst_val=gene_Vst0_log2, list_fam=l_ind_fam02, nam_plot=HIST_VST03_L2, brks=seq(-1.8, 1, by=0.2), y_lim=c(0,3))
plot_dble_hist(vst_val=gene_Vst0_log2, list_fam=l_ind_fam02, nam_plot=HIST_VST04_L2, brks=seq(-2, 1, by=0.3), y_lim=c(0,3))


cat("\n")
print("#### 29.4) non parametric tests")
families <- c("Control", "Gr", "Or", "P450")
vec_families.0 <- sapply(names(gene_Vst), get_elements)
vec_families0.0 <- sapply(names(gene_Vst0), get_elements)

cat("\n", "## 29.4.0) all", sep="")
print(test_dif_NoParam(vec_families.0, families, gene_Vst, toprint="Non transformed rounded values", wilcox=F))
print(test_dif_NoParam(vec_families0.0, families, gene_Vst0, toprint="Non transformed raw values", wilcox=F))
print(test_dif_NoParam(vec_families.0, families, gene_Vst_log2, toprint="Transformed rounded values", wilcox=F))
print(test_dif_NoParam(vec_families0.0, families, gene_Vst0_log2, toprint="Transformed raw values", wilcox=F))

cat("\n", "## 29.4.1) Gr", sep="")
print(test_dif_NoParam(vec_families.0, c("Control", "Gr"), gene_Vst, toprint="Non transformed rounded values"))
print(test_dif_NoParam(vec_families0.0, c("Control", "Gr"), gene_Vst0, toprint="Non transformed raw values"))
print(test_dif_NoParam(vec_families.0, c("Control", "Gr"), gene_Vst_log2, toprint="Transformed rounded values"))
print(test_dif_NoParam(vec_families0.0, c("Control", "Gr"), gene_Vst0_log2, toprint="Transformed raw values"))

cat("\n", "## 29.4.2) Or", sep="")
print(test_dif_NoParam(vec_families.0, c("Control", "Or"), gene_Vst, toprint="Non transformed rounded values"))
print(test_dif_NoParam(vec_families0.0, c("Control", "Or"), gene_Vst0, toprint="Non transformed raw values"))
print(test_dif_NoParam(vec_families.0, c("Control", "Or"), gene_Vst_log2, toprint="Transformed rounded values"))
print(test_dif_NoParam(vec_families0.0, c("Control", "Or"), gene_Vst0_log2, toprint="Transformed raw values"))

cat("\n", "## 29.4.3) P450", sep="")
print(test_dif_NoParam(vec_families.0, c("Control", "P450"), gene_Vst, toprint="Non transformed rounded values"))
print(test_dif_NoParam(vec_families0.0, c("Control", "P450"), gene_Vst0, toprint="Non transformed raw values"))
print(test_dif_NoParam(vec_families.0, c("Control", "P450"), gene_Vst_log2, toprint="Transformed rounded values"))
print(test_dif_NoParam(vec_families0.0, c("Control", "P450"), gene_Vst0_log2, toprint="Transformed raw values"))


pdf()
layout(matrix(1:4, nrow=2, ncol=2, byrow=T))
get_boxplot(gene_Vst, families)
get_boxplot(gene_Vst0, families)
get_boxplot(gene_Vst_log2, families)
get_boxplot(gene_Vst0_log2, families)
dev.off()








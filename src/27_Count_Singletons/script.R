#!/bin/Rscript
rm(list=ls())

cat("\n")
print("##### load functions")
source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")
source("./params.R")
source("./functions.R")

cat("\n")
print(" #### 1) load data and colors")
load(PREVIOUS_DATA)
set.seed(0)

    # 1.1) reads count
alpha_matrix <- PrePro_roundToZeroFive(alpha_matrix)

    # 1.2) remove bad cytisus
ind_bad <- PrePro_findIndex(BAD_CYTISUS, colnames(alpha_matrix))
ind_bad2 <- PrePro_findIndex(c(BAD_CYTISUS, RACE_UNKNOWN), colnames(alpha_matrix))
alpha_matrix <- alpha_matrix[,-ind_bad]
print(dim(alpha_matrix))

# 2) remove uniformative loci
inval_rows <- PrePro_findMeaninglessRows(alpha_matrix) # loci with no variation
alpha_matrix2 <- alpha_matrix[!inval_rows,]
print(dim(alpha_matrix2))

test_singletons <- apply(alpha_matrix2, 1, test_singleton)
N_singletons <- table(test_singletons)

ratio <- N_singletons[2]/nrow(alpha_matrix2)

cat("\n")
print("##### save data")
dummy <- NULL
argv <- commandArgs(TRUE)[1]
outFileName <- argv[1]
ver(sprintf("Saving data to %s",outFileName))
save.image(file=outFileName)


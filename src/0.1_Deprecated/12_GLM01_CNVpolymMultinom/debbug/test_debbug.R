#!/bin/Rscript
library(nnet)
library(MuMIn)
library(foreign)
library(methods)
library(multicore)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./functions.R")
source("./params.R")

# MOD_ALL1 <- Polymorphism ~ LnGeneLength + LnExonLength + Family + trimmed + LnGeneLength * Family + LnExonLength * Family + Family * trimmed

main <- function(){
    
    load(PREVIOUS_DATA)
    Pre_GLMtab <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
    , Genes_Info, CATEG_FOR_GLM, CompGenes)
    GLMtab_all2 <- set_GLMtab(Pre_GLMtab, NonTrim_only=F, covar="LnExonLength")
    MOD_ALL1 <- Polymorphism ~ LnGeneLength + LnExonLength + Family + trimmed + LnGeneLength * Family + LnExonLength * Family + Family * trimmed
    FIXED_TERMS <- c("trimmed", "Family")


#     GLMtab_all2 <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
#     MOD_ALL1 <- prog ~ ses + schtyp + write + science + female + female * schtyp
#     FIXED_TERMS <- c("ses", "write")

    M_MAX <- 5

    print(MOD_ALL1)

    print("fit the main model")
    fm1 <- multinom(formula=MOD_ALL1, data = GLMtab_all2)
    print("start all tests")
    print(test_trimmed <- dredge(fm1, fixed=FIXED_TERMS, m.max=M_MAX))

}

main()

#!/bin/Rscript
rm(list=ls())
#~library('multicore')
library(lme4)
library(MuMIn)
library(methods)
library(ggplot2)

# 0) set up parallel dredge
library(parallel)
    # Set up the cluster
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("./params.R")
source("./functions.R")
source("./Models_trimmed.R")

argv <- commandArgs(TRUE)[1]

# main <- function(argv){

    load(PREVIOUS_DATA)
	set.seed(0)

######################
cat("
##########################################
13_Binom2Steps-polym_Diverg-Related: test the effect of gene features on Pr(duplication)
##########################################\n")
groups <- names(CLUSTERS)
cat(groups)

# 13.1) Set up GLM data table
cat("\n\n# 13.1) Set up GLM data table\n")
Pre_GLMtab0 <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
, Genes_Info, CATEG_FOR_GLM, NonTrimGenes, test_blocks=T, list_groups=SUBCLUSTERS, fqcy=F)

Pre_GLMtab0 <- add_Phylog_lvl(Pre_GLMtab0, CLUSTERS)
Pre_GLMtab0 <- add_DupField(Pre_GLMtab0)
Pre_GLMtab0 <- add_GeneField(Pre_GLMtab0)

res_all_groups <- list()
v_samples <- paste("sample_size_", groups, sep=""); list_samples <- list()


for (ite in seq(length(groups)))
{
    gp <- groups[[ite]]
    cat("\n\n\n     ##################\n\tgroupe ", gp,
"\n     ##################\n", sep="")


    Pre_GLMtab <- subset(Pre_GLMtab0, Phylog_lvl==gp)

    # 13.2) Check data
    cat("\n\n     # 13.2) Check data distrisbution for group '", gp, "'\n", sep="")
    GLMtab_all1 <- set_GLMtab(Pre_GLMtab, NonTrim_only=F, covar="LnExonLength")

    print(sample_size1 <- table(GLMtab_all1$Family, GLMtab_all1$trimmed))
    print(sample_size2 <- table(GLMtab_all1$Family, GLMtab_all1$Dup))
    print(sample_size3 <- table(GLMtab_all1$Family, with(GLMtab_all1, interaction(Dup, trimmed))))
    assign(v_samples[ite], list(sample_size1=sample_size1, sample_size2=sample_size2, sample_size3=sample_size3))
    list_samples[[ite]] <- get(v_samples[ite])
    nompdf <- sub(".pdf", paste("_", gp, ".pdf", sep=""), PAIRS_ALL_EXON_LENGTH)
    Draw_pdf(pairs_glm(MODP2, data=GLMtab_all1), nompdf)

    # 13.3) Effect of variables on duplication events
        # 13.3.1) fit the models
    cat("\n\n     # 13.3) Effect of variables on duplication events for group '", gp, "'\n", sep="")
    cat("\n         # 13.3.1) Fit the maximal model for group", gp, "\n")
    fm1 <- glmer(formula=MOD_ALL1, data = GLMtab_all1, family=FAMILY)
    sum_fm1 <- summary(fm1)
    print(sum_fm1, corr=F)
    nomfil <- sub(".txt", paste("_", gp, ".txt", sep=""), GLM_DUP_MAX)
    output_glm(sum_fm1, nomfil)

        # 13.3.2) Fit all other models & model averaging
    cat("\n         # 13.3.2) Fit all other models & model averaging for group '", gp, "'\n", sep="")
    clusterExport(clust, c("GLMtab_all1", "FAMILY"))
    clusterEvalQ(clust, library(lme4))
    print(system.time(test_trimmed <- pdredge(fm1, cluster=clust, fixed=FIXED_TERMS1, m.max=M_MAX)))
    
    nomfil <- sub(".txt", paste("_", gp, ".txt", sep=""), GLM_DUP_DREDGE)
    cat(capture.output(test_trimmed), file=nomfil, sep="\n")
    
    mdl_avg <- model.avg(test_trimmed, subset = delta < DELTA1)
    sum_avg1 <- summary(mdl_avg)
    print(sum_avg1) # save in a file
    nomfil <- sub(".txt", paste("_", gp, ".txt", sep=""), GLM_DUP_AVG)
    cat(capture.output(sum_avg1, nomfil), file=nomfil, sep="\n")

        # 13.3.3) draw the results
    cat("\n         # 13.3.3) Draw predicted probability of duplication for group '", gp, "'\n", sep="")
    nompdf <- sub(".pdf", paste("_", gp, ".pdf", sep=""), DUP_PDF)
    Draw_pdf(drw_pred(test_trimmed, GLMtab_all1, DELTA1), nompdf)


    # 13.4) Effect of variables on complete duplication events
    cat("\n\n     # 13.4) Effect of variables on complete duplication events for group '", gp, "'\n", sep="")
        # 13.4.1) data
    cat("\n         # 13.4.1) Remove non polymorphic genes for group '", gp, "'\n", sep="")
    GLMtab_all2 <- GLMtab_all1[GLMtab_all1$Dup==T,]

        # 13.4.2)  Fit the maximal model
    cat("\n         # 13.4.2) Fit the maximal model for group '", gp, "'\n", sep="")
    fm2 <- glmer(formula=MOD_ALL2, data = GLMtab_all2, family=FAMILY)
    sum_fm2 <- summary(fm2)
    print(sum_fm2, corr=F)
    nomfil <- sub(".txt", paste("_", gp, ".txt", sep=""), GLM_CPDUP_MAX)
    output_glm(sum_fm2, nomfil)

        # 13.4.3) Fit al other models & model averaging
    cat("\n         # 13.4.3) Fit al other models & model averaging for group '", gp, "'\n", sep="")
    print("Test all terms")
    clusterExport(clust, c("GLMtab_all2", "FAMILY"))
    clusterEvalQ(clust, library(lme4))
    test_trimmed2 <- pdredge(fm2, cluster=clust, fixed=FIXED_TERMS2, m.max=M_MAX)
    nomfil <- sub(".txt", paste("_", gp, ".txt", sep=""), GLM_CPDUP_DREDGE)
    cat(capture.output(test_trimmed2), file=nomfil, sep="\n")

    mdl_avg2 <- model.avg(test_trimmed2, subset = delta < DELTA2[ite])
    sum_avg2 <- summary(mdl_avg2)
    print(sum_avg2)
    nomfil <- sub(".txt", paste("_", gp, ".txt", sep=""), GLM_CPDUP_AVG)
    cat(capture.output(sum_avg2, nomfil), file=nomfil, sep="\n")

        # 13.4.4) draw the results
    cat("\n         # 13.4.4) Draw predicted probability of complete duplication for group '", gp, "'\n", sep="")
    nompdf <- sub(".pdf", paste("_", gp, ".pdf", sep=""), CPDUP_PDF)
    Draw_pdf(drw_pred(test_trimmed2, GLMtab_all2, DELTA2[ite], draw_CpDup=T), nompdf)

    # 13.5) record results
    assign(gp, list(ModelDup=list(max_mdl_Dup=fm1, all_mdl_Dup=test_trimmed, mdl_avg_Dup=mdl_avg),
                ModelCpDup=list(max_mdl_CpDup=fm2, all_mdl_CpDup=test_trimmed2, mdl_avg_CpDup=mdl_avg2)))
    res_all_groups[[ite]] <- get(gp)

}
names(res_all_groups) <- groups
names(list_samples) <- groups
######################
outFileName <- argv[1]
ver(sprintf("Saving data to %s",outFileName))
#     dummy <- numeric()
save(res_all_groups, list_samples, file=outFileName)
# }

# if(DEBUG)
# 	traceback(main(argv));
# if(!DEBUG)
# 	main(argv);










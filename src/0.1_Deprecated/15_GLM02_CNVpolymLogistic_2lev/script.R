#!/bin/Rscript
library('multicore')
library(nnet)
# library(ggplot2)
# library(reshape2)
#QQ~ library(fields)
library(MuMIn)
library("epicalc")
library(methods)

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
    
	Pre_GLMtab <- set_Pre_GLMtable(ListGenes, BaitsGeneNames, alpha_matrix
	, Genes_Info, CATEG_FOR_GLM, CompGenes)

	# 1) Check data
		# 1.c) draw histograms per class of interactions
	vec_int <- with(Pre_GLMtab, interaction(trimmed, Family))
	tab_draw <- cbind(Pre_GLMtab[,c(1:2, 6, 3)], LnGeneLength=log(Pre_GLMtab[,3]), TotExonLength=Pre_GLMtab[,4], LnExonLength=log(Pre_GLMtab[,4]))
	Draw_pdf(Draw_distrib(tab_draw, Vec_int1=vec_int, Vec_int2=tab_draw$trimmed, numeric_var=4:7, nrow=4), INTERACTION_HIST)

	# 2) Investigate the effect of trimmed on CNV
		# 2.1) set up the table
	GLMtab_all <- set_GLMtab(Pre_GLMtab, NonTrim_only=F)
	print(sample_size1 <- table(GLMtab_all$Polymorphism, GLMtab_all$Family))
	print(sample_size2 <- table(GLMtab_all$Family, GLMtab_all$trimmed))
	print(sample_size3 <- table(GLMtab_all$Polymorphism, interaction(GLMtab_all$Family, GLMtab_all$trimmed)))

		# 2.2) run tests
	print("fit the main model")
	fm1 <- multinom(formula=MOD_ALL1, data = GLMtab_all)
	print("start all tests")
	test_trimmed <- dredge(fm1, fixed=FIXED_TERMS, m.max=M_MAX)
	model_avg <- model.avg(test_trimmed, subset = delta < 4)

		# 2.3) best
	best_mod <- get_best(test_trimmed)
	test_best <-  get.models(test_trimmed, best_mod)[[1]]
	summ_best <- summary(test_best)
	summ_avg <- summary(model_avg)
	Output_glm_res(summ_best, GLM_RES_BEST_TRIMMED)
	output_tab <- cbind(Model_Terms=rownames(summ_avg$coefmat), round(summ_avg$coefmat, 4))
	write.table(output_tab, file=GLM_RES_AVG_TRIMMED, sep="\t", quote=F, row.names=F)

		# 2.4) without interactions (4th best)
# 	summ_best <- summary(multinom(formula=MODELS[[8]], data=GLMtab_all))
	test <- multinom(formula=Polymorphism ~ LnGeneLength + LnExonLength + Family + trimmed, data=GLMtab_all)
	summ_best <- summary(test)
	Output_glm_res(summ_best, GLM_RES_TRIMMED_NoINTERACT)

	# 3) draw the results?
# 	pp_dpolym <- get_prediction(test_best, GLMtab_all)
# 	lpp <- melt(pp_dpolym$probaTab, id.vars = c("Family", "trimmed", "LnGeneLength", "LnExonLength"), value.name = "probability")
# 	lpp$AllInter <- interaction(lpp$Family, lpp$trimmed)
# 	ggplot(lpp, aes(x = LnGeneLength, y = probability, colour = AllInter)) + geom_line() + facet_grid(variable ~ ., scales = "free")	ggplot(lpp, aes(x = LnGeneLength, y = probability, colour = Family)) + geom_line() + facet_grid(variable ~ ., scales = "free")
# 	Draw_pdf(plotmod_predic(pp_dpolym$heat_map_Complete, GLMtab_all, "Yes"), HEAT_MAP_COMP)
# 	Draw_pdf(plotmod_predic(pp_dpolym$heat_map_Partial, GLMtab_all, "No"), HEAT_MAP_PART)
    
    ######################
	outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
#     dummy <- numeric()
	
    save(test_trimmed, summ_best, sample_size1, sample_size2, sample_size3, file=outFileName)
	print(sort(test_trimmed$AIC))
# }

# if(DEBUG)
# 	traceback(main(argv));
# if(!DEBUG)
# 	main(argv);










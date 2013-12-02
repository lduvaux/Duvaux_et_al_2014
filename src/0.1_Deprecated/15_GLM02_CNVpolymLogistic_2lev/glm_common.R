###############
set_Pre_GLMtable <- function(list_gen, Baits_Gn_Name, alpha_mat,
	Tab_Genes_Info, Categ_4GLM, CompGenes, test_blocks=F, list_groups){
	# 1) remove genes from bad gene families
	good <- unlist(mclapply(list_gen, testGoodFam, Tab_Genes_Info, Categ_4GLM))
	list_gen <- list_gen[good]

	Tab_Genes_Info <- Tab_Genes_Info[Tab_Genes_Info$GeneCateg%in%Categ_4GLM, ]

	if (test_blocks) {
		# 2) estimate CNV polymorphism
		tab_gen_polym0 <- unlist(mclapply(list_gen, pro_polymGene, Baits_Gn_Name, alpha_mat, fqcy=F, raw_scor=T, blocks=T, groups=list_groups))
		group <- rep(names(list_groups), length(list_gen))
		repgens <- rep(list_gen, each=length(list_groups))
		names(tab_gen_polym0) <- paste(repgens, group, sep="_")
		
		# 3) set up other variables: family, gene length, exons length
		gene_Res <- getAllGeneLengths(Tab_Genes_Info)
		gene_Res <- gene_Res[rownames(gene_Res)%in%list_gen,]

		ind <- PrePro_findIndex(rep(list_gen, each=length(list_groups)), rownames(gene_Res))
		gene_Res <- gene_Res[ind, ]
		
		ratiolength <- gene_Res$TotExonLength/gene_Res$GeneLength
		ftab <- cbind(Polymorphism=tab_gen_polym0, gene_Res, ratioLength=ratiolength)
		ftab$trimmed <- as.factor(ifelse(repgens%in%CompGenes, "No", "Yes"))
		ftab$Race <- as.factor(group)
		ftab <- ftab [, c(1, 2, 7, 3:6)]
		}
	else {
		# 2) estimate CNV polymorphism
		tab_gen_polym0 <- unlist(mclapply(list_gen, pro_polymGene, Baits_Gn_Name, alpha_mat, fqcy=F, raw_scor=T, blocks=F))
		names(tab_gen_polym0) <- list_gen
		
		# 3) set up other variables: family, gene length, exons length
		gene_Res <- getAllGeneLengths(Tab_Genes_Info)
		gene_Res <- gene_Res[rownames(gene_Res)%in%names(tab_gen_polym0),]

		ind <- PrePro_findIndex(names(tab_gen_polym0), rownames(gene_Res))
		gene_Res <- gene_Res[ind, ]
		
		ratiolength <- gene_Res$TotExonLength/gene_Res$GeneLength
		ftab <- cbind(Polymorphism=tab_gen_polym0, gene_Res, ratioLength=ratiolength)
		ftab$trimmed <- as.factor(ifelse(rownames(ftab)%in%CompGenes, "No", "Yes"))
		}
	return(ftab)
}

testGoodFam <- function(gene, Tab_Genes_Info, Categ4GLM){
	sstab <- subset(Tab_Genes_Info, GeneAlias==gene)
	fam <- sstab$GeneCateg[1]
	res <- fam%in%Categ4GLM
	return(res)
}

getAllGeneLengths <- function(Tab_Genes_Info){
	genes <- as.character(unique(Tab_Genes_Info$GeneAlias))
	geneRes0 <- unlist(mclapply(genes, getGeneLength, Tab_Genes_Info))

	geneRes <- data.frame(family=geneRes0[seq(1, length(geneRes0), by=3)], GeneLength=as.numeric(geneRes0[seq(2, length(geneRes0), by=3)]), TotExonLength=as.numeric(geneRes0[seq(3, length(geneRes0), by=3)]))
	rownames(geneRes) <- genes
	colnames(geneRes) <- c("Family","GeneLength","TotExonLength")
	
	return(geneRes)
}

getGeneLength <- function(gene, Tab_Genes_Info){
	Gene_Info <- subset(Tab_Genes_Info, GeneAlias==gene)
	geneLength <- computeGeneLength(Gene_Info)
	return(geneLength)
}

computeGeneLength <- function(Gene_Info){
	gene <- Gene_Info$GeneAlias[1]
	famille <- Gene_Info$GeneCateg[1]
# 	print(gene)
	if (length(unique(Gene_Info$contigV1))>1) {
		GeneSize <- NA
		ExonSize <- NA
		print(paste("Warnings: gene", gene, "on several contigs. NA is used"))}
	else {
		vec_coord <- unlist(Gene_Info[, c("startV1","stopV1")])
		rg <- range(vec_coord)
		GeneSize <- abs(rg[1]-rg[2])+1
		
		tab <- Gene_Info[, c("startV1", "stopV1")]
		ExonSize0 <- sapply(seq(nrow(tab)), function(x) abs(tab[x,1]-tab[x,2])+1)
		ExonSize <- sum(ExonSize0)}
	res <- c(famille, GeneSize, ExonSize)
	return(res)
}

###############
panel.hist <- function(x, ...)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(usr[1:2], 0, 1.5) )
	h <- hist(x, plot = FALSE)
	breaks <- h$breaks; nB <- length(breaks)
	y <- h$counts; y <- y/max(y)
	rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", meth="s", cex.cor, ...)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- abs(cor(x, y, method=meth))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex.cor * r)
}

#############
Draw_distrib <- function(tab, Vec_int1, Vec_int2, numeric_var, nrow=5, ncol=4){

	if (length(Vec_int2)==1){
		Vec_int2 <- as.factor(rep(Vec_int2, nrow(tab)))
		Vec_int1 <- as.factor(paste(Vec_int2,Vec_int1, sep="."))
	}

	# set up all interactions
	lev <- levels(Vec_int1)
	Ncell <- nrow*ncol
	layout(matrix(1:Ncell, nrow=nrow, ncol=ncol, byrow=T))
	noms <- colnames(tab)

	for (p in levels(Vec_int2)){
		for (ite in numeric_var)
		{
			vec <- tab[Vec_int2==p,ite]; n <- length(vec)
			hist(vec, main=paste(p,".AllGenes", sep=""), xlab=paste(noms[ite], " (", n, ")", sep=""), breaks=12)
		}
	}

	for (p in lev){
		for (ite in numeric_var)
		{
			vec <- tab[Vec_int1==p,ite]; n <- length(vec)
			hist(vec, main=p, xlab=paste(noms[ite], " (", n, ")", sep=""), breaks=12)
		}
	}
}

############
Output_glm_res <- function(Sum, outfile, sep="\t")
{
	z <- Sum$coefficients/Sum$standard.errors
	p <- (1 - pnorm(abs(z), 0, 1))*2
	tab <- round(rbind(Sum$coefficients, Sum$standard.errors, z, p),5)
	nlevel <- length(Sum$lev)
	if (nlevel==2)
		rownames(tab) <- paste(Sum$lev[2], c("coef", "sd", "z", "p-val"), sep="_")
	else
		rownames(tab) <- paste(rownames(tab), rep(c("coef", "sd", "z", "p-val"),each=nlevel-1), sep="_")
	tab  <- redo_table(t(tab))
	write.table(tab, file=outfile, sep=sep, row.names=F, quote=F)
}

redo_table <- function(tab)
{
	return(tab <- cbind(Model_Terms=rownames(tab), tab))
}

Fit_models <- function (all_models, data_tab)
{
	a <- list()
	vAIC <- numeric()
	for (ite in 1:length(all_models))
	{
		print(names(all_models)[ite])
		mod <- all_models[[ite]]
		test <- multinom(formula=mod, data=data_tab)
		# sum_test <- summary(test)	# deprecated: summary (test) doesn't work in the function whereas it works as a simple loop!!! Probably come from the object stocked in test$call
		a[[ite]] <- test
		vAIC[ite] <- test$AIC
	}
	names(vAIC) <- names(all_models)
	names(a) <- names(all_models)
	res <- list(test_details=a, AIC=vAIC)
	return(res)
}

############
get_best <- function(dredge_object, Delta=5)
{
	good_mod0 <- subset(dredge_object, delta < Delta)
	good <- !c(F, sapply(2:length(good_mod0$df), function(x) T%in%(good_mod0$df[x]>good_mod0$df[1:(x-1)])))
	good_mod <- subset(good_mod0, good)
	df_order <- order(good_mod$df)
	if (df_order[1]==1)
		res <- attributes(good_mod)$row.names[1]
	else {
		df_order2 <- df_order[-which(df_order==1)]
		inc <- 1
		mod1 <- get.models(test_trimmed, inc)[[1]]
		for (ite in 2:nrow(good_mod))
		{
			mod2 <- get.models(good_mod, ite)[[1]]
			lrt_test <- lrtest(mod1, mod2)
			print(paste("P-val:", round(lrt_test$p.value, 4)))
			if (lrt_test$p.value>=0.06) {
				mod1 <- mod2
				inc <- ite}
		}
		res <- attributes(good_mod)$row.names[inc]
	}
	fres <- which(attributes(test_trimmed)$row.names==res)
	return(fres)
}



















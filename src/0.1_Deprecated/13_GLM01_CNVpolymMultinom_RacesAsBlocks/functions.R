source("./10_GeneralStats_functions.R")
source("./glm_common.R")

###############
set_GLMtab <- function(tab0, NonTrim_only=F, covar)
{
	if (covar!="LnExonLength" & covar!="ratioLength")
		stop("Covariate is not given")

	if (NonTrim_only)
		tab0 <- tab0[tab0$trimmed=="Yes",]

	if (covar=="LnExonLength") {
		tab <- with(tab0, data.frame(Polymorphism, Family, Race, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), LnExonLength=log(TotExonLength)))
		rownames(tab) <- rownames(tab0)
		
		bad <- which(is.na(tab$LnGeneLength) | is.na(tab$LnGeneLength)| is.na(tab$LnExonLength)| is.na(tab$LnExonLength))
		nambad <- rownames(tab)[bad]
		tab <- tab[-bad,]}

	if (covar=="ratioLength") {
		tab <- with(tab0, data.frame(Polymorphism, Family, Race, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), ratioLength=qlogis(ratioLength)))
		rownames(tab) <- rownames(tab0)

		bad <- which(is.na(tab$LnGeneLength) | is.na(tab$LnGeneLength)| is.na(tab$ratioLength)| is.na(tab$ratioLength))
		nambad <- rownames(tab)[bad]
		tab <- tab[-bad,]}

	print(paste ("The following locus have been removed from the analysis:", nambad, sep=" "))
	return(tab)
}

###############
get_prediction <- function(test_mod, GLMtab, n=30, heat_mat=T, gene_categ = c("Control", "Gr", "Or", "P450"))
{

	# 0) nrow table
	nn <- n^2 * 4 * 2	# all combinations of possible parameters
	dpolym <- as.data.frame(matrix(data=NA, nrow=nn, ncol=4))
	colnames(dpolym) <- c("Family", "trimmed", "LnGeneLength", "LnExonLength")
	
	# 1) set up data frame
		# 1.1) gene length
	rg_LnGeneLength <- range(GLMtab$LnGeneLength)
	LnGeneLength <- round(seq(rg_LnGeneLength[1], rg_LnGeneLength[2], length.out = n),3)
	dpolym[,3] <- LnGeneLength
		# 1.2) exon length
	rg_LnExonLength <- range(GLMtab$LnExonLength)
	LnExonLength <- round(rep(seq(rg_LnExonLength[1], rg_LnExonLength[2], length.out = n), each=n),3)
	dpolym[,4] <- LnExonLength
		# 1.3) gene family
	Family <- as.factor(rep(gene_categ, each = (n^2)))
	dpolym[,1] <- Family
		# 1.4) trimmed
	trimmed <- as.factor(rep(c("Yes", "No"), each = (n^2) * 4))
	dpolym[,2] <- trimmed

	dpolym <- dpolym[order(dpolym$Family, dpolym$LnGeneLength, dpolym$trimmed, dpolym$LnExonLength),]

	# 2) finalize
	ftab <- cbind(dpolym, predict(test_mod, newdata = dpolym, type = "probs", se = TRUE))
	bad <- which(ftab$LnGeneLength<ftab$LnExonLength)
	ftab[bad, c("1_NoDup", "2_PtDup", "3_CpDup")]=NA

	# 3) set up matrices for heat maps
	if (heat_mat) {
		heat_maps1 <- set_heat_mat(ftab, gene_categ, n, "Yes")
		heat_maps2 <- set_heat_mat(ftab, gene_categ, n, "No")
		res <- list(probaTab=ftab, heat_map_Complete=heat_maps1, heat_map_Partial=heat_maps2)}
	else
		res <- list(probaTab=ftab)
	return(res)
}

set_heat_mat <- function(ftab, levs1, n, trimmed)
{
	res <- list()
	inc <- 1
	for (i in levs1)
	{
		tab <- subset(ftab, ftab[,"Family"]==i & trimmed==trimmed)
		mat1 <- matrix(tab[,"1_NoDup"], nrow=n, ncol=n)
			colnames(mat1) <- tab[seq(1, nrow(tab), by=n),"LnGeneLength"]
			rownames(mat1) <- tab[1:n,"LnExonLength"]
		mat2 <- matrix(tab[,"2_PtDup"], nrow=n, ncol=n)
			colnames(mat2) <- tab[seq(1, nrow(tab), by=n),"LnGeneLength"]
			rownames(mat2) <- tab[1:n,"LnExonLength"]
		mat3 <- matrix(tab[,"3_CpDup"], nrow=n, ncol=n)
			colnames(mat3) <- tab[seq(1, nrow(tab), by=n),"LnGeneLength"]
			rownames(mat3) <- tab[1:n,"LnExonLength"]
		leg <- "For each matrix, Rownames: LnExonLength ; Colnames: LnGeneLength"
		sublist <- list(mat1, mat2, mat3, leg); names(sublist) <- c("1_NoDup", "2_PtDup", "3_CpDup", "Legend")
		res [[i]] <- sublist
	}
	return(res)
}

plotmod_predic <- function(list_heat_map, GLMtab, trimmed, mat_layout = matrix(1:9, 3, 3, byrow=T))
{
	layout(mat_layout)
	for (ite in 1:length(list_heat_map)) {
		i <- list_heat_map[[ite]]
		print(categ <- names(list_heat_map)[ite])
		for (j in 1:3)
		{
			print(polymo <- names(i)[j])
			tab <- i[[j]]
			y_exons <- as.numeric(rownames(tab))
			x_genes <- as.numeric(colnames(tab))
			image.plot(x_genes, y_exons, t(tab), xlab="Gene length (log)", ylab="Exon length (log)", main=paste(categ, polymo, sep="_"), lwd=0)

			ssGLMtab <- subset(GLMtab, Family==categ & trimmed==trimmed)
			good <- which(ssGLMtab$Polymorphism==polymo)
			ssGLMtab <- rbind(ssGLMtab[-good,], ssGLMtab[good,])
			pchs <- as.numeric(ssGLMtab$Polymorphism)
			print(nrow(ssGLMtab))
			colo <- rep("black", length(pchs))
			wh <- which(pchs==j) ; colo[wh] <- "white"
			points(ssGLMtab$LnGeneLength, ssGLMtab$LnExonLength, pch=pchs, col=colo)
			abline(0,1)
		}
	}
}

####################
pairs_glm <- function(model, data_tab)
{
	pairs(model, data_tab, upper.panel = panel.smooth, lower.panel = panel.cor, diag.panel =panel.hist)
}









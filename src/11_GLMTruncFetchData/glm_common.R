###############
set_Pre_GLMtable <- function(list_gen, Baits_Gn_Name, alpha_mat,
	Tab_Genes_Info, Categ_4GLM, NonTrimGenes, test_blocks=F, list_groups, fqcy=F){
	# 1) remove genes from bad gene families
	good <- unlist(mclapply(list_gen, testGoodFam, Tab_Genes_Info, Categ_4GLM))
	list_gen <- list_gen[good]

	Tab_Genes_Info <- Tab_Genes_Info[Tab_Genes_Info$GeneCateg%in%Categ_4GLM, ]

	if (test_blocks) {
		# 2) estimate CNV polymorphism
        if (fqcy) {
            vec_gen_polym0 <- unlist(mclapply(list_gen, pro_polymGene, Baits_Gn_Name, alpha_mat, fqcy=T, blocks=T, groups=list_groups, CpDup_Only=F))
            vec_gen_polym1 <- unlist(mclapply(list_gen, pro_polymGene, Baits_Gn_Name, alpha_mat, fqcy=T, blocks=T, groups=list_groups, CpDup_Only=T))}
        else
            vec_gen_polym0 <- unlist(mclapply(list_gen, pro_polymGene, Baits_Gn_Name, alpha_mat, fqcy=F, blocks=T, groups=list_groups))

		group <- rep(names(list_groups), length(list_gen))
		repgens <- rep(list_gen, each=length(list_groups))
		names(vec_gen_polym0) <- paste(repgens, group, sep="_")
		
		# 3) set up other variables: family, gene length, exons length
		gene_Res <- getAllGeneLengths(Tab_Genes_Info)
		gene_Res <- gene_Res[rownames(gene_Res)%in%list_gen,]

		ind <- PrePro_findIndex(rep(list_gen, each=length(list_groups)), rownames(gene_Res))
		gene_Res <- gene_Res[ind, ]
		
		ratiolength <- gene_Res$TotExonLength/gene_Res$GeneLength
        IntronLength <- gene_Res$GeneLength-gene_Res$TotExonLength

        # 4) prepare the table
        if (fqcy) {
            ftab <- cbind(Fqcy_all=vec_gen_polym0, Fqcy_CpDup=vec_gen_polym1, gene_Res, ratioLength=ratiolength, IntronLength=IntronLength)
            ftab$trimmed <- as.factor(ifelse(repgens%in%NonTrimGenes, "No", "Yes"))
            ftab$Race <- as.factor(group)}
        else {
            ftab <- cbind(Duplication=vec_gen_polym0, gene_Res, ratioLength=ratiolength, IntronLength=IntronLength)
            ftab$trimmed <- as.factor(ifelse(repgens%in%NonTrimGenes, "No", "Yes"))
            ftab$Race <- as.factor(group)}
		}
	else {
		# 2) estimate CNV polymorphism
        if (fqcy) {
            vec_gen_polym0 <- unlist(mclapply(list_gen, pro_polymGene, Baits_Gn_Name, alpha_mat, fqcy=T, blocks=F, CpDup_Only=F))
            vec_gen_polym1 <- unlist(mclapply(list_gen, pro_polymGene, Baits_Gn_Name, alpha_mat, fqcy=T, blocks=F, CpDup_Only=T))}
        else
            vec_gen_polym0 <- unlist(mclapply(list_gen, pro_polymGene, Baits_Gn_Name, alpha_mat, fqcy=F, blocks=F))
		names(vec_gen_polym0) <- list_gen
		
		# 3) set up other variables: family, gene length, exons length
		gene_Res <- getAllGeneLengths(Tab_Genes_Info)
		gene_Res <- gene_Res[rownames(gene_Res)%in%names(vec_gen_polym0),]

		ind <- PrePro_findIndex(names(vec_gen_polym0), rownames(gene_Res))
		gene_Res <- gene_Res[ind, ]

		ratiolength <- gene_Res$TotExonLength/gene_Res$GeneLength
        IntronLength <- gene_Res$GeneLength-gene_Res$TotExonLength

        # 4) prepare the table
        if (fqcy) {
            ftab <- cbind(Fqcy_all=vec_gen_polym0, Fqcy_CpDup=vec_gen_polym1, gene_Res, ratioLength=ratiolength, IntronLength=IntronLength)
            ftab$trimmed <- as.factor(ifelse(rownames(ftab)%in%NonTrimGenes, "No", "Yes"))}
        else {
            ftab <- cbind(Duplication=vec_gen_polym0, gene_Res, ratioLength=ratiolength, IntronLength=IntronLength)
            ftab$trimmed <- as.factor(ifelse(rownames(ftab)%in%NonTrimGenes, "No", "Yes"))}
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
set_GLMtab <- function(tab0, NonTrim_only=F)
{
	if (NonTrim_only)
		tab0 <- tab0[tab0$trimmed=="Yes",]

    vzero <- which(tab0$IntronLength==0)
    tab0$IntronLength[vzero] <- tab0$IntronLength[vzero] + 1.1  # remove zero length introns from data sets

    tab <- with(tab0, data.frame(Polymorphic, Duplication, Phylog_lvl, Race, Gene, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), LnExonLength=log(TotExonLength), LnIntronLength=log(IntronLength), ratioLength=qlogis(ratioLength)))
    bad <- is.na(tab$LnGeneLength)

    rownames(tab) <- rownames(tab0)
    nambad <- rownames(tab)[bad]
    tab <- tab[-bad,]

	print(paste ("Because of Nas, the following locus have been removed from the analysis:", nambad, sep=" "))
	return(tab)
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

#############
output_glm <- function(sum_glm, nom_fil)
{
    cat(capture.output(sum_glm), file=nom_fil, sep="\n")
    tab_only <- sum_glm$coefficients
    tab_only <- cbind(rownames(tab_only), round(tab_only,5))
    nom_fil <- sub(".txt", "_coef.txt", nom_fil)
    write.table(tab_only, file=nom_fil, sep="\t", quote=F, row.names=F)
}

############### add fields to Pre_GLMtab0
add_Phylog_lvl <- function(tab, clusters)
{
    Phylog_lvl <- rep('divergent', nrow(tab))
    ind <- tab$Race%in%clusters$related
    Phylog_lvl[ind] <- "related"
    tab <- cbind(tab, Phylog_lvl)
    return(tab)
}

add_PolymField <- function (tab)
{
    Polymorphic <- tab$Duplication!="1_NoDup"
    tab <- cbind(tab, Polymorphic)
    return(tab)
}

add_GeneField <- function (tab)
{
    getGeneName <- function(stg) {
        ve <- unlist(strsplit(stg, "_"))[1:2]
        nom <- paste(ve, collapse="_")
        return(nom)}

    vec <- rownames(tab)
    Gene <- sapply(vec, getGeneName)
    tab <- cbind(tab, Gene)
    return(tab)
}











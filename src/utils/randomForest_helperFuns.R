# individuals not assigned to races
# names ; criteria to let unknown ; STRUCTURE assignment (µsat)
	# "Medicago_204_T60"; CNV-NJ ; Medicago sativa 
	# "Pisum_248_T98" ; CNV-NJ/µsat ; Trifolium pratense/Medicago sativa
	# "Pisum_193_T94" ; CNV-NJ/µsat ; Trifolium pratense/Medicago sativa
	# "Lathyrus_167_T18" ; CNV-NJ/µsat ; Trifolium pratense
	# "L.corn._243_T36" ; CNV-NJ/µsat ; Trifolium pratense
	# "Trifolium_N319_T118" ; CNV-NJ ; Trifolium pratense
	# "L.ped._N184_T54" ; CNV-NJ/µsat ; Pisum/Lotus pedunculatus
	# "Ononis_277_T85" ; CNV-NJ ; ?
	# "Pisum_246_T96" ; CNV-NJ/µsat ; Medicago sativa/Trifolium pratense
	# "Medicago_221_T69" ; CNV-NJ/µsat ; Medicago sativa/Trifolium pratense
	# "Medicago_138_T57" ; CNV-NJ ; Medicago sativa
	# "Lathyrus_282_T25" ; CNV-NJ ; ?
	# "Pisum_247_T97" ; CNV-NJ/µsat ; Trifolium pratense/Medicago sativa
	# "Pisum_192_T93" ; CNV-NJ/µsat ; Trifolium pratense
	# "Medicago_N143_T74" ; CNV-NJ ; Medicago sativa
	# "L.corn._239_T35" ; CNV-NJ/µsat ; Medicago lupulina/Lotus pedunculatus/Lotus corniculatus
	# "Trifolium_136_T110" ; CNV-NJ/µsat ; Medicago sativa/Trifolium pratense
	# "Pisum_249_T99" ; CNV-NJ/µsat ; Medicago sativa
	# "Medicago_N341_T78" ; CNV-NJ ; Medicago sativa
	# "Ononis N123 T89" ; µsat/mid distance (Ononis,weird group) on the tree ; Ononis/Trifolium

# individuals weirdly placed on the tree but not tagged as unknown
	# "Medicago 140 T58" ; good µsat assignment/mid distance between Medicago and pisum
	#  ; good µsat assignment
	# some trifolium hard to distangle

library(parallel)
Imptce_randomForest_superv <- function(N_TREES, tdata_mat, races, samp_size1)
{
	res <- randomForest(x=tdata_mat, y=as.factor(races), ntree=N_TREES, proximity = T, importance=T, sampsize=samp_size1)
	return(res)
}

Imptce_Stability_superv <- function(n, N_TREES, tdata_mat, races, samp_size1){
	multi_rf <- mclapply(rep(N_TREES, n), Imptce_randomForest_superv, tdata_mat, races, samp_size1)
	names(multi_rf) <- paste("rf", 1:n, sep="")
	return(multi_rf)
}

Imptce_randomForest_unsuperv <- function(N_TREES, x, y, xtest, ytest, samp_size2)
{
	res <- randomForest(x=x, y=y, xtest=xtest, ytest=ytest, ntree=N_TREES, proximity = T, importance=T, sampsize=samp_size2)
	return(res)
}

Imptce_Stability_unsuperv <- function(n, N_TREES, x, y, xtest, ytest, samp_size2)
{
	multi_rf <- mclapply(rep(N_TREES, n), Imptce_randomForest_unsuperv, x, y, xtest, ytest, samp_size2)
	names(multi_rf) <- paste("rf", 1:n, sep="")
	return(multi_rf)
}

merge_nRF <- function(n, lrf)
{
	samp <- sample(1:length(lrf), n)
	sp1 <- samp[1]
	comb <- lrf[[sp1]]
	
	if (n>1) {
		for (ite in 2:length(samp))
		{
			sp <- samp[ite]
			comb <- combine(comb, lrf[[sp]])
		}
	}
	return(comb)
}

set_combis_nRF <- function(lrf, lg, by=5)
{
	seqq_n <- c(seq(1, lg, by=by), lg)
	print(seqq_n)
	res <- mclapply(seqq_n, merge_nRF, lrf)
	names(res)=seqq_n
	return(res)
}

gradNtree2Imp <- function(gradNtree)
{
	Imp <- sapply(gradNtree, function(x) x$importance[,9])
	colnames(Imp)=names(gradNtree)
	return(Imp)
}

gradNtree2ImpSD <- function(gradNtree)
{
	ImpSD <- sapply(gradNtree, function(x) x$importanceSD[,9])
	colnames(ImpSD)=names(gradNtree)
	return(ImpSD)
}

sort_VarPerImp <- function(gradNtree)
{
	mat_imp <- gradNtree2Imp(gradNtree)
	mat_SD <- gradNtree2ImpSD(gradNtree)
	seqq_n <- as.numeric(colnames(mat_imp))
	sol_more_trees <- mat_imp[,ncol(mat_imp)]
	var_sorted <- sort(sol_more_trees, decreasing=T)
	return(var_sorted)
}

plot_10_bestvar <- function(gradNtree, lrf, ysup=0.01, pdf_name="Rplots.pdf")
{
    seqq_n <- as.numeric(names(gradNtree))
    mat_imp <- gradNtree2Imp(gradNtree)
    mat_SD <- gradNtree2ImpSD(gradNtree)
    sol_more_trees <- mat_imp[,ncol(mat_imp)]
	var_sorted <- order(sol_more_trees, decreasing=T)
	pdf(pdf_name)
    for (ite in 1:10)
    {
		toplot <- var_sorted[ite]
		if (ite!=1) par(new=T)
		plot(seqq_n, mat_imp[toplot,], ylim=c(0, ysup),col=ite, type="l", lwd=1.5, xlab="number of trees (.1e2)", ylab="MeanDecreaseAccuracy")
		points(seqq_n,mat_imp[toplot,] - mat_SD[toplot,]/2,col=ite, type="l", lwd=0.5, lty=3)
		points(seqq_n,mat_imp[toplot,] + mat_SD[toplot,]/2,col=ite, type="l", lwd=0.5, lty=2)
	}
	dev.off()
}

#################
Keep_MaxBait <- function(gini_NmaxTreesf, new_xf)
{
	all_names <- names(gini_NmaxTreesf)
	exon_name <- unique(sapply(all_names, Pro_ExonName))
	vec_bad <- unlist(sapply(exon_name, rm_NoMaxBait, all_names, gini_NmaxTreesf))
	
	gini_NmaxTreesf <- gini_NmaxTreesf[-vec_bad]
	new_xf <- new_xf[,-vec_bad]
	
	if (identical(colnames(new_xf), names(gini_NmaxTreesf))) {
		ll <- list(gini_NmaxTreesf=gini_NmaxTreesf, new_xf=new_xf)
		return(ll)}
	else
		stop("names are not the same in mat_xf and gini_NmaxTreesf")
}

rm_NoMaxBait <- function(DNAunit, all_names, vec_sol)
{
	ind <- grep(DNAunit, all_names)
	vec <- vec_sol[ind]
	bad <- ind[-which.max(vec)]
	return(bad)
}

#################
Keep_MaxBait_PerGene <- function(gini_small_rf, mat_x)
{
	all_names <- names(gini_small_rf)
	gene_name <- unique(sapply(all_names, Pro_geneName))
	vec_bad <- unlist(sapply(gene_name, rm_NoMaxBait, all_names, gini_small_rf))
	
	gini_small_rf <- gini_small_rf[-vec_bad]
	mat_x <- mat_x[,-vec_bad]
	
	if (identical(colnames(mat_x), names(gini_small_rf))) {
		ll <- list(gini_small_rf=gini_small_rf, mat_x=mat_x)
		return(ll)
		}
	else
		stop("names are not the same in mat_x and gini_small_rf")
}

####################
all_exon_names <- function(gini_gene_rf)
{
	exon_names_raw <- grep("PMT", names(gini_gene_rf), value=T, invert=T)
	exon_names <- sort(sapply(exon_names_raw, Pro_ExonName))
	return(exon_names)
}

all_PMT_names <- function(gini_gene_rf)
{
	bestPMT <- grep("PMT", names(gini_gene_rf), value=T)
	bestPMT <- unlist(sapply(bestPMT, function(x) strsplit(x, "|", fixed=T)))
	PMT_names <- sapply(bestPMT, Pro_ExonName)
	PMT_names <- sort(as.character(PMT_names))
	return(PMT_names)
}

get_info_AllTargGene <- function(info_TarGene_file, all_targ)
{
	infoTargGene <- PrePro_fetchInfo_TarGene(info_TarGene_file)
	ind <- PrePro_findIndex(all_targ, infoTargGene$NewTargetName)
	infoTargGene <- infoTargGene[ind,]
	if (length(all_targ)==nrow(infoTargGene))
		return(infoTargGene)
	else
		stop("problem indexing")
}

get_1baitPerContig <- function(infoTargGene, gini_gene_rf, verbose=F)
{  
    vcontigs <- unique(infoTargGene$contigV2)
	if (verbose) print(paste("Number of contigs:", length(vcontigs)))
	all_bests <- sapply(vcontigs, get_BestBaitperContig, infoTargGene, gini_gene_rf, verb=F)
	if (verbose) print(paste("Number of best baits:", length(all_bests)))
	return(all_bests)
}

get_BestBaitperContig <- function(contig, infoTargGene, gini_gene_rf, verb=F)
{
    if (verb) print(contig)
    tab <- subset(infoTargGene, contigV2==contig)
	if (nrow(tab)>1) {
		baits <- tab$NewTargetName
		ind <- sapply(baits, grep_bait, names(gini_gene_rf))
		if (!is.integer(ind) | length(ind) != length(baits)) stop("Indexing problem")
#		best <- baits[which.min(ind)]
		best <- names(gini_gene_rf)[min(ind)]
		}
	else {
#~		ind <- grep(paste("^", tab$NewTargetName, "_", sep=""), names(gini_gene_rf))
		ind <- grep_bait(tab$NewTargetName, names(gini_gene_rf))
		if (!is.integer(ind)) stop("Indexing problem")
		best <- names(gini_gene_rf)[ind]}
	return(best)
}

grep_bait <- function(pattern, x)
{
	ind <- which(pattern==x)
	if (length(ind)==1)
		return(ind)
	else {
		ind <- grep(paste("^", pattern, "[_\\|]", sep=""), x)
		if (length(ind)==0) ind <- grep(paste("\\|", pattern, "$", sep=""), x)
	}
	return(ind)
}

plot_CN_contigous_baits <- function(bait,x_prim, y_prim, races_3, races_uniq, ylim=NULL, mat_imptce, best_score, nom, test=F, locmain=NULL)
{
	v_col <- rep("black", 8); names(v_col) <- races_uniq
	v_col[races_3[1]] <- "red"
	v_col[races_3[2]] <- "orange"
	v_col[races_3[3]] <- "green"
	
	if (bait==nom) {
		score <- round(best_score[-c(9:10)],4)
		colmain <- "red"
		}
	else {
		score <- round(mat_imptce[bait,-c(9:10)],4)
		colmain <- "black"}
	gini <- round(mat_imptce[bait,10],4)
	
	if (is.null(locmain)) 
		main <- paste(bait, " (", gini, ")", sep="")
	else
		main <- paste(locmain, " (", gini, ")", sep="")
	
	if (test) {
		mp <- boxplot(x_prim[, bait]~y_prim, plot=F)
		return(max(c(mp$stats, mp$conf, mp$out)))
	}
	else {
		if (is.null(ylim))
			mp <- boxplot(x_prim[, bait] ~ y_prim, border=v_col, xlab="", ylab="Nber of copies", main=main, names = rep("",8), col.main=colmain)
		else
			mp <- boxplot(x_prim[, bait] ~ y_prim, border=v_col, xlab="", ylab="Nber of copies", main=main, names = rep("",8), ylim=ylim, col.main=colmain)
		text(1:8, par("usr")[3]- 0.03, labels = paste(races_uniq, "\n", score, sep=""), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9, col=v_col)
	}
}

get_bait_spec <- function(nom, x_prim, contig_rf, ind_races, do_print=F, print_CN_info=F)
{
	best_score <- contig_rf$importance[nom,-c(9:10)]
	if (print_CN_info) {
		print(tab <- table(x_prim[, nom]))
		print(rk <- sort(best_score, decreasing=T))}
	else {
		tab <- table(x_prim[, nom])
		rk <- sort(best_score, decreasing=T)}
	races_3 <- names(rk)[1:3]
	if (do_print){
		print(x_prim[ind_races[[races_3[1]]],nom])
		print(x_prim[ind_races[[races_3[2]]],nom])
		print(x_prim[ind_races[[races_3[3]]],nom])
	}
	ll <- list(best_score=best_score, races_3=races_3)
	return(ll)
}



check_baits4CN_perRace <- function(contig_rf_imp, x_prim, y_prim, contig_rf, mat_imp_NmaxTrees, info_TarGene_file, ind_races, outdir="./", n_best=4,races_uniq)
{
	racs <- colnames(contig_rf_imp)
	for (ite in 1:length(racs))
	{
#		cat("\n");cat("\n")
		rac <- racs[ite]
#		print(rac)
		vec_imp <- sort(contig_rf_imp[,ite], decreasing=T)
		best_n_best <- names(vec_imp)
		sapply(1:n_best, function(x) check_baits4CN(best_n_best[x], x_prim, y_prim, contig_rf, mat_imp_NmaxTrees, info_TarGene_file, ind_races=ind_races, outdir=outdir, rang=x, pref_nompdf=rac,races_uniq=races_uniq))
	}
}

check_baits4CN <- function(nom, x_prim, y_prim, contig_rf, mat_imp_NmaxTrees, info_TarGene_file, ind_races, outdir="./", pref_nompdf="GblBest", rang,races_uniq)
{
# 	print(nom)
	pref_nompdf <- paste(pref_nompdf, "_", rang, "-", sep="")
	test_PMT <- grepl("^PMT_", nom)
	if (test_PMT)
		check_baits4CN_PMT(nom, x_prim, y_prim, contig_rf, mat_imp_NmaxTrees, info_TarGene_file, ind_races=ind_races, outdir=outdir, pref_nompdf,races_uniq=races_uniq)
	else
		check_baits4CN_gene(nom, x_prim, y_prim, contig_rf, mat_imp_NmaxTrees, ind_races=ind_races, outdir=outdir, pref_nompdf,races_uniq=races_uniq) 
}

check_baits4CN_gene <- function(nom, x_prim,y_prim, contig_rf, mat_imp_NmaxTrees, ind_races, outdir="./", pref_nompdf,races_uniq){
#	print("test gene")
	# the bait itself
	a  <- 0
	ll <- get_bait_spec(nom, x_prim, contig_rf, ind_races)
	races_3 <- ll$races_3
	best_score <- ll$best_score
	# contiguous baits
	gen <- Pro_geneName(nom)
    other_baits_ind <- grep(paste(gen, "_", sep=""), colnames(x_prim))	# contiguous only present in x_prim
  
    other_baits_imp <- sort(mat_imp_NmaxTrees[other_baits_ind, "MeanDecreaseGini"], decreasing=T)
    if (length(other_baits_ind)==1)
		all_baits <- nom
	else
		all_baits <- names(other_baits_imp)

    rg <- max(sapply(all_baits, plot_CN_contigous_baits,x_prim, y_prim, races_3, races_uniq, mat_imptce=mat_imp_NmaxTrees, best_score=best_score, nom=nom, test=T))
    pdf(paste(outdir, "Res_", pref_nompdf, "CNV_DistrPerRace-", nom, ".pdf", sep=""))
	layout(matrix(1:6, 3, 2, byrow = TRUE))

	par(las=2)
	sapply(all_baits, plot_CN_contigous_baits,x_prim, y_prim, races_3, races_uniq, mat_imptce=mat_imp_NmaxTrees, best_score=best_score, nom=nom, ylim=c(0, rg))
	dev.off()
}

grep_bait2 <- function(pattern, x)
{
	ind <- which(pattern==x)
	if (length(ind)==0)
		ind <- grep(paste("^", pattern, "_", sep=""), x)
	return(ind)
}

check_baits4CN_PMT <- function(nom, x_prim,y_prim, contig_rf, mat_imp_NmaxTrees, info_TarGene_file, ind_races, outdir="./", pref_nompdf,races_uniq)
{
	# the bait itself
	ll <- get_bait_spec(nom, x_prim, contig_rf, ind_races)
	races_3 <- ll$races_3
	best_score <- ll$best_score
	
	# contiguous baits
	nom0 <- Pro_ExonName(nom)
	other_baits0 <- get_allPMTfor1gene(nom0, info_TarGene_file)
	other_baits <- other_baits0$all_alias
	other_baits_ind <- unlist(sapply(other_baits, grep_bait2, colnames(x_prim)))	# contiguous only present in x_prim, so need to use unlist
	if (is.matrix(other_baits_ind))
		other_baits_ind <- as.numeric(other_baits_ind)

    other_baits_imp <- sort(mat_imp_NmaxTrees[other_baits_ind, "MeanDecreaseGini"], decreasing=T)
    
    if (length(other_baits_ind)==1) {
		corresp_other_baits <- PrePro_findIndex(names(other_baits_ind), other_baits)
		all_baits <- nom}
    else {
		tt <- unique(sapply(names(other_baits_imp), Pro_ExonName))	# needed for PMT splitted in several baits
		corresp_other_baits <- unlist(PrePro_findIndex(tt, other_baits))
		all_baits <- names(other_baits_imp)}
    corresp_alias3 <- other_baits0$all_alias3[corresp_other_baits]
    all_baits2 <- paste(all_baits, corresp_alias3, sep="/")
    
    rg <- max(sapply(all_baits, plot_CN_contigous_baits,x_prim, y_prim, races_3, races_uniq, mat_imptce=mat_imp_NmaxTrees, best_score=best_score, nom=nom, test=T))
    pdf(paste(outdir, "Res_", pref_nompdf, "CNV_DistrPerRace-", nom, ".pdf", sep=""))
	layout(matrix(1:6, 3, 2, byrow = TRUE))
	par(las=2)
	sapply(1:length(all_baits), function(x) plot_CN_contigous_baits(bait=all_baits[x],x_prim, y_prim, races_3, races_uniq, ylim=c(0, rg), mat_imptce=mat_imp_NmaxTrees, best_score=best_score, nom=nom, locmain=all_baits2[x]))
	dev.off()
}

get_allPMTfor1gene <- function(nom0, info_TarGene_file)
{
	infoTargGene_full <- PrePro_fetchInfo_TarGene(info_TarGene_file)
	ind <- PrePro_findIndex(nom0, infoTargGene_full$NewTargetName)
	nom_alias3 <- infoTargGene_full$Aliases3[ind]
	suf <- paste(unlist(strsplit(nom_alias3, "PMT"))[1], "PMT", sep="")
	all_ind <- grep(suf, infoTargGene_full$Aliases3)
	all_alias <- infoTargGene_full$NewTargetName[all_ind]
	all_alias3 <- infoTargGene_full$Aliases3[all_ind]
	ll <- list(all_alias=all_alias, all_alias3=all_alias3)
	return(ll)
}

##########
ctrl_good_exonNber <- function(gene, infoTargGene, all_exons_CNVdata)
{
	exon_infoTargGene <- sort(grep(paste(gene, "\\.", sep=""), infoTargGene$NewTargetName, value=T))
	exon_CNVdata <- sort(grep(paste(gene, "\\.", sep=""), all_exons_CNVdata, value=T))
	
	return(identical(exon_infoTargGene, exon_CNVdata))
}

#############
check_noSNP <- function(contig_exon, start_exon, end_exon, tab_SNPs)
{
	tab <- subset(tab_SNPs, CHROM==contig_exon)
	pos_SNPs <- tab$POS
	if (nrow(tab)==0)
		 test=T
	else 
		test <- all(pos_SNPs<start_exon | pos_SNPs>end_exon)
	return(test)
}











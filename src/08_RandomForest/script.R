#!/bin/Rscript

library(ape)
library(randomForest)
library(parallel)

source("../utils/functions.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

source("../utils/randomForest_helperFuns.R")
source("./params.R")
source("./functions.R")

#~main <- function(argv){

	load(PREVIOUS_DATA)

    races <- PrePro_fetchRaces(RAW_DATA, CTRL_GUYS , BAD_GUYS)
    races[races=="Medicago_ctrl"] <- "Medicago"
	races_uniq <- PrePro_fetchUniqRaces(races)

	set.seed(1)
	pure_indiv <- names(y)
	ind_pure_indiv <- PrePro_findIndex(pure_indiv, indiv)

    # 1) run a RF by removing all noisy loci
        cat("\n")
	print("###### 1) run a RF by removing all noisy loci ######")
			# 1.1) extract locus importance from the biggest RF
    temp_tab <- sapply(gradNtree_lrf_unsuperv, function(x) x$importance[,"MeanDecreaseGini"])
    gini_NmaxTrees <- temp_tab[, ncol(temp_tab)]
    mat_imp_NmaxTrees <- gradNtree_lrf_unsuperv[[length(gradNtree_lrf_unsuperv)]]$importance

        # 1.2) remove uninformative variables (i.e. baits)
	bad <- which(gini_NmaxTrees==0)

    gini_NmaxTrees1 <- gini_NmaxTrees[-bad]
    new_x <- x[,-bad]

        # 1.3) keep only the better bait of each informative exon
    data_PerExons <- Keep_MaxBait(gini_NmaxTrees1, new_x)

        # 1.4) run the RF
	new_x1 <- data_PerExons$new_xf
	new_xtest1 <- xtest[,PrePro_findIndex(colnames(new_x1), colnames(xtest))]
	
	print(sprintf("Growing random forest of %i trees...", NTREES))
    exon_rf <- randomForest(x=new_x1, y=y, xtest=new_xtest1, ytest=ytest, ntree=NTREES, proximity = TRUE, importance=T, sampsize=SAMP_SIZE2)
    print(exon_rf)

        # 1.5) extract best baits per exon
    gini_exon_rf0 <- exon_rf$importance[,"MeanDecreaseGini"]
    N_cate_rfEx <- table(sapply(names(gini_exon_rf0), get_elements))

    gini_exon_rf <- sort(gini_exon_rf0, decreasing=T)
    gini_exon_rf_ind <- order(gini_exon_rf0, decreasing=T)
    best50_ExonPromot <- names(gini_exon_rf[1:50])
    tab_best50 <- (exon_rf$importance[gini_exon_rf_ind,])[1:50,]
    geneOfbest50 <- unique(sapply(best50_ExonPromot, get_elements, sep="\\."))
	
	tem_tab <- cbind(rownames(tab_best50), round(tab_best50,4))
	write.table(tem_tab, TAB50EXON, quote=T, sep="\t", row.names=F)
	tem_tab <- exon_rf$importance[gini_exon_rf_ind,]
	tem_tab <- cbind(rownames(tem_tab), round(tem_tab,4))
	write.table(tem_tab, TABIMPEXON, quote=T, sep="\t", row.names=F)


    # 2) run a RF with one bait per gene
        cat("\n")
	print("###### 2) run a RF with one bait per gene ######")
			# 2.1) rm non interesting baits
	data_PerGene <- Keep_MaxBait_PerGene(gini_exon_rf0, new_x1)
    
        # 2.2) run the RF
    new_x2 <- data_PerGene$mat_x
	new_xtest2 <- xtest[,PrePro_findIndex(colnames(new_x2), colnames(xtest))]
	
    print(sprintf("Growing random forest of %i trees...", NTREES))
    gene_rf <- randomForest(x=new_x2, y=y, xtest=new_xtest2, ytest=ytest, ntree=NTREES, proximity = TRUE, importance=T, sampsize=SAMP_SIZE2)
    print(gene_rf)
    
        # 2.3) extract best baits per gene
    gini_gene_rf <- sort(gene_rf$importance[,10], decreasing=T)
    N_cate_rfGn <- table(sapply(names(gini_gene_rf), get_elements)) # give the total number of gene baits and promoter baits to be used in the bootstrap analysis

    gini_gene_rf_ind <- order(gene_rf$importance[,10], decreasing=T)
    best20_GenePromot <- names(gini_gene_rf[1:20])
    tab_best20_genes <- (gene_rf$importance[gini_gene_rf_ind,])[1:20,]
	
	tem_tab <- cbind(rownames(tab_best20_genes), round(tab_best20_genes,4))
	write.table(tem_tab, TAB20EXON, quote=T, sep="\t", row.names=F)
	tem_tab <- gene_rf$importance[gini_gene_rf_ind,]
	tem_tab <- cbind(rownames(tem_tab), round(tem_tab,4))
	write.table(tem_tab, TABIMPGn, quote=T, sep="\t", row.names=F)


    # 3) run a RF with one bait per contig (either PMT either exon)
        cat("\n")
	print("###### 3) run a RF with one bait per contig (either PMT either exon) ######")
        # 3.1) rm non interesting baits
	all_targ <- sapply(names(gini_gene_rf), Pro_ExonName)
	infoTargGene <- get_info_AllTargGene(INFO_TARGENE_FILE, all_targ)
	bestBait_PerContig <- as.character(get_1baitPerContig(infoTargGene, gini_gene_rf))
	
        # 3.2) run the RF
	ind <- PrePro_findIndex(bestBait_PerContig, colnames(new_x2))
	new_x3 <- new_x2[,ind]
	new_xtest3 <- xtest[,PrePro_findIndex(colnames(new_x3), colnames(xtest))]
    print(sprintf("Growing random forest of %i trees...", NTREES))
    contig_rf <- randomForest(x=new_x3, y=y, xtest=new_xtest3, ytest=ytest, ntree=NTREES, proximity = TRUE, importance=T, sampsize=SAMP_SIZE2)
    print(contig_rf)

        # 3.3) extract best baits per contig
    gini_contig_rf <- sort(contig_rf$importance[,10], decreasing=T)
    gini_contig_rf_ind <- order(contig_rf$importance[,10], decreasing=T)
    best20_contigPromot <- names(gini_contig_rf[1:20])
    tab_best20_contig <- (contig_rf$importance[gini_contig_rf_ind,])[1:20,]
	
	tem_tab <- cbind(rownames(tab_best20_contig), round(tab_best20_contig,4))
	write.table(tem_tab, TAB20CONTIG, quote=F, sep="\t", row.names=F)
	tem_tab <- contig_rf$importance[gini_contig_rf_ind,]
	tem_tab <- cbind(rownames(tem_tab), round(tem_tab,4))
	write.table(tem_tab, TABIMPCONTIG, quote=F, sep="\t", row.names=F)
    
    indiv_votes <- round(contig_rf$votes,2)
    test_indiv_votes <- round(contig_rf$test$votes,2)
    tvotes <- rbind(cbind(rownames(indiv_votes),indiv_votes), cbind(rownames(test_indiv_votes), test_indiv_votes))
    indiv_predic <- as.character(contig_rf$predicted)
    test_indiv_predic <- as.character(contig_rf$test$predicted)
    tab <-  cbind(tvotes, c(indiv_predic, test_indiv_predic))
    write.table(tab, ASSIGN_TAB, quote=F, sep="\t", row.names=F)


    # 4) perform test to detect gene category with significant effect to distinguish races
    cat("\n")
    print("###### 4) perform test to detect gene category with significant effect to distinguish races ######")
        # 4.1) observed sum of ranks
	genes <- addZeroImpGenes(gini_contig_rf)
    group <- sapply(names(genes), get_elements)
	df <- data.frame(grp = group, rnk = rank(genes))
	sum_ranks <- aggregate(rnk ~ grp, df, sum)
    rownames(sum_ranks) <- sum_ranks[,1]
	print(sum_ranks)

        # 4.2) P being same contig
    nn <- paste(sapply(names(genes), collapse_elements, sep="_", what=1:2, colla="_"), "_", sep="")
    bait_nam <- get_1rdom_bait_per_gn(nn)

    data4plot_same_cont <- get_data4plot_same_contig(gini_gene_rf, INFO_TARGENE_FILE, bait_nam)

    P_Gn_obs <- data4plot_same_cont$P_Gn_obs
    P_PMT_obs <- data4plot_same_cont$P_PMT_obs
    P_Gn_PMT_obs <- data4plot_same_cont$P_Gn_PMT_obs
    ds_Gn <- data4plot_same_cont$ds_Gn
    ds_PMT <- data4plot_same_cont$ds_PMT
    ds_Gn_PMT <- data4plot_same_cont$ds_Gn_PMT
    distr_rdom_P_Gn_cont_Gn <- data4plot_same_cont$P_Gn_sim
    distr_rdom_P_PMT_cont_PMT <- data4plot_same_cont$P_PMT_sim
    distr_rdom_P_Gn_cont_PMT <- data4plot_same_cont$P_Gn_PMT_sim

        # 4.3) P being same contig new random algo
    system.time(
        sims2 <- mclapply(1:1000, function (x)
            get_baits_per_pairs(x, bait_names=bait_nam, P_PMT=P_PMT_obs, P_Gn=P_Gn_obs, P_Gn_PMT=P_Gn_PMT_obs, info_TargGene_fil=INFO_TARGENE_FILE, gini_gene_rf=gini_gene_rf, ds_pmt=ds_PMT, ds_gns=ds_Gn, ds_gns_pmt=ds_Gn_PMT, inc=10, inc2=5, inc3=5, verbose=0)
        , mc.cores=8)
    )

    P_PMTs <- sapply(sims2, function(sol) sol$P_PMT)
    P_Gns <- sapply(sims2, function(sol) sol$P_Gn)
    P_Gns_PMTs <- sapply(sims2, function(sol) sol$PGn_PMT)

    draw_P_same_contig(P_Gn_obs,
        distr_rdom_P_Gn_cont_Gn, P_Gns, P_PMT_obs, distr_rdom_P_PMT_cont_PMT, P_PMTs, P_Gn_PMT_obs, distr_rdom_P_Gn_cont_PMT, P_Gns_PMTs)

        # 4.4) compute the expected distribution of ranks per gene category
            # 4.4.1) random drawing (Gns and PMTs distinguished)
#~    N_bait_alpMat <- count_categ()
    distr_rk_rdom_mat <- get_rdom_rk_mat(bait_nam, N_cate_rfGn, gini_gene_rf, INFO_TARGENE_FILE, gini_contig_rf)

            # 4.4.2) random_LD drawing (Gns and PMTs distinguished)
    distr_rk_rdom_LD_mat <- get_rdom_LD_rk_mat(sims2, bait_nam, gini_gene_rf, INFO_TARGENE_FILE, gini_contig_rf)

            # 4.4.3) plots
                # 4.4.3.1) for all genes in gini_contig_rf
    distr_rdom_rk <- apply(distr_rk_rdom_mat, 2, get_rk_sum, length(gini_contig_rf))
    distr_rdom_LD_rk <- apply(distr_rk_rdom_LD_mat, 2, get_rk_sum, length(gini_contig_rf))
    draw_rk_distrib(sum_ranks, distr_rdom_rk, distr_rdom_LD_rk)

                # 4.4.3.1) for a subset of genes
    dat_obs <- rownames(df)
    for (i in c(30, 50))
    {
        sum_ranks0 <- get_rk_sum(dat_obs, i)
        rownames(sum_ranks0) <- sum_ranks0[,1]
        distr_rdom_rk0 <- apply(distr_rk_rdom_mat, 2, get_rk_sum, i)
        distr_rdom_LD_rk0 <- apply(distr_rk_rdom_LD_mat, 2, get_rk_sum, i)
        draw_rk_distrib(sum_ranks0, distr_rdom_rk0, distr_rdom_LD_rk0)
    }


    
        # 4.5) draw the expected number of gene per categ in the top x
    for (i in c(100, 50 , 30))
    {
        categ <- sum_ranks[,1]
        obs_count <- get_count_top(names(gini_contig_rf), top=i, categ)
        rdom_count <- apply(distr_rk_rdom_mat, 2, get_count_top,
            top=i, categ)
        rdom_LD_count <- apply(distr_rk_rdom_LD_mat, 2, get_count_top,
            top=i, categ)
        draw_count_distrib(top=i, categ, obs_count, rdom_count, rdom_LD_count)
    }













#~    # 5) check CN for best baits and their contigous baits
#~    print("###### 5) check CN for best baits and their contigous baits ######")
#~		# 5.1) chose relevant individuals
#~    good_indiv <- which(as.character(y)==contig_rf$predicted)
#~    x_prim <- x[good_indiv,]
#~    y_prim <- y[good_indiv]
    
#~		# 5.2) index for races
#~	ind_races <- mclapply(races_uniq, grep, rownames(x_prim))
#~	names(ind_races) <- races_uniq
#~		# 5.3) check the best 20 baits
#~	bests20 <- rownames(tab_best20_contig)

#~    system("mkdir -p Res_RF_bestGenes/")
#~    sapply(1:length(bests20), function(x) check_baits4CN(bests20[x], x_prim=x_prim, y_prim=y_prim, contig_rf=contig_rf, mat_imp_NmaxTrees=mat_imp_NmaxTrees, info_TarGene_file=INFO_TARGENE_FILE, ind_races=ind_races, outdir="./Res_RF_bestGenes/", rang=x,races_uniq = races_uniq))

#~		# 5.4) the 4 best loci per race
#~	contig_rf_imp <- contig_rf$importance[,-(9:10)]
#~	check_baits4CN_perRace(contig_rf_imp, x_prim, y_prim=y_prim, contig_rf, mat_imp_NmaxTrees, INFO_TARGENE_FILE, ind_races, outdir="./Res_RF_bestGenes/", n_best=4,races_uniq)
	

#~		# 6) chisq on most important contig
#~	load(RAW_DATA)

#~	N_gene_categ <- getGenePerCategTable()
	
#~	bait_name_info <- rownames(tem_tab[1:30,])
#~	categ_bait_info <- sapply(bait_name_info, function(x) unlist(strsplit(x, "_"))[1])
#~	N_baitInfo_categ <- table(categ_bait_info)

#~	bad <- !names(N_gene_categ)%in%names(N_baitInfo_categ)
#~	bad <- names(N_gene_categ)[bad]
#~	N_baitInfo_categ <- c(N_baitInfo_categ, rep(0, length(bad)))
#~	names(N_baitInfo_categ)[(length(N_baitInfo_categ)+1-length(bad)):length(N_baitInfo_categ)] <- bad

#~	ind <- PrePro_findIndex(names(N_gene_categ), names(N_baitInfo_categ))

#~	N_baitInfo_categ <- N_baitInfo_categ[ind]

#~	tab_chi <- rbind(All_gene=N_gene_categ, Informative_gene=N_baitInfo_categ)
#~	chi_info <- chisq.test(tab_chi)

#~	ratio <- sum(N_baitInfo_categ)/sum(N_gene_categ)
#~	vect_exp <- N_gene_categ*ratio

#~		# barplot
#~	tab_barp <- rbind(vect_exp,N_baitInfo_categ)
#~	tab_barp <- tab_barp[,c(1:3, 5:9, 4, 10:11)]

#~	pdf(PLOTDIAGGn)
#~	colos <- rep(c("black", "darkmagenta", rep("blue",6), rep("darkgreen", 3)), each=2)
#~	dens <- c(22, -1)
#~	barp <- barplot(tab_barp, main="Best genes to distinguish races", axisnames = FALSE, cex.axis=1.5, beside=T, col=colos, density=dens, cex.main=2, legend=c("# expected", "# observed"), args.legend=list(density=dens, cex=1.8, fill="darkgrey"))
#~	labells <- rep("", length(tab_barp))
#~	labells[seq(1, length(tab_barp), by=2)] <-  colnames(tab_barp)
#~	text(0.5+barp, par("usr")[3], labels = labells, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=1.4, font=2)
#~	dev.off()
	
#~	outFileName <- argv[1]
#~    ver(sprintf("Saving data to %s",outFileName))
#~    save(contig_rf,file=outFileName)
#~}


#~argv <- commandArgs(TRUE)[1]
#~if(DEBUG)
#~	traceback(main(argv));
#~if(!DEBUG)
#~	main(argv);










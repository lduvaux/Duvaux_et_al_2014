#~ source("./02_ProcessScripts/transformDataByFit.R")
#~ source("./00_Accessory_Scripts/getter_functions.R")
#~ library(parallel)

################ new functions to test sgnificativity ################
get_elements <- function(x, sep="_", what=c(1)){
    return(strsplit(x, sep)[[1]][what])
}

collapse_elements <- function(nom, sep="_", what=1:3, colla="_"){
    res <- paste(get_elements(nom, sep, what), collapse=colla)
    return(res)
}

count_categ <- function(){
    load(PREVIOUS_DATA2)
    noms <- rownames(alpha_matrix)
    print(system.time(count_categ <- table(sapply(noms, get_elements))))
    return(count_categ)
}

get_contig <- function(target, info_TargGene){
    return(info_TargGene$contigV2[info_TargGene$NewTargetName==target])
}

get_P_same_contig <- function(noms1, noms2, info_TargGene_fil, verbose=T){

    v1 <- sapply(noms1, collapse_elements, sep="_", what=1:3, colla="_") 
    v2 <- sapply(noms2, collapse_elements, sep="_", what=1:3, colla="_") 
    info_TargGene <- read.delim(info_TargGene_fil,stringsAsFactors=F)
    
    vcontigs1 <- subset(info_TargGene, NewTargetName%in%v1)[,"contigV2"]
    vcontigs2 <- subset(info_TargGene, NewTargetName%in%v2)[,"contigV2"]
    test_contig <- identical(vcontigs1, vcontigs2)
    if (verbose) {
        if (test_contig) print('Symetrical matrix')
        else print('Asymetrical matrix')}

    res0 <- outer(vcontigs1, vcontigs2, "==")
    res <- table(res0)

    if (length(res)==1)
        P <- 0
    else {
        if (test_contig) {
            res[2] <- res[2]-length(vcontigs1) # remove diagonal
            res <- res/2}    # discard one half of the matrix (since it is symetric)
        P <- as.numeric(res[2]/sum(res))
    }
    return(P)
}

get_distr_rdom_P_same_contig <- function(ite, bait_names, Gn=T, PMT=T, info_TargGene_fil) {
    PMTs_ind <- grep("^PMT_", bait_names)

    # 1) randomly chose the remaining genes before contig selection
    if (Gn){
    Gns <- bait_names[-PMTs_ind]
    noms_a <- get_rdom_Gn_rk(bait_names, N_cate_rfGn)}

    # 2) randomly chose the remaining PMTs before contig selection
    if (PMT){
    PMTs <- bait_names[PMTs_ind]
    noms_b <- get_rdom_PMT_rk(bait_names, N_cate_rfGn)}

    # 3) compute P
    if (Gn & PMT)
        {noms1 <- noms_a ; noms2 <- noms_b }
    else if (!PMT)
        {noms1 <- noms_a ; noms2 <- noms_a }
    else if (!Gn)
        {noms1 <- noms_b ; noms2 <- noms_b }

    if (ite%%100==0) print(paste(ite, "iterations done"))
    return(get_P_same_contig(noms1, noms2, info_TargGene_fil, verbose=F))
}

get_data4plot_same_contig <- function(gini_gene_rf, INFO_TARGENE_FILE, bait_nam)
{
            # 4.2.1) gene bait on same contig as another gene bait
    # best bait per genes, all uninformative genes removed
    print("##### Compute P Gn same contig other Gn (random drawing)")  
    Gns_gene_rf <- names(gini_gene_rf)[-grep("^PMT_", names(gini_gene_rf))] 
  P_Gn_cont_Gn_obs <- as.numeric(get_P_same_contig(Gns_gene_rf, Gns_gene_rf, INFO_TARGENE_FILE))
  
    print(system.time(distr_rdom_P_Gn_cont_Gn <- unlist(mclapply(1:1000, get_distr_rdom_P_same_contig, bait_nam, Gn=T, PMT=F, INFO_TARGENE_FILE, mc.cores=8))))
    ds_Gn <- sd(distr_rdom_P_Gn_cont_Gn)

            # 4.2.2) pmt bait on same contig as another pmt bait
    print("##### Compute P PMT same contig other PMT (random drawing)")
    PMTs_gene_rf <- names(gini_gene_rf)[grep("^PMT_", names(gini_gene_rf))] 
  P_PMT_cont_PMT_obs <- as.numeric(get_P_same_contig(PMTs_gene_rf, PMTs_gene_rf, INFO_TARGENE_FILE))

    print(system.time(distr_rdom_P_PMT_cont_PMT <- unlist(mclapply(1:1000, get_distr_rdom_P_same_contig, bait_nam, Gn=F, PMT=T, INFO_TARGENE_FILE, mc.cores=8))))
    ds_PMT <- sd(distr_rdom_P_PMT_cont_PMT)

        # 4.2.3) gene bait on same contig as a pmt bait
    print("##### Compute P Gn same contig other PMT (random drawing)") 
  P_Gn_cont_PMT_obs <- as.numeric(get_P_same_contig(Gns_gene_rf, PMTs_gene_rf, INFO_TARGENE_FILE))

    print(system.time(distr_rdom_P_Gn_cont_PMT <- unlist(mclapply(1:1000, get_distr_rdom_P_same_contig, bait_nam, Gn=T, PMT=T, INFO_TARGENE_FILE, mc.cores=8))))
    ds_Gn_PMT <- sd(distr_rdom_P_Gn_cont_PMT)

    # results
    res <- list(
            P_Gn_obs=P_Gn_cont_Gn_obs, P_Gn_sim=distr_rdom_P_Gn_cont_Gn, ds_Gn=ds_Gn,
            P_PMT_obs=P_PMT_cont_PMT_obs, P_PMT_sim=distr_rdom_P_PMT_cont_PMT, ds_PMT=ds_PMT,
            P_Gn_PMT_obs=P_Gn_cont_PMT_obs, P_Gn_PMT_sim=distr_rdom_P_Gn_cont_PMT, ds_Gn_PMT=ds_Gn_PMT)

    return(res)
}


draw_P_same_contig <- function(P_Gn_cont_Gn_obs, distr_rdom_P_Gn_cont_Gn, P_Gns, P_PMT_cont_PMT_obs, distr_rdom_P_PMT_cont_PMT, P_PMTs, P_Gn_cont_PMT_obs, distr_rdom_P_Gn_cont_PMT, P_Gns_PMTs)
{
    pdf("distrib_Proba_sameContig_3.pdf")
    layout(matrix(1:4, 2,2, byrow=T))
    rg <- range(c(P_Gn_cont_Gn_obs, distr_rdom_P_Gn_cont_Gn, P_Gns))
    pval_Gn_cont_Gn_obs <- get_pval(P_Gn_cont_Gn_obs,distr_rdom_P_Gn_cont_Gn, two_sided=T)
        plot(density(distr_rdom_P_Gn_cont_Gn), xlim=c(rg[1], rg[2]), col="blue", ylim=c(0,700), xlab="P-same contig", main="")
    pval_Gns <- get_pval(P_Gn_cont_Gn_obs, P_Gns, two_sided=T)
        par(new=T)
        plot(density(P_Gns), xlim=c(rg[1], rg[2]), col="red", ylim=c(0,700), xlab="", main=paste("Gns: Prdom=", pval_Gn_cont_Gn_obs, " ; P-new=", pval_Gns, sep=""))
        abline(v=P_Gn_cont_Gn_obs, col='green')
    
    rg <- range(c(P_PMT_cont_PMT_obs, distr_rdom_P_PMT_cont_PMT, P_PMTs))
    pval_PMT_cont_PMT_obs <- get_pval(P_PMT_cont_PMT_obs, distr_rdom_P_PMT_cont_PMT, two_sided=T)
        plot(density(distr_rdom_P_PMT_cont_PMT), xlim=c(rg[1], rg[2]), col="blue", ylim=c(0,300), xlab="P-same contig", main="")
    pval_PMTs <- get_pval(P_PMT_cont_PMT_obs, P_PMTs, two_sided=T)
        par(new=T)
        plot(density(P_PMTs), xlim=c(rg[1], rg[2]), col="red", ylim=c(0,300), xlab="", main=paste("PMTs: Prdom=", pval_PMT_cont_PMT_obs, " ; P-new=", pval_PMTs, sep=""))
        abline(v=P_PMT_cont_PMT_obs, col='green')

    rg <- range(c(P_Gn_cont_PMT_obs, distr_rdom_P_Gn_cont_PMT, P_Gns_PMTs))
    pval_Gn_cont_PMT_obs <- get_pval(P_Gn_cont_PMT_obs, distr_rdom_P_Gn_cont_PMT, two_sided=T)
        plot(density(distr_rdom_P_Gn_cont_PMT), xlim=c(rg[1], rg[2]), col="blue", ylim=c(0,600), xlab="P-same contig", main="")
    pval_Gns_PMTs <- get_pval(P_Gn_cont_PMT_obs, P_Gns_PMTs, two_sided=T)
        par(new=T)
        plot(density(P_Gns_PMTs), xlim=c(rg[1], rg[2]), col="red", ylim=c(0,600), xlab="", main=paste("Gns_PMTs: Prdom=", pval_Gn_cont_PMT_obs, " ; P-new=", pval_Gns_PMTs, sep=""))
        abline(v=P_Gn_cont_PMT_obs, col='green')
    dev.off()
}

get_rdom_PMT_rk <- function(bait_names, N_cate_rfGn)
# bait_names: set of bait representing all genes (1 bait per gene)
# N_cate_rfGn: count of bait per category at step 2.3
{
    ind <- grep("^PMT_", bait_names)
    pmts <- bait_names[ind]
    N_inf_pmts <- as.numeric(N_cate_rfGn["PMT"])
    res <- sample(pmts, N_inf_pmts)
    return(res)
}

get_rdom_Gn_rk <- function(bait_names, N_cate_rfGn)
# bait_names: set of bait representing all genes (1 bait per gene)
# N_cate_rfGn: count of bait per category at step 2.3
{
    ind <- grep("^PMT_", bait_names)
    gns <- bait_names[-ind]
    N_inf_gns <- sum(N_cate_rfGn) - as.numeric(N_cate_rfGn["PMT"])
    res <- sample(gns, N_inf_gns)
    return(res)
}

get_baits_per_pairs <- function(indic, bait_names, P_PMT, P_Gn, P_Gn_PMT, info_TargGene_fil, gini_gene_rf, ds_pmt, ds_gns, ds_gns_pmt, verbose=2, inc=10, inc2=5, inc3=5, coef_sd1=2, coef_sd2=-0.5, coef_sd3=-1)
#draw pairs of baits randomly but by respecting the probabilities of being on the same contig
# verbose can be 0 (no comments), 1(only comments outside the loops), or 2 (all comments), note that verbose==F and verbose==T results in verbose==0 and verbose==1 respectively
{
#~    system.time({
    if (!verbose%in%0:2) stop("bad value of verbose")
    # I) set up the number of PMTs and Gns
    info_TargGene <- read.delim(info_TargGene_fil, stringsAsFactors=F)
        # informative
    inf_PMT <- names(gini_gene_rf)[grep("^PMT_", names(gini_gene_rf))]
    inf_PMT_exon <- sapply(inf_PMT, collapse_elements, sep="_", what=1:3, colla="_")
    Ninf_PMT <- length(grep("^PMT_", names(gini_gene_rf)))
    PMTs <- character(Ninf_PMT)
    Ninf_Gns <- length(gini_gene_rf) - Ninf_PMT
    Gns <- character(Ninf_Gns)

        # total
    all_PMT <- bait_names[grep("^PMT_", bait_names)]
    all_PMT_exon <- sapply(all_PMT, collapse_elements, sep="_", what=1:3, colla="_")
    N_all_PMT <- length(all_PMT)
    all_Gns <- bait_names[-grep("^PMT_", bait_names)]
    all_Gn_exon <- sapply(all_Gns, collapse_elements, sep="_", what=1:3, colla="_")
    N_all_Gns <- length(all_Gns)


    # II) chose the PMTs
        # II.1) draw some initial random PMTs
    n_init <- round(1/5*Ninf_PMT)
    ite <- n_init
    PMTs[1:n_init] <- all_PMT[sample(1:N_all_PMT, n_init)]
    

        # II.2) check % same contigs
#~    system.time(
    while (ite < Ninf_PMT) {
        PMT_exons <- sapply(PMTs[1:ite], collapse_elements, sep="_", what=1:3, colla="_")
        test <- get_P_same_contig(PMTs[1:ite], PMTs[1:ite], info_TargGene_fil, verbose=F)
        if (verbose==2) print(round(test, 4))

        if (test < (P_PMT + ds_pmt * coef_sd1)) {
            if (verbose==2) print(paste("# test=", round(test,4), ", so draw baits from common contigs"))
            # pick 5 baits from the same contigs the baits already selected
            contig_test <- unique(subset(info_TargGene, NewTargetName%in%PMT_exons)[,"contigV2"])
            tab <- info_TargGene[info_TargGene$contigV2%in%contig_test,]
            tab <- subset(tab, NewTargetName%in%all_PMT_exon & !NewTargetName%in%PMT_exons)

            if (nrow(tab)!=0){
                new_exons <- tab$NewTargetName[sample(1:nrow(tab), inc)]
#~                new_baits <- all_PMT[sapply(new_exons, function(x) grep(paste("^", x, "_", sep=""), all_PMT))]
                new_baits <- all_PMT[pmatch(new_exons, all_PMT)]    # take care that only one match is possible here
                if (NA%in%new_baits)
                    stop("Several matches of new_exons in all_PMT")
                
            } else {
                if (verbose==2) print("#### no other bait available on common contigs, so draw baits at random")
                veve <- all_PMT[!all_PMT%in%PMTs[1:ite]]
                new_baits <- veve[sample(1:length(veve), inc)]
            }
        } else {
            if (verbose==2) print(paste("# test=", round(test,4), ", so draw baits at random",sep=""))
            veve <- all_PMT[!all_PMT%in%PMTs[1:ite]]
            new_baits <- veve[sample(1:length(veve), inc)]}
        
        if ((ite+inc) > Ninf_PMT){
            ins <- Ninf_PMT-ite
            PMTs[(ite+1):length(PMTs)] <- new_baits[1:ins]
        } else {
        if (length((ite+1):(ite+inc))!= length(new_baits)) stop ("problem increment inc")
        PMTs[(ite+1):(ite+inc)] <- new_baits}
        ite <- ite + inc
    }
#~    )

    # III) chose the Gns
        # III.1) draw some initial random Gns
    n_init <- round(1/5*Ninf_Gns)
    ite <- n_init
    Gns[1:n_init] <- all_Gns[sample(1:N_all_Gns, n_init)]

        # III.2) check % same contigs
    while (ite < Ninf_Gns) {
            # III.2.1) test towards other genes
        Gn_exons <- sapply(Gns[1:ite], collapse_elements, sep="_", what=1:3, colla="_")
        test2 <- get_P_same_contig(Gns[1:ite], Gns[1:ite], info_TargGene_fil, verbose=F)
        if (verbose==2) print(round(test2, 4))

        if (test2 < (P_Gn + ds_gns * coef_sd2)){
            if (verbose==2) print(paste("# test2=", round(test2,4), ", so draw baits from common contigs"))
            # pick 5 baits from the same contigs the baits already selected
            contig_test <- unique(subset(info_TargGene, NewTargetName%in%Gn_exons)[,"contigV2"])
            tab <- info_TargGene[info_TargGene$contigV2%in%contig_test,]
            tab <- subset(tab, NewTargetName%in%all_Gn_exon & !NewTargetName%in%Gn_exons)

            if (nrow(tab)!=0){
                new_exons <- tab$NewTargetName[sample(1:nrow(tab), inc2)]
#~                new_baits <- all_Gns[sapply(new_exons, function(x) grep(paste("^", x, "_", sep=""), all_Gns))]
                new_baits <- all_Gns[pmatch(new_exons, all_Gns)]
                if (NA%in%new_baits)
                    stop("Several matches of new_exons in all_PMT")
            } else {
                if (verbose==2) print("#### no other bait available on common contigs, so draw baits at random")
                veve <- all_Gns[!all_Gns%in%Gns[1:ite]]
                new_baits <- veve[sample(1:length(veve), inc2)]
            }
        } else {
            if (verbose==2) print(paste("# test2=", round(test2,4), ", so draw baits at random",sep=""))
            veve <- all_Gns[!all_Gns%in%Gns[1:ite]]
            new_baits <- veve[sample(1:length(veve), inc2)]}

            # III.2.2) update Gns
        if ((ite+inc2) >= Ninf_Gns){
            ins <- Ninf_Gns-ite
            Gns[(ite+1):length(Gns)] <- new_baits[1:ins]
        } else {
        if (length((ite+1):(ite+inc2))!= length(new_baits)) stop ("problem increment inc2")
        Gns[(ite+1):(ite+inc2)] <- new_baits}
        ite <- ite + inc2

        if (ite < Ninf_Gns){
            # III.2.3) test towards PMTs
            Gn_exons <- sapply(Gns[1:ite], collapse_elements, sep="_", what=1:3, colla="_")
            test3 <- get_P_same_contig(Gns[1:ite], PMTs, info_TargGene_fil, verbose=F)
#~            print(ite)
#~            print(Gns[1:ite])
#~            print(PMTs)
#~            print(test3)
            if (verbose==2) print(round(test3, 4))  

            if (test3 < (P_Gn_PMT + ds_gns_pmt * coef_sd3)){
                if (verbose==2) print(paste("# test3=", round(test3,4), ", so draw baits from common contigs"))
                # pick 5 baits from the same contigs the baits already selected
                contig_test <- unique(subset(info_TargGene, NewTargetName%in%PMT_exons)[,"contigV2"])
                tab <- info_TargGene[info_TargGene$contigV2%in%contig_test,]
                tab <- subset(tab, NewTargetName%in%all_Gn_exon & !NewTargetName%in%Gn_exons)

                if (nrow(tab)!=0){
                    new_exons <- tab$NewTargetName[sample(1:nrow(tab), inc3)]
    #~                new_baits <- all_Gns[sapply(new_exons, function(x) grep(paste("^", x, "_", sep=""), all_Gns))]
                    new_baits <- all_Gns[pmatch(new_exons, all_Gns)]
                    if (NA%in%new_baits)
                        stop("Several matches of new_exons in all_PMT")
                } else {
                    if (verbose==2) print("#### no other bait available on common contigs, so draw baits at random")
                    veve <- all_Gns[!all_Gns%in%Gns[1:ite]]
                    new_baits <- c(new_baits, veve[sample(1:length(veve), inc3)])
                }
            } else {
                if (verbose==2) print(paste("# test3=", round(test3,4), ", so draw baits at random",sep=""))
                veve <- all_Gns[!all_Gns%in%Gns[1:ite]]
                new_baits <- veve[sample(1:length(veve), inc3)]}

                # III.2.4) update Gns
            if ((ite+inc3) > Ninf_Gns){
                ins <- Ninf_Gns-ite
                Gns[(ite+1):length(Gns)] <- new_baits[1:ins]
            } else {
            if (length((ite+1):(ite+inc3))!= length(new_baits)) stop ("problem increment inc3")
            Gns[(ite+1):(ite+inc3)] <- new_baits}
            ite <- ite + inc3
        }
    }

    if (verbose%in%1:2) print(P_PMT)
    test <- get_P_same_contig(PMTs, PMTs, info_TargGene_fil, verbose=F)
    if (verbose%in%1:2) print(test)
    if (verbose%in%1:2) print(P_Gn)
    test2 <- get_P_same_contig(Gns, Gns, info_TargGene_fil, verbose=F)
    if (verbose%in%1:2) print(test2)
    if (verbose%in%1:2) print(P_Gn_PMT)
    test3 <- get_P_same_contig(Gns, PMTs, info_TargGene_fil, verbose=F)
    if (verbose%in%1:2) print(test3)
    if (verbose%in%1:2) print(P_Gn_PMT)
#~    })
    ll <- list(PMTs=PMTs, Gns=Gns, P_PMT=test, P_Gn=test2, PGn_PMT=test3)
    if (indic%%50==0) print(paste(indic, " iterations done")) 
    return(ll)
}

get_null_draw <- function(ite, Gns, PMTs, bait_names, n_info, info_TargGene_fil, kept, verbose=F)
# Gns: a list of genes drawn at random or at random accounting for LD
# PMTs: a list of promoters drawn at random or at random accounting for LD
# bait_names: vector bait names representing the set of possible target (only one bait per target). A target can be a Gn or a PMT.
# n_info: nber of informative baits after step 2.3
# kept: number of gene to keep in the final ranking correspond to the number of genes shown in the final analysis)
{
    # 1) set up the importance 'imp' vector
    ind <- match(c(Gns, PMTs), bait_names)
    imp <- rep(0, length(bait_names))
    names(imp) <- bait_names
    imp [ind] <- sample(seq(n_info), n_info)   # should stay ~240 genes with some importance
    imp <- sort(imp, decreasing=T)  # has to be sorted for the function 'get_1baitPerContig'

    # 2) select best bait per contig
    imp1 <- imp[imp!=0]
    if (verbose) cat(imp1, sep="\n", file="importance.txt")
    if (verbose) cat(names(imp1), sep="\n", file="genes.txt")
        
    all_targ <- sapply(names(imp1), collapse_elements, sep="_", what=1:3, colla="_")
    info_TargGene <- get_info_AllTargGene(info_TargGene_fil, all_targ)

    bestBait_PerContig <- as.character(get_1baitPerContig(info_TargGene, imp1, verbose=F))
    if (verbose) print(bestBait_PerContig)
    ind <- PrePro_findIndex(bestBait_PerContig, names(imp))

    # 3) keep only a rankable importance for the best bait best bait per contig
        # attribute an EQUAL low importance of 0 for all other baits (step 3.1 in the main script)
    imp[-ind] <- 0
        # keep only the 'kept' first loci
    bad <- order(imp, decreasing=T)[(kept+1):length(bait_names)]
    imp[bad] <- 0
    if (ite%%100==0) print(paste(ite, "iterations done"))
    return(names(sort(imp, decreasing=T)))
}

get_rk_sum <- function(rked_gns, kept)
{
    g <- sapply(rked_gns, get_elements)
    imp <- rep(0, length(rked_gns))
    imp[1:kept] <- kept:1
    df_fc <- data.frame(grp = g, rnk = rank(imp))
    res <- aggregate(rnk ~ grp, df_fc, sum)
    v <- res$rnk
    names(v) <- res$grp
    return(v)
}

get_rdom_rk_mat <- function(bait_nam, N_cate_rfGn, gini_gene_rf, INFO_TARGENE_FILE, gini_contig_rf)
{
    print(system.time(distr_rdom <- mclapply(
        1:1000, get_null_draw,
            Gns=get_rdom_Gn_rk(bait_nam, N_cate_rfGn),
            PMTs=get_rdom_PMT_rk(bait_nam, N_cate_rfGn),
            bait_names=bait_nam, n_info=length(gini_gene_rf),
            info_TargGene_fil=INFO_TARGENE_FILE,
            kept=length(gini_contig_rf), verbose=F, 
        mc.cores=8)))
    distr_rdom_mat <- sapply(seq(length(distr_rdom)), function(x) distr_rdom[[x]])
    rownames(distr_rdom_mat) <- 1:length(bait_nam) ; colnames(distr_rdom_mat) <- paste("Simul_", 1:1000, sep="")
    return(distr_rdom_mat)
}

get_rdom_LD_rk_mat <- function(sims2, bait_nam, gini_gene_rf, INFO_TARGENE_FILE, gini_contig_rf){
    print(system.time(distr_rdom_LD <- mclapply(
        seq(length(sims2)), function(x) get_null_draw(
            x, 
            Gns=sims2[[x]]$Gns, PMTs=sims2[[x]]$PMTs,
            bait_names=bait_nam, n_info=length(gini_gene_rf),
            info_TargGene_fil=INFO_TARGENE_FILE,
            kept=length(gini_contig_rf), verbose=F), 
        mc.cores=8)))
    distr_rdom_LD_mat <- sapply(seq(length(distr_rdom_LD)), function(x) distr_rdom_LD[[x]])
    rownames(distr_rdom_LD_mat) <- 1:length(bait_nam) ; colnames(distr_rdom_LD_mat) <- paste("Simul_", 1:1000, sep="")
    return(distr_rdom_LD_mat)
}

get_count_top <- function(rked_gns, top=50, categ)
{
    res <- rep(0, length(categ))
    names(res) <- categ
    rked_gns <- rked_gns[1:top]
    g <- sapply(rked_gns, get_elements)
    tab <- table(g)
    ind <- match(names(tab), names(res))
    res[ind] <- tab
    return(res)
}

addZeroImpGenes <- function(xx){

	import_genes <- sapply(names(xx),function(x){
		st <- strsplit(x,"_")[[1]][1:2]
		return(paste(st[1],st[2],sep="_"))
		})

	all_genes <- getAllGeneNames()
	zero_imp <- all_genes[!(all_genes %in% import_genes)]
	v <- rep(0,length(zero_imp))
	names(v) <- zero_imp
	out <- c(xx,v)
	return(out)
}

getAllGeneNames <- function(){
    load(PREVIOUS_DATA2)
    data_mat <- alpha_matrix
	out <- unique(sapply(rownames(data_mat),function(x){
		st <- strsplit(x,"_")[[1]][1:2]
		return(paste(st[1],st[2],sep="_"))
		}))
	return(out)
}

get_pval <- function(observ, dist_exp, two_sided=F)
# if two_sided=F, P is always P(observ>expected)
{
    pval <- round(1-length(subset(dist_exp, observ>dist_exp))/length(dist_exp), 3)
    if (two_sided) {
        pval2 <- abs(1-pval)
        pval <- ifelse(pval<=pval2, pval, pval2)*2}
    return(pval)
}

get_1rdom_bait_per_gn <- function (gene_names){
    
    load(PREVIOUS_DATA2)
    noms <- rownames(alpha_matrix)
    res <- sapply(gene_names, grep_first, noms, value=T)
    return(res)
}

grep_first <- function(patt, vec, value=T){
    return(grep(paste("^", patt, sep=""), vec, value=value)[1])
}

draw_rk_distrib <- function(sum_ranks, distr_rdom_rk, distr_rdom_LD_rk)
{
    pdf("Distrib-rk_rdom.pdf")
	for(i in sum_ranks$grp){
        obs <- sum_ranks$rnk[which(sum_ranks$grp==i)]
        rg <- range(c(obs, distr_rdom_rk[i,]))
        p_val <- get_pval(obs, distr_rdom_rk[i,], two_sided=TWOSIDED)
        tit <- paste(i, ifelse(TWOSIDED, " (two sided", " (one sided"), " P= ", p_val, ")", sep="")
		hist(distr_rdom_rk[i,], main=tit, breaks=50, xlab="Sum of the ranks", cex.main=0.9, xlim=c(rg[1], rg[2]))
		srtd <- sort(distr_rdom_rk[i,])
		abline(v = obs,col="red",lwd=3)
	}
	dev.off()
    
	pdf("Distrib-rk_rdom-LD.pdf")
	for(i in sum_ranks$grp){
        obs <- sum_ranks$rnk[which(sum_ranks$grp==i)]
        rg <- range(c(obs, distr_rdom_LD_rk[i,]))
        p_val <- get_pval(obs, distr_rdom_LD_rk[i,], two_sided=TWOSIDED)
        tit <- paste(i, ifelse(TWOSIDED, " (two sided", " (one sided"), " P= ", p_val, ")", sep="")
		hist(distr_rdom_LD_rk[i,], main=tit, breaks=50, xlab="Sum of the ranks", cex.main=0.9, xlim=c(rg[1], rg[2]))
		srtd <- sort(distr_rdom_LD_rk[i,])
		abline(v = obs,col="red",lwd=3)
	}
	dev.off()

}

draw_count_distrib <- function(top, categ, obs_count, rdom_count, rdom_LD_count)
{
    
	pdf(paste("Distrib-count", top, "_rdom.pdf", sep=""))
	for(i in 1:length(obs_count)){
        obs <- obs_count[i]
        rg <- range(c(obs, rdom_count[i,]))
        p_val <- get_pval(obs, rdom_count[i,], two_sided=TWOSIDED)
        tit <- paste(categ[i], ifelse(TWOSIDED, " (two sided", " (one sided"),
            " P= ", p_val, ")", sep="")
        tab <- table(rdom_count[i,])
		aa <- plot(tab, main=tit, xlab="Counts", cex.main=0.9, xlim=c(rg[1], rg[2]), ylab="Frequency")
        test <- as.numeric(names(tab))==obs
        y <- ifelse(T%in%test, tab[test], 0)
		points(obs, y, col="red",lwd=5)
	}
	dev.off()
    
	pdf(paste("Distrib-count", top, "_rdom-LD.pdf", sep=""))
	for(i in 1:length(obs_count)){
        obs <- obs_count[i]
        rg <- range(c(obs, rdom_LD_count[i,]))
        p_val <- get_pval(obs, rdom_LD_count[i,], two_sided=TWOSIDED)
        tit <- paste(categ[i], ifelse(TWOSIDED, " (two sided", " (one sided"),
            " P= ", p_val, ")", sep="")
        tab <- table(rdom_LD_count[i,])
		aa <- plot(tab, main=tit, xlab="Counts", cex.main=0.9, xlim=c(rg[1], rg[2]), ylab="Frequency")
        test <- as.numeric(names(tab))==obs
        y <- ifelse(T%in%test, tab[test], 0)
		points(obs, y, col="red",lwd=5)
	}
	dev.off()

}

############################################################# 

getGenePerCategTable <- function(){
	load(PREVIOUS_DATA2)
	N_gene_categ <- as.vector(gene_per_categ[,2])
	names(N_gene_categ) <- as.vector(as.character(gene_per_categ[,1]))
	return(N_gene_categ)
}


PreProNJ_GetGeneCateg <-  function(All_bait_categ){
	cate <- unique(All_bait_categ)
	cate <- c("All", cate)
	return(cate)
}


define_alpha_matrix <- function(alpha_matrix, gene_categ, i_categ, baits)
{
	categ <- gene_categ[i_categ]
	if (i_categ==1)
		talpha <- alpha_matrix
	else {
		ind <- grep(paste("^", categ, "_", sep=""), baits)
		talpha <- alpha_matrix[ind,]
	}
	return(talpha)
}

compute_ludo_dtce <- function(talpha)
{
	nind <- ncol(talpha)
	nb_baits=nrow(talpha)
	
	vdis <- numeric(sum(1:(nind-1)))
	forbid <- NULL
	inc <- 1
	for (ite in 1:(nind-1))
	{
		forbid <- c(forbid, ite)
		vec <- talpha[,ite]
		tab_temp <- talpha[,-forbid]
		if (is.matrix(tab_temp)) 
			dis <- sapply(1:ncol(tab_temp), function(x) nb_baits-sum(vec==tab_temp[,x])) 
		else 
			dis <- nb_baits-sum(vec==tab_temp)
		vdis[inc:(inc+nind-length(forbid)-1)] <- dis
		inc <- inc+nind-length(forbid)
	}
	return(vdis)
}

put_vector_into_matrix <- function(nind, vdis, indiv_names)
{
	M <- matrix(0, nind, nind)
	M[lower.tri(M)] <- vdis
	M <- t(M)
	M[lower.tri(M)] <- vdis
	dimnames(M) <- list(indiv_names)
	return(M)
}

compute_dstce <- function(talpha, method)
{
	if (method=="ludo") {
		vdis <- compute_ludo_dtce(talpha)
		M <- put_vector_into_matrix(ncol(talpha), vdis, colnames(talpha))
	}
	else if (method=="manhattan")
		M <- as.matrix(dist(t(talpha), method="manhattan"))
	else if ("wrong distance method")
	return(M)
}

get_tree <- function(distance_matrix, racine)
{
	tree <- nj(distance_matrix)
	tree <- ladderize(root (tree, racine), right=T)
	return(tree)
}

get_tree_edges <- function(tree)
# select terminal branches
{
	edges <-  which(tree$edge[,2] <= Ntip(tree))
	edge_nb <- tree$edge[edges,2]
	tree_edges <- list(edges=edges, edge_nb=edge_nb)
	return(tree_edges)
}

apply_colors_2_edges <- function(tree, tree_edges, lane, col_races, col_lib)
{
	edge_colors <- rep("black", length(tree$edge[,1]))
	if (lane) 
		edge_colors[tree_edges$edges] <- col_lib[tree_edges$edge_nb] 
	else
		edge_colors[tree_edges$edges] <- col_races[tree_edges$edge_nb]
	return(edge_colors)
}

plot_tree <- function(tree, col_races, edge_colors, lane, libraries, races_uniq, col_lib_raw, race_colo_raw, categ, nbaits)
{
	plot(tree, show.tip.label = T, tip.color=col_races, edge.color=edge_colors, edge.width=1.5, root.edge=F, cex = 0.7, font=4)
#	plot(tree, show.tip.label = T, tip.color=col_races, edge.color=edge_colors, edge.width=1.5, root.edge=F, cex = 0.45, font=4)

	if (lane) 
		legend("topright", legend=c(libraries, races_uniq), col = col_lib_raw, lty = c(rep(1, length(libraries)), rep(-1, length(race_colo_raw))), lwd=2, bg = 'gray92', cex=0.8, text.col=c(rep("black", length(libraries)), race_colo_raw))
	else 
		legend("topright", legend=races_uniq, col = race_colo_raw, lwd=2, bg = 'gray92', cex=0.8)
	title(paste("NJ tree of CNV distribution across aphid individuals\n", categ, " genes (", nbaits, " baits)", sep=""), cex.main=1.5, font=4)
}

processing_NJ <- function(i_categ, gene_categ, alpha_matrix, method, racine, lane, col_races, col_lib, libraries, races_uniq, col_lib_raw, race_colo_raw, baits)
{

	categ <- gene_categ[i_categ]
	talpha <- define_alpha_matrix(alpha_matrix, gene_categ, i_categ, baits)
	nbaits <- nrow(talpha)
	distance_matrix <-  compute_dstce(talpha, method)
	tree <- get_tree(distance_matrix, racine)
	
	tree_edges <- get_tree_edges(tree)
	edge_colors <- apply_colors_2_edges(tree, tree_edges, lane, col_races, col_lib)
	plot_tree(tree, col_races, edge_colors, lane, libraries, races_uniq, col_lib_raw, race_colo_raw, categ, nbaits)
	
	print(categ)
	print(dim(talpha))

}

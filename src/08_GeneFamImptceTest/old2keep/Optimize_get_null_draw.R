get_null_draw <- function(ite, Gns, PMTs, bait_names, n_info, info_TargGene0, verbose=F){

    ind <- match(c(Gns, PMTs), bait_names)
    imp <- rep(0, length(bait_names))
    names(imp) <- bait_names
    imp [ind] <- sample(seq(n_info), n_info)   # should stay ~240 genes with some importance
    imp <- sort(imp, decreasing=T)  # has to be sorted for the function 'get_BestBaitperContig'

    # 3) select best bait per contig
    imp1 <- imp[imp!=0]
    if (verbose) cat(imp1, sep="\n", file="importance.txt")
    all_targ <- sapply(names(imp1), collapse_elements, sep="_", what=1:3, colla="_")
    ind <- PrePro_findIndex(all_targ, info_TargGene0$NewTargetName)
#~    info_TargGene <- info_TargGene[info_TargGene$NewTargetName%in%all_targ,]
    info_TargGene <- info_TargGene0[ind,]

    if (verbose) cat(names(imp1), sep="\n", file="genes.txt")
    bestBait_PerContig <- as.character(get_1baitPerContig(info_TargGene, imp1, verbose=F))
    if (verbose) print(bestBait_PerContig)
    ind <- PrePro_findIndex(bestBait_PerContig, names(imp))
    
    # 4) keep only a rankable importance for the best bait best bait per contig
        # attribute an EQUAL low importance of 0 for all other baits
        # corresponds to step 3.1) in the main script
    imp[-ind] <- 0
        # kepp on ly the 'kept' first loci
    bad <- order(imp, decreasing=T)[(kept+1):length(bait_names)]
    imp[bad] <- 0

    # 5) compute sum of ranks per grp
    df_fc <- data.frame(rnk = rank(imp), grp = g)
    res <- aggregate(rnk ~ grp,df_fc,sum)
    v <- res$rnk
    names(v) <- res$grp
    if (ite%%100==0) print(paste(ite, "iterations done"))
    return(v)
}

get_null_draw2 <- function(ite, Gns, PMTs, bait_names, n_info, info_TargGene_fil, verbose=F){

    ind <- match(c(Gns, PMTs), bait_names)
    imp <- rep(0, length(bait_names))
    names(imp) <- bait_names
    imp [ind] <- sample(seq(n_info), n_info)   # should stay ~240 genes with some importance
    imp <- sort(imp, decreasing=T)  # has to be sorted for the function 'get_BestBaitperContig'

    # 3) select best bait per contig
    imp1 <- imp[imp!=0]
    if (verbose) cat(imp1, sep="\n", file="importance.txt")
    all_targ <- sapply(names(imp1), collapse_elements, sep="_", what=1:3, colla="_")
    info_TargGene <- get_info_AllTargGene(info_TargGene_fil, all_targ)

    if (verbose) cat(names(imp1), sep="\n", file="genes.txt")
    bestBait_PerContig <- as.character(get_1baitPerContig(info_TargGene, imp1, verbose=F))
    if (verbose) print(bestBait_PerContig)
    ind <- PrePro_findIndex(bestBait_PerContig, names(imp))
    
    # 4) keep only a rankable importance for the best bait best bait per contig
        # attribute an EQUAL low importance of 0 for all other baits
        # corresponds to step 3.1) in the main script
    imp[-ind] <- 0
        # kepp on ly the 'kept' first loci
    bad <- order(imp, decreasing=T)[(kept+1):length(bait_names)]
    imp[bad] <- 0

    # 5) compute sum of ranks per grp
    df_fc <- data.frame(rnk = rank(imp), grp = g)
    res <- aggregate(rnk ~ grp,df_fc,sum)
    v <- res$rnk
    names(v) <- res$grp
    if (ite%%100==0) print(paste(ite, "iterations done"))
    return(v)
}


# set up the data set
library(microbenchmark)
microbenchmark(
                get_null_draw(1, sims2[[1]]$Gns, sims2[[1]]$PMTs, bait_nam, n_info=length(gini_gene_rf), info_TargGene0=read.delim(INFO_TARGENE_FILE), verbose=F),
                get_null_draw2(1, sims2[[1]]$Gns, sims2[[1]]$PMTs, bait_nam, n_info=length(gini_gene_rf), info_TargGene=INFO_TARGENE_FILE, verbose=F)
)

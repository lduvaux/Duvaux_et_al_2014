test_CDD <- function(x, gd_genes, GLMtab_all1, tab_best_obs, v_races){

    gene <- gd_genes[x]
    sub_GLM <- subset(GLMtab_all1, GLMtab_all1[,"Gene"]==gene)
    ind <- grep(paste(gene, "_", sep=""), tab_best_obs$Subtarget)
    vec_best_obs <- tab_best_obs[ind, v_races]
    
    best_disc_race_ind <- which.max(vec_best_obs[,])
    best_disc_race <- colnames(vec_best_obs)[best_disc_race_ind]

    CDD_best_ind <- grep(best_disc_race, rownames(sub_GLM))
    CDD_best <- sub_GLM$Duplication[CDD_best_ind]

    CDD_test_all <- "3_CpDup"%in%sub_GLM$Duplication

    v <- c(as.character(CDD_best), CDD_test_all)
    return(v)
}

get_info_rking <- function(tab_best_obs, GLMtab_all1, v_races){
    subtargs <- as.character(tab_best_obs$Subtarget)
    genes <- sapply(subtargs, collapse_elements, what=1:2)
    info_CDD <- matrix(data=NA, nrow=nrow(tab_best_obs), ncol=2, dimnames=list(rownames(tab_best_obs), c("CDD_best", "CDD_test_all")))
    
    inside <- genes%in%GLMtab_all1$Gene
    gd_genes <- genes[inside]

    fetch_CDD <- t(sapply(seq(gd_genes), test_CDD, gd_genes, GLMtab_all1, tab_best_obs, v_races))
    info_CDD[inside,] <- as.character(fetch_CDD)

        # N_CDD & N_CDD_50
    N_CDD <- table(fetch_CDD)[3]
    N_CDD_50 <- table(fetch_CDD[1:50,])[3]

    # rks CDD in all the final subtargets
    rks <- nrow(tab_best_obs):1
    ind <- ifelse(info_CDD[,"CDD_best"]!="3_CpDup" | is.na(info_CDD[,"CDD_best"]), F, T)
    sum_rk <- sum(rks[ind])

    # rks CDD in all the 50 best subtargets
    rks_50 <- 50:1
    ind <- ifelse(info_CDD[1:50,"CDD_best"]!="3_CpDup" | is.na(info_CDD[1:50,"CDD_best"]), F, T)
    sum_rk_50 <- sum(rks_50[ind])

    ll <- list(N_CDD=N_CDD, N_CDD_50=N_CDD_50, sum_rank_CDD=sum_rk, sum_rank_CDD_50=sum_rk_50, info_CDD=info_CDD)
}

get_real_name <- function(target, tab_target){

    tab_target <- subset(tab_target, NewTargetName==target)

    res <- ifelse(length(grep("ApisSNMP|Control", tab_target[1,"NewTargetName"]))==1, collapse_elements(tab_target[1, "NewTargetName"], what=1:2), tab_target[1,1])
    return(res)
}

get_CNmedian_race <- function(vec, l_races){
    res <- sapply(l_races, function(x) median(vec[x]))
    return(res)
}


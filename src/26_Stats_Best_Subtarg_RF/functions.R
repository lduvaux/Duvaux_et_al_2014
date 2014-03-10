test_CDD <- function(x, gd_genes, GLMtab_all1, tab_best_obs){

    gene <- gd_genes[x]
    sub_GLM <- subset(GLMtab_all1, GLMtab_all1[,"Gene"]==gene)
    ind <- grep(paste(gene, "_", sep=""), tab_best_obs$Subtarget)
    vec_best_obs <- tab_best_obs[ind,4:11]
    
    best_disc_race_ind <- which.max(vec_best_obs[,])
    best_disc_race <- colnames(vec_best_obs)[best_disc_race_ind]

    CDD_best_ind <- grep(best_disc_race, rownames(sub_GLM))
    CDD_best <- sub_GLM$Duplication[CDD_best_ind]

    CDD_test_all <- "3_CpDup"%in%sub_GLM$Duplication

    v <- c(as.character(CDD_best), CDD_test_all)
    return(v)
}

get_info_rking <- function(tab_best_obs, GLMtab_all1){
    subtargs <- as.character(tab_best_obs$Subtarget)
    genes <- sapply(subtargs, collapse_elements, what=1:2)
    info_CDD <- matrix(data=NA, nrow=nrow(tab_best_obs), ncol=2, dimnames=list(rownames(tab_best_obs), c("CDD_best", "CDD_test_all")))
    
    inside <- genes%in%GLMtab_all1$Gene
    gd_genes <- genes[inside]

    fetch_CDD <- t(sapply(seq(gd_genes), test_CDD, gd_genes, GLMtab_all1, tab_best_obs))
    info_CDD[inside,] <- as.character(fetch_CDD)
    N_CDD <- table(fetch_CDD)[3]
    rks <- 111:1
#~    ind0 <- info_CDD[,"CDD_best"]=="3_CpDup" 
    ind <- ifelse(info_CDD[,"CDD_best"]!="3_CpDup" | is.na(info_CDD[,"CDD_best"]), F, T)
    sum_rk <- sum(rks[ind])

    ll <- list(N_CDD=N_CDD, sum_rank_CDD=sum_rk, info_CDD=info_CDD)
}













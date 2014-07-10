set_infoSubtarg <- function(prev_tab){
    t0 <- matrix(ncol=10, nrow=nrow(prev_tab), data=as.numeric(prev_tab[,2:ncol(prev_tab)]))
    colnames(t0) <- colnames(prev_tab)[2:ncol(prev_tab)]
    tab <- data.frame(prev_tab[,1], stringsAsFactors=F)
    colnames(tab) <- "Subtarget"
    tab <- cbind(tab, t0)
    return(tab)
}

# check CDD for best discriminating genes
set_alpha_mat <- function(alpha_mat, last_rf_mat){
    alpha_mat_round <- PrePro_roundToZeroFive(alpha_mat)
    ind <- match(rownames(last_rf_mat), colnames(alpha_mat_round))
    mat <- alpha_mat_round[,ind]
    return(mat)
}

set_submatrix <- function(race, gene, alp_mat){
    x_ind <- grep(paste("^", gene, "_", sep=""), rownames(alp_mat))
    y_ind <- grep(paste("^", race, "_", sep=""), colnames(alp_mat))
    mat <- alp_mat[x_ind, y_ind]
}

check_CDD_ind <- function(v_alpha_ij){
# check if a loci is completely duplicated for an idindividual
    # v_alpha_i: vector of alpha values of a locus i for individual j

    test <-
    ifelse(length(unique(v_alpha_ij))==1 & v_alpha_ij[1]==1,
        F,
            ifelse(length(unique(v_alpha_ij))==1 & v_alpha_ij[1]!=1,
            T, F)
    )
    return(test)
}

get_CDD_MDR <- function(race, gene, alp_mat, doprint=F){
    if (doprint) print(gene)
    submat <- set_submatrix(race, gene, alp_mat)
    if (class(submat)=="matrix"){
        mat_CDD <- apply(submat, 1, check_CDD_ind)
        res <- ifelse(T%in%mat_CDD, T, F)
    } else res <- NA
    return(res)
}

# count median Copy number per race
get_MCN <- function(race, alp_mat, doprint=F){
    submat <- alp_mat[grep(paste("^", race, "_", sep=""), rownames(alp_mat)),]
    res <- apply(submat, 2, median)
    return(res)
}



















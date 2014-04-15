get_vs <- function(vec_val, vec_gp){
    gps <- unique(vec_gp)
    vec_var <- sapply(gps, function(x) var(vec_val[which(vec_gp==x)]))
    res <- mean(vec_var)
    return(res)
}

get_Vst <- function(vec_val, vec_gp){
    vt <- var(vec_val)
    vs <- get_vs(vec_val, vec_gp)
    vst <- (vt-vs)/vt
    return(vst)
}


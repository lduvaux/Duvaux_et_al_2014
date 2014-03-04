test_singleton <- function(vec){
    tab <- table(vec !=1)
    res <- ifelse(tab[2]==1, T, F)
    return(res)
}

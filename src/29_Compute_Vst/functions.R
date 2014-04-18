get_vs <- function(vec_val, vec_gp){
    gps <- unique(vec_gp)
    vec_var <- sapply(gps, function(x) var(vec_val[which(vec_gp==x)]))
    res <- mean(vec_var)
    return(res)
}

get_Vst <- function(vec_val, vec_gp, log_2=F){
    if (log_2) vec_val <- log2(vec_val+0.0001)
    vt <- var(vec_val)
    vs <- get_vs(vec_val, vec_gp)
    vst <- (vt-vs)/vt
    return(vst)
}

get_max_contig <- function(contig, vector_vst, vec_cont){
    ind <- which(vec_cont==contig)
    sub_vector_vst <- vector_vst[ind]
    res <- which.max(sub_vector_vst)
    res <- sub_vector_vst[res]
    return(res)    
}

get_all_maxVst <- function(vector_vst, vec_cont){
    contigs <- unique(vec_cont)
    vec_maxVst <- sapply(contigs, get_max_contig, vector_vst, vec_cont)
    return(vec_maxVst)
}

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

plot_dble_hist <- function(famille, list_ind, xlim=c(-1.8, 1), ylim=c(0,2), breaks=seq(-1.8, 1, by=0.1), alpha, color){

    main <- paste("Vst Distribution of control\nand ", famille, " genes", sep="")
    hist(vec_vst[list_ind[["Control"]]],freq=F, xlim=xlim, ylim=ylim, breaks=breaks, xlab="Vst", main=main)
    par(new=T)
    colo <- makeTransparent(color, alpha)
    hist(vec_vst[list_ind[[famille]]],freq=F, xlim=xlim, ylim=ylim, xlab="", ylab="",col=colo, border=colo, breaks=breaks, main="")
}

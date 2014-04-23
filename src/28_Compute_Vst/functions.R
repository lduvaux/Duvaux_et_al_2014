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

dble_hist <- function(vec_vst, famille, list_ind, xlim, ylim, breaks, alpha, color){
    N_control <- length(list_ind[["Control"]])
    N_other <- length(list_ind[[famille]])
    main <- paste("Vst Distribution of control (N=", N_control, ")\nand ", famille, " genes (N=", N_other, ")", sep="")
    hist(vec_vst[list_ind[["Control"]]],freq=F, xlim=xlim, ylim=ylim, breaks=breaks, xlab="Vst", main=main)
    par(new=T)
    colo <- makeTransparent(color, alpha)
    hist(vec_vst[list_ind[[famille]]],freq=F, xlim=xlim, ylim=ylim, xlab="", ylab="",col=colo, border=colo, breaks=breaks, main="")
}

plot_dble_hist <- function(vst_val, famm=c("Gr", "Or", "P450"), list_fam, colos=c("blue", "purple", "green"), x_lim=c(-1.8, 1), y_lim=c(0,2), brks=seq(-1.8, 1, by=0.1), alpha=70, nam_plot){
    pdf(nam_plot)
    layout(matrix(1:4, nrow=2, ncol=2, byrow=T))
    sapply(1:length(famm), function(x) dble_hist(vst_val, famm[x], list_fam, xlim=x_lim, ylim=y_lim, breaks=brks, color=colos[x], alpha=alpha))
    dev.off()
}

get_segments_gene <- function(gen, v_genes, tab_alpha){
#~    print(gene)
    ind <- which(v_genes==gen)
    tab <- tab_alpha[ind,]
    if (length(ind)==1)
        res <- t(as.matrix(tab))
    else
        res <- as.matrix(unique(tab))
    return(res)
}

get_alpha_seg <- function(alpha_mat){
    vec_genes <- sapply(rownames(alpha_mat), collapse_elements, what=1:2)
    genes <- unique(vec_genes)
    lis_alpha_seg <- sapply(genes, get_segments_gene, vec_genes, alpha_mat)
    n_seg <- sapply(lis_alpha_seg, nrow)
    genes2 <- rep(genes, n_seg)
    vec_alpha_seg <- unlist(lis_alpha_seg)
    m_alph_seg <- matrix(data=vec_alpha_seg, ncol=ncol(alpha_mat), nrow=length(vec_alpha_seg)/ncol(alpha_mat), byrow=T, dimnames=list(genes2, colnames(alpha_mat)))
    return(m_alph_seg)
}

compute_gene_Vst <- function(m_alph_seg, races, log_2=F){
    v_genes <- rownames(m_alph_seg)
        # average Vs per gene
    if (log_2)
        vec_Vs <- apply(log2(m_alph_seg+0.0001), 1, get_vs, races)
    else
        vec_Vs <- apply(m_alph_seg, 1, get_vs, races)
    mean_Vs <- multi_mean(vec_Vs, v_genes)
        # average Vt per gene
    if (log_2)
        vec_Vt <- apply(log2(m_alph_seg+0.0001), 1, var)
    else
        vec_Vt <- apply(m_alph_seg, 1, var)
    mean_Vt <- multi_mean(vec_Vt, v_genes)

    Vst <- (mean_Vt - mean_Vs)/mean_Vt
    return(Vst)
}

multi_mean <- function(vec, v_grps){
    grps <- unique(v_grps)
    v_means <- sapply(grps, function(x) mean(vec[which(v_grps==x)]))
    names(v_means) <- grps
    return(v_means)
}

sort_name_per_categ <- function(vect){
    vec_categ <- sapply(vect, get_elements)
    categs <- unique(vec_categ)
    list_ind <- lapply(categs, function(x) which(vec_categ==x))
    names(list_ind) <- categs
    return(list_ind)
}

get_boxplot <- function(vec_vst, families, main=NULL){

    if (is.null(main))
        main <- deparse(substitute(vec_vst))

    vec_families <- sapply(names(vec_vst), get_elements)
    good <- vec_families%in%families
    vec_families <- vec_families[good]
    vec_vst_new <- vec_vst[good]
    res <- boxplot(vec_vst_new~vec_families, main=main)
    return(res)
}

test_dif_NoParam <- function(v_factors, sub_factors, v_values, toprint, wilcox=T){
    cat("\n")
    cat("# ", toprint, sep="")
    good <- v_factors%in%sub_factors
    sub_v_factors <- as.factor(v_factors[good])
    sub_val <- v_values[good]
    if (wilcox)
        test <- wilcox.test(sub_val ~ sub_v_factors)
    else
        test <- kruskal.test(sub_val ~ sub_v_factors)
    return(test)
}




















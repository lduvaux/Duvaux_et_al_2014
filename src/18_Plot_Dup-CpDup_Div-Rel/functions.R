draw_plot <- function(obsprob, colo, categ, yli=c(0, 0.6), xla="", yla="", sep_ind="def", cex_lab=1.2, m_gp=c(4, 1, 0))
{

#~    if (length(categ)!=length(obsprob))

    lg <- length(obsprob)
    v_space <- rep(0.2, lg)
    if (sep_ind[1]=="def")
        sep_ind <- lg/2+1
    if (length(sep_ind)==0)
        stop("No index for separations given")
    v_space[sep_ind] <- 1.5

    par(mar=c(7, 6,4,2))
    mp <- barplot(obsprob, ylim=yli, axisnames = FALSE, space=v_space, col=colo, xlab=xla, ylab="", font.lab=2, cex.lab=cex_lab, mgp=m_gp)
    mtext(yla, side=2, line=2.5, font=2, cex.lab=cex_lab)
    text(mp, par("usr")[3], labels = categ, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
}

get_exdon_number <- function(gene, targets){
    ind <- grep(paste(gene, "_", sep=""), targets)
    res <- length(ind)
    if (res==0) stop(paste("No exons present for gene", gene))
    return(res)
}

get_length <- function(gene, GLMtab_all2){
    ind <- pmatch(gene, GLMtab_all2$Gene)
    lg <- GLMtab_all2$LnExonLength[ind[1]]
    return(lg)
}



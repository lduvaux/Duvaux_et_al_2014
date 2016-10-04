draw_plot <- function(obsprob, colo, categ, yli=c(0, 0.6), xla="", yla="", sep_ind="def", m_gp=c(4, 1, 0), cex_lab=1.2)
{

#~    if (length(categ)!=length(obsprob))

    lg <- length(obsprob)
    v_space <- rep(0.2, lg)
    if (sep_ind[1]=="def")
        sep_ind <- lg/2+1
    if (length(sep_ind)==0)
        stop("No index for separations given")
    v_space[sep_ind] <- 1.5
    
    mp <- barplot(obsprob, ylim=yli, axisnames = FALSE, space=v_space, col=colo, xlab=xla, ylab="", mgp=m_gp, font.lab=2, cex.lab=cex_lab)
    text(mp, par("usr")[3], labels = categ, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
    mtext(yla, side=2, line=2.5, font=2, cex.lab=cex_lab)
}

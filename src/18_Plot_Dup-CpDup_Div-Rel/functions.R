draw_plot <- function(obsprob, colo, categ, yli=c(0, 0.6), yla="", sep_ind="def")
{

#~    if (length(categ)!=length(obsprob))

    lg <- length(obsprob)
    v_space <- rep(0.2, lg)
    if (sep_ind[1]=="def")
        sep_ind <- lg/2+1
    if (length(sep_ind)==0)
        stop("No index for separations given")
    v_space[sep_ind] <- 1.5
    
    mp <- barplot(obsprob, ylim=yli, axisnames = FALSE, space=v_space, col=colo, ylab=yla)
    text(mp, par("usr")[3], labels = categ, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
}

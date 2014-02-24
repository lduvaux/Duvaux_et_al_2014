draw_plot <- function(x, y, colo, categ, yli=c(0, 1), yla="", lab)
{

    par(mar=c(7, 6,4,2))
    mp <- boxplot(y~x, ylim=yli, xaxt="n", col=colo, ylab=yla)
    text(seq(mp$names), par("usr")[3], labels = lab, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)

}

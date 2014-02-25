draw_plot <- function(x, y, colo, categ, yli=c(0, 1), yla="", lab, main="", new_mar=c(7, 6,4,2), att, adjust=0)
{

    par(mar=new_mar)
    mp <- boxplot(y~x, ylim=yli, xaxt="n", col=colo, ylab=yla, main=main, at=att)
    text(att-adjust, par("usr")[3], labels = lab, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)

    return(mp)
}

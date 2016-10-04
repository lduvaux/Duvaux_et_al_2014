draw_plot <- function(x, y, colo, categ, yli=c(0, 1), xla="", yla="", lab, main="", new_mar=c(7, 6,4,2), att, adjust=0, cex_lab=1.2, xline=4)
{

    par(mar=new_mar)
    mp <- boxplot(y~x, ylim=yli, xaxt="n", col=colo, xlab="", ylab="", main=main, at=att, font.lab=2, cex.lab=cex_lab)
    mtext(xla, side=1, line=xline, font=2, cex.lab=cex_lab)
    mtext(yla, side=2, line=2.5, font=2, cex.lab=cex_lab)
    text(att-adjust, par("usr")[3], labels = lab, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)

    return(mp)
}

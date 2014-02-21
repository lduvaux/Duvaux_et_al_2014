draw_plot <- function(obs_prob, colo, categ, yli=c(0, 0.6))
{

    v_space <- rep(0.2, 8); v_space[5] <- 1.5
    mp <- barplot(obs_prob, ylim=yli, axisnames = FALSE, space=v_space, col=colo)
    text(mp, par("usr")[3], labels = categ, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
    text(x=c(2.5, 9), y=c(0.2, 0.5), labels=c("Non\ntruncated", "Truncated"), lheight=1.5)
}

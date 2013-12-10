Draw_prediction <- function(fitted_mdl, tab_data)
{

    tmpdat <- tab_data[,c("LnExonLength", "Family", "Race", "trimmed", "Duplication")]

    # calculate predicted probabilities and store in a list
    LnExon <- with(tmpdat, seq(from = min(LnExonLength), to = max(LnExonLength), length.out = 100))
    biprobs <- lapply(levels(tab_data$Family), function(fam) {
        tmpdat$Family[] <- fam
        lapply(LnExon, function(j) {
            tmpdat$LnExonLength <- j
            predict(fitted_mdl, newdata = tmpdat, type = "response")})
    })

    test_predict <- predict(fitted_mdl, newdata = tmpdat, type = "response")

    # get means and quartiles for all jvalues for each level of Family
    plotdat2 <- lapply(biprobs, function(X) {
        temp <- t(sapply(X, function(x) {
            c(M=mean(x), quantile(x, c(.25, .75))) }))
    
        temp <- as.data.frame(cbind(temp, LnExon))
        colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "LnExonLength")
        return(temp)})

    # collapse to one data frame
    plotdat2 <- do.call(rbind, plotdat2)

    # add families
    plotdat2$Family <- factor(rep(levels(tab_data$Family), each = length(LnExon)))

    # graph it
    ggplot(plotdat2, aes(x = LnExonLength, y = PredictedProbability)) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Family), alpha = .15) +
        geom_line(aes(colour = Family), size = 2) + ylim(c(0, 1)) +
        facet_wrap(~  Family)
}

drw_pred <- function(fitted_mdl, tab_data)
{
    fac0 <- with(tab_data, interaction(Dup, Family))
    fac <- tab_data$Family
    # 1) prediction
    pp <- fitted(fitted_mdl)

    # 2) raw data
    fq_Gn <- table(tab_data$Family)
    fq_Dup0 <- table(fac0)
    fq_Dup <- fq_Dup0[seq(2, length(fq_Dup0), by=2)]
    PrDup_data <- fq_Dup/fq_Gn

    # 3) drawing
    vcol <- rep(c("grey80", "red", "green", "blue"), each=1)
    vcol2 <- rep(c("black", "red", "green", "blue"), each=1)
    bp <- boxplot(pp~fac, xaxt="n", lwd=1.5, col=vcol)
    text(seq(bp$names), par("usr")[3], labels = bp$names, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9, col=vcol2, font=2)
    points(1:4, PrDup_data, lwd=2, col=1, pch=21, bg=5)
}

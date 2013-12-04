source("./10_GeneralStats_functions.R")
source("./glm_common.R")

###############
set_GLMtab <- function(tab0, NonTrim_only=F, covar)
{
	if (covar!="LnExonLength" & covar!="ratioLength")
		stop("Covariate is not given")

	if (NonTrim_only)
		tab0 <- tab0[tab0$trimmed=="Yes",]

	if (covar=="LnExonLength") {
		tab <- with(tab0, data.frame(Polymorphism, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), LnExonLength=log(TotExonLength)))
        bad <- which(is.na(tab$LnGeneLength) | is.na(tab$LnExonLength))}

	if (covar=="ratioLength") {
		tab <- with(tab0, data.frame(Polymorphism, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), ratioLength=qlogis(ratioLength)))
		bad <- which(is.na(tab$LnGeneLength) | is.na(tab$ratioLength))}

    rownames(tab) <- rownames(tab0)
    nambad <- rownames(tab)[bad]
    tab <- tab[-bad,]

	print(paste ("The following locus have been removed from the analysis:", nambad, sep=" "))
	return(tab)
}

############### draw predictions
drw_pred <- function(fitted_mdl, tab_data)
{
    fac <- with(tab_data, interaction(Family, trimmed))
    fac1 <- with(tab_data, interaction(Polymorphism, Family, trimmed))
    fac2 <- interaction(rep(c("1_NoDup", "2_PtDup", "3_CpDup"), length(fac)), rep(fac, each=ncol(tab_pred)))

    # 1) prediction
    tab_pred <- fitted(fitted_mdl)
    vec_pred <- as.vector(t(tab_pred))
    
    # 2) raw data
    tb <- table(fac1)
    fq_tb <- tb/sum(tb)
    densities <- rep(sapply(seq(1, length(fq_tb), by=ncol(tab_pred)), function (x) sum(fq_tb[x:(x+(ncol(tab_pred)-1))])), each=ncol(tab_pred))
    fq_raw_data <- fq_tb/densities

    # 3) drawing
#~    for (i in seq(ncol(tab_pred))) {
    par(mar=c(7, 4, 4, 2))
#~    vlty <- rep(c(1,2,3,4), each=3)
    vcol <- rep(c("grey85", "red", "green", "blue"), each=3)
    bp <- boxplot(vec_pred~fac2, xaxt="n", lwd=1.5, col=vcol)
    abline(v=12.5, lwd=2, lty=3)
    text(seq(bp$names), par("usr")[3], labels = bp$names, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
    points(1:length(bp$n), fq_raw_data, col="black", lwd=3)
#~    }
}

####################
pairs_glm <- function(model, data_tab)
{
	pairs(model, data_tab, upper.panel = panel.smooth, lower.panel = panel.cor, diag.panel =panel.hist)
}









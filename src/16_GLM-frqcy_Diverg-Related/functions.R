source("./10_GeneralStats_functions.R")
source("./glm_common.R")

############### add fields to Pre_GLMtab0
add_Phylog_lvl <- function(tab, clusters)
{
    Phylog_lvl <- rep('divergent', nrow(tab))
    ind <- tab$Race%in%clusters$related
    Phylog_lvl[ind] <- "related"
    tab <- cbind(tab, Phylog_lvl)
    return(tab)
}

add_CpDupField <- function (tab)
{
    CpDup <- tab$Fqcy_CpDup!=0
    tab <- cbind(tab, CpDup)
    return(tab)
}

add_raceField <- function (tab)
{
    v <- rownames(tab)
    race <- sapply(v, function (x) unlist(strsplit(x, "_"))[3])
    tab <- cbind(tab, race)
    return(tab)
}

############### set_GLMtab
set_GLMtab <- function(tab0, NonTrim_only=F, covar)
{
	if (covar!="LnExonLength" & covar!="ratioLength")
		stop("Covariate is not given")

	if (NonTrim_only)
		tab0 <- tab0[tab0$trimmed=="Yes",]

	if (covar=="LnExonLength") {
		tab <- with(tab0, data.frame(Fqcy_all, Fqcy_CpDup, CpDup, Phylog_lvl, race, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), LnExonLength=log(TotExonLength)))
		bad <- which(is.na(tab$LnGeneLength) | is.na(tab$LnExonLength))}

	if (covar=="ratioLength") {
		tab <- with(tab0, data.frame(Fqcy_all, Fqcy_CpDup, CpDup, Phylog_lvl, race, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), ratioLength=qlogis(ratioLength)))
		bad <- which(is.na(tab$LnGeneLength) | is.na(tab$ratioLength))}

    rownames(tab) <- rownames(tab0)
    nambad <- rownames(tab)[bad]
    tab <- tab[-bad,]

	print(paste ("Because of Nas, the following locus have been removed from the analysis:", nambad, sep=" "))
	return(tab)
}

############### draw predictions 
drw_pred <- function(fitted_mdl, tab_data)
{
    fac <- with (tab_data, interaction(CpDup, Family, trimmed))
#~    fac <- with (tab_data, interaction(CpDup, Family))

    # 1) prediction
#~    tab_pred <- fitted(fitted_mdl)
#~    fac2 <- interaction(rep(c("1_NoDup", "2_PtDup", "3_CpDup"), length(fac)), rep(fac, each=ncol(tab_pred)))
#~    vec_pred <- as.vector(t(tab_pred))
    
    # 2) raw data
#~    tb <- table(fac1)
#~    fq_tb <- tb/sum(tb)
#~    densities <- rep(sapply(seq(1, length(fq_tb), by=ncol(tab_pred)), function (x) sum(fq_tb[x:(x+(ncol(tab_pred)-1))])), each=ncol(tab_pred))
#~    fq_raw_data <- fq_tb/densities

    # 3) drawing
#~    for (i in seq(ncol(tab_pred))) {
    par(mar=c(7, 4, 4, 2))
#~    vlty <- rep(c(1,2,3,4), each=3)
    vcol <- rep(c("grey80", "red", "green", "blue"), each=2)
    vcol2 <- rep(c("black", "red", "green", "blue"), each=2)
    bp <- boxplot(tab_data$Fqcy_all~fac, xaxt="n", lwd=1.5, col=vcol)
    abline(v=seq(2, 14, by=2)+0.5, lwd=c(rep(1.5,3), 2, rep(1.5,3)), lty=c(rep(3,3),2,rep(3,3)))
    text(seq(bp$names), par("usr")[3], labels = bp$names, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9, col=vcol2, font=c(1,2))
#~    }
}


####################
pairs_glm <- function(model, data_tab)
{
	pairs(model, data_tab, upper.panel = panel.smooth, lower.panel = panel.cor, diag.panel =panel.hist)
}







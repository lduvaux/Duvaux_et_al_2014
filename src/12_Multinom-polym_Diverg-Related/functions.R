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
		tab <- with(tab0, data.frame(Duplication, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), LnExonLength=log(TotExonLength)))
        bad <- which(is.na(tab$LnGeneLength) | is.na(tab$LnExonLength))}

	if (covar=="ratioLength") {
		tab <- with(tab0, data.frame(Duplication, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), ratioLength=qlogis(ratioLength)))
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
    fac1 <- with(tab_data, interaction(Duplication, Family, trimmed))


    # 1) prediction
    tab_pred <- fitted(fitted_mdl)
    fac2 <- interaction(rep(c("1_NoDup", "2_PtDup", "3_CpDup"), length(fac)), rep(fac, each=ncol(tab_pred)))
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
    vcol <- rep(c("grey80", "red", "green", "blue"), each=3)
    vcol2 <- rep(c("black", "red", "green", "blue"), each=3)
    bp <- boxplot(vec_pred~fac2, xaxt="n", lwd=1.5, col=vcol)
    abline(v=c(3.5, 6.5, 9.5, 12.5, 15.5, 18.5, 21.5), lwd=c(rep(1.5,3), 1.5,rep(1.5,3)), lty=c(rep(3,3),2,rep(3,3)))
    text(seq(bp$names), par("usr")[3], labels = bp$names, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9, col=vcol2, font=c(1,3,2))
    points(1:length(bp$n), fq_raw_data, col=1, lwd=2, pch=21, bg=6)
#~    }
}

####################
pairs_glm <- function(model, data_tab)
{
	pairs(model, data_tab, upper.panel = panel.smooth, lower.panel = panel.cor, diag.panel =panel.hist)
}

############
Output_glm_res <- function(Sum, outfile, sep="\t")# for nnet object only
{
	z <- Sum$coefficients/Sum$standard.errors
	p <- (1 - pnorm(abs(z), 0, 1))*2
	tab <- round(rbind(Sum$coefficients, Sum$standard.errors, z, p),5)
	nlevel <- length(Sum$lev)
	if (nlevel==2)
		rownames(tab) <- paste(Sum$lev[2], c("coef", "sd", "z", "p-val"), sep="_")
	else
		rownames(tab) <- paste(rownames(tab), rep(c("coef", "sd", "z", "p-val"),each=nlevel-1), sep="_")
	tab  <- redo_table(t(tab))
	write.table(tab, file=outfile, sep=sep, row.names=F, quote=F)
}

redo_table <- function(tab)
{
	return(tab <- cbind(Model_Terms=rownames(tab), tab))
}

############
get_best_nnet_multinom <- function(dredge_object, Delta=5)
{
	good_mod0 <- subset(dredge_object, delta < Delta)
	good <- !c(F, sapply(2:length(good_mod0$df), function(x) T%in%(good_mod0$df[x]>good_mod0$df[1:(x-1)])))
	good_mod <- subset(good_mod0, good)
	df_order <- order(good_mod$df)
	if (df_order[1]==1)
		res <- attributes(good_mod)$row.names[1]
	else {
		df_order2 <- df_order[-which(df_order==1)]
		inc <- 1
		mod1 <- get.models(test_trimmed, inc)[[1]]
		for (ite in 2:nrow(good_mod))
		{
			mod2 <- get.models(good_mod, ite)[[1]]
			lrt_test <- lrtest(mod1, mod2)
			print(paste("P-val:", round(lrt_test$p.value, 4)))
			if (lrt_test$p.value>=0.06) {
				mod1 <- mod2
				inc <- ite}
		}
		res <- attributes(good_mod)$row.names[inc]
	}
	fres <- which(attributes(test_trimmed)$row.names==res)
	return(fres)
}









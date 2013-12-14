source("./10_GeneralStats_functions.R")
source("./glm_common.R")

####################
pairs_glm <- function(model, data_tab)
{
	pairs(model, data_tab, upper.panel = panel.smooth, lower.panel = panel.cor, diag.panel =panel.hist)
}

####################
add_PolymField_fqcy <- function (tab)
{
    Polymorphic <- tab$Fqcy_all!=0
    tab <- cbind(tab, Polymorphic)
    return(tab)
}

add_CpDupField_fqcy <- function (tab)
{
    CpDup <- tab$Fqcy_CpDup!=0
    tab <- cbind(tab, CpDup)
    return(tab)
}

####################
set_GLMtab_fqcy <- function(tab0, NonTrim_only=F)
{
	if (NonTrim_only)
		tab0 <- tab0[tab0$trimmed=="Yes",]

    vzero <- which(tab0$IntronLength==0)
    tab0$IntronLength[vzero] <- tab0$IntronLength[vzero] + 1.1  # remove zero length introns from data sets

    tab <- with(tab0, data.frame(Polymorphic, CpDup, Fqcy_all, Fqcy_CpDup, Phylog_lvl, Race, Gene, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), LnExonLength=log(TotExonLength), LnIntronLength=log(IntronLength), ratioLength=qlogis(ratioLength)))
    bad <- is.na(tab$LnGeneLength)

    rownames(tab) <- rownames(tab0)
    nambad <- rownames(tab)[bad]
    tab <- tab[-bad,]

	print(paste ("Because of Nas, the following locus have been removed from the analysis:", nambad, sep=" "))
	return(tab)
}

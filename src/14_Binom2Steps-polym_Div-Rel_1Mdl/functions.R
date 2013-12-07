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

add_DupField <- function (tab)
{
    Dup <- tab$Polymorphism!="1_NoDup"
    tab <- cbind(tab, Dup)
    return(tab)
}

add_GeneField <- function (tab)
{
    getGeneName <- function(stg) {
        ve <- unlist(strsplit(stg, "_"))[1:2]
        nom <- paste(ve, collapse="_")
        return(nom)}

    vec <- rownames(tab)
    Gene <- sapply(vec, getGeneName)
    tab <- cbind(tab, Gene)
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
		tab <- with(tab0, data.frame(Polymorphism, Dup, Phylog_lvl, Race, Gene, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), LnExonLength=log(TotExonLength)))
		bad <- which(is.na(tab$LnGeneLength) | is.na(tab$LnExonLength))}

	if (covar=="ratioLength") {
		tab <- with(tab0, data.frame(Polymorphism, Dup, Phylog_lvl, Race, Gene, Family, trimmed=as.factor(trimmed), LnGeneLength=log(GeneLength), ratioLength=qlogis(ratioLength)))
		bad <- which(is.na(tab$LnGeneLength) | is.na(tab$ratioLength))}

    rownames(tab) <- rownames(tab0)
    nambad <- rownames(tab)[bad]
    tab <- tab[-bad,]

	print(paste ("Because of Nas, the following locus have been removed from the analysis:", nambad, sep=" "))
	return(tab)
}

############### draw predictions 
drw_pred <- function(dredge_objt, tab_data, delta_dredge, draw_CpDup=F)
{
    # 1) temporary data
    print("Model averaging for prediction")
    fitted_mdl <- model.avg(dredge_objt, subset = delta < delta_dredge, fit=T)
    tmpdat <- tab_data[,c("LnGeneLength","LnExonLength", "Family", "Race", "trimmed", "Dup")]
    v_cont <- c("LnGeneLength","LnExonLength")

#~    g <- list()
    for (i in seq(v_cont)) {
        print(v_cont[i])
        # 2) prepare the prediction table
        p <- predict(fitted_mdl, newdata = tmpdat, se.fit=T, qtype = "response", backtransform=T)
        p <- do.call(cbind, p)
        ymin <- p[,1] - p[,2]*1.96
        ymax <- p[,1] + p[,2]*1.96
        p <- cbind(fit=p[,1], ymin, ymax)


        if (draw_CpDup) {
            CpDup <- ifelse(tab_data$Polymorphism=="3_CpDup",1,0)
            tab_p <- cbind(tmpdat[,c(v_cont[i],"Family")], CpDup, p)
        } else {
            Dup <- as.numeric(tab_data$Dup)
            tab_p <- cbind(tmpdat[,c(v_cont[i],"Family")], Dup, p)
        }
        
        # 3) plot
#~        gplot <- ggplot(tab_p, aes(x = get(v_cont[i]), y = Dup)) + 
#~                geom_point() + # plot the real data
#~                stat_smooth(aes(colour = Family), method = 'glm', family = 'binomial', size = 2) + # plot the fit
#~                facet_wrap( ~ Family)
        y_var <- ifelse(draw_CpDup, "CpDup", "Dup")
        gplot <- ggplot(tab_p, aes_string(x = v_cont[i], y = y_var)) +
                facet_wrap( ~ Family) + # split by families
                geom_point() + # plot the real data
                xlab(v_cont[i]) +
                geom_ribbon(data = tab_p, aes(y = fit, ymin = ymin, ymax = ymax), alpha = 0.25) +    # plot CI
                geom_line(data = tab_p, aes(y = fit, colour = Family)) # plot the fit

        print(gplot)
#~        g[[i]] <- gplot
        }
#~        return(g)
}

####################
pairs_glm <- function(model, data_tab)
{
	pairs(model, data_tab, upper.panel = panel.smooth, lower.panel = panel.cor, diag.panel =panel.hist)
}

############







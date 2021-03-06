source("../10_Functions_GLMMs/functions.R")
source("./glm_common.R")

####################
pairs_glm <- function(model, data_tab)
{
	pairs(model, data_tab, upper.panel = panel.smooth, lower.panel = panel.cor, diag.panel =panel.hist)
}

############### draw predictions 
#~drw_pred <- function(dredge_objt, tab_data, delta_dredge, draw_CpDup=F)
#~{
#~    # 1) temporary data
#~    print("Model averaging for prediction")
#~    fitted_mdl <- model.avg(dredge_objt, subset = delta < delta_dredge, fit=T)
#~    tmpdat <- tab_data[,c("LnGeneLength","LnExonLength", "Family", "Race", "trimmed", "Duplication")]
#~    v_cont <- c("LnGeneLength","LnExonLength")

#~#~    g <- list()
#~    for (i in seq(v_cont)) {
#~        print(v_cont[i])
#~        # 2) prepare the prediction table
#~        p <- predict(fitted_mdl, newdata = tmpdat, se.fit=T, qtype = "response", backtransform=T)
#~        p <- do.call(cbind, p)
#~        ymin <- p[,1] - p[,2]*1.96
#~        ymax <- p[,1] + p[,2]*1.96
#~        p <- cbind(fit=p[,1], ymin, ymax)


#~        if (draw_CpDup) {
#~            CpDuplication <- ifelse(tab_data$Duplicationlication=="3_CpDup",1,0)
#~            tab_p <- cbind(tmpdat[,c(v_cont[i],"Family")], CpDuplication, p)
#~        } else {
#~            Duplication <- as.numeric(tab_data$Duplication)
#~            tab_p <- cbind(tmpdat[,c(v_cont[i],"Family")], Duplication, p)
#~        }
        
#~        # 3) plot
#~#~        gplot <- ggplot(tab_p, aes(x = get(v_cont[i]), y = Dup)) + 
#~#~                geom_point() + # plot the real data
#~#~                stat_smooth(aes(colour = Family), method = 'glm', family = 'binomial', size = 2) + # plot the fit
#~#~                facet_wrap( ~ Family)
#~        y_var <- ifelse(draw_CpDuplication, "CpDup", "Duplication")
#~        gplot <- ggplot(tab_p, aes_string(x = v_cont[i], y = y_var)) +
#~                facet_wrap( ~ Family) + # split by families
#~                geom_point() + # plot the real data
#~                xlab(v_cont[i]) +
#~                geom_ribbon(data = tab_p, aes(y = fit, ymin = ymin, ymax = ymax), alpha = 0.25) +    # plot CI
#~                geom_line(data = tab_p, aes(y = fit, colour = Family)) # plot the fit

#~        print(gplot)
#~#~        g[[i]] <- gplot
#~        }
#~#~        return(g)
#~}








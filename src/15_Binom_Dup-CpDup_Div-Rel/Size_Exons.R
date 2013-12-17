ind <- seq(1, nrow(GLMtab_all1), by=8)
boxplot(exp(GLMtab_all1$LnExonLength)[ind]~GLMtab_all1$Family[ind])
x11()
boxplot(exp(GLMtab_all1$LnGeneLength)[ind]~GLMtab_all1$Family[ind])

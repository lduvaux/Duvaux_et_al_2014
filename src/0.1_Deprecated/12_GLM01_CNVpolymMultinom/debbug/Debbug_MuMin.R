#### MuMIn example 
library(MuMIn)
data(Cement)
fm1 <- lm(y ~ ., data = Cement)
dd <- dredge(fm1)
summary(model.avg(dd, subset = delta < 4))

# in this case, the 'model.avg' function return the model average
# coefficients for the explanatory variables (as expected)

#### application to multinomial regression
library(foreign)
library(nnet)
library(MASS)
ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
ml$prog2 <- relevel(ml$prog, ref = "academic")
test <- multinom(prog2 ~ ses + write + schtyp, data = ml)
dd <- dredge(test)
summary(get.models(dd, 1)[[1]])
summary(model.avg(dd, subset = delta < 4))

# here, the function returns "model average coefficients" for the different
# levels of the response variable rather than for the model average
# coefficients for the explanatory variables


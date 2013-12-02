#     # load data and libraries
#     library(methods)
#     library(MuMIn)
#     data(Cement)

#     # options
#     MOD_ALL1 <- y ~ X1 + X2 +X3 + X4
#     FIXED_TERMS <- c("X1", "X2")

#     # run analyses
#     print("Fit the main model")
#     fm1 <- lm(formula=MOD_ALL1, data = Cement)  

#     print("Test all terms")
#     test_all_terms <- dredge(fm1, fixed=FIXED_TERMS)
    
#     print("Print results")
#     print(test_all_terms)


    # load data and libraries
    library(methods)
    library(MuMIn)
    library(foreign)
    library(nnet)
#     ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")

#     # options
#     MOD_ALL1 <- prog ~ ses + schtyp + write + science + female + female * schtyp
#     FIXED_TERMS <- c("ses", "write")

#     print("Fit the main model")
#     fm1 <- multinom(formula=MOD_ALL1, data = ml)  

#     print("Test all terms")
#     test_all_terms <- dredge(fm1, fixed=FIXED_TERMS)

#     print("Print results")
#     print(test_all_terms)

    main <- function()
    {
        ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
        MOD_ALL1 <- prog ~ ses + schtyp + write + science + female + female * schtyp
        FIXED_TERMS <- c("ses", "write")

        print("Fit the main model")
        fm1 <- multinom(formula=MOD_ALL1, data = ml)  

        print("Test all terms")
        test_all_terms <- dredge(fm1, fixed=FIXED_TERMS)

        print("Print results")
        print(test_all_terms)
    }
    main()










    

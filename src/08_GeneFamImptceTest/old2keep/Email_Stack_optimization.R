I tested out  on my own data the different solutions proposed by @flodel, @Sven Hohenstein and @Martin Morgan. (***IMPORTANT NOTE: altough all methods give the same result in my specific case, remind thatthey all have their own way, and thus can give different results depending on the question***)

Here is a quick summary (the results are shown below):

 1. In my tests, ```length(a)``` and ```length(b)``` are set to 200 or 400 and 3,500 and 10,500 respectively
 2. there is only a single match of each value of ```a``` in ```b```
 3. ```pmatch``` always performs better than other methods (notably for small length of vectors ```a``` and ```b```, say less than 100 and 1,000 respectively - not shown below), 
 4. ```sapply(a, grep, b, fixed=T)```, ```f3``` (Sven's method) and ```reduced.match``` (flodel's method) functions always perform better than ```sapply(a, grep, b))``` and ```sapply(paste0("^", a), grep, b)```.

However, altough ```pmatch``` is always the quickest in my tests, the best method may really **depends of the problem and all deserves to be tested for each specific cases** (for instance, note that ```f3``` results doesn't change while increasing ```length(a)```):

***Here are the tests***

    # length(a)==200 & length(b)==3500
    system.time(mclapply(seq(100), function(x) sapply(a, grep, b), mc.cores=8))
       user  system elapsed 
     60.308   0.270   9.418 
    system.time(mclapply(seq(100), function(x) sapply(paste0("^", a), grep, b), mc.cores=8))
       user  system elapsed 
     57.325   0.310   7.892 
    system.time(mclapply(seq(100), function(x) sapply(a, grep, b, fixed=T), mc.cores=8))
       user  system elapsed 
     14.506   0.251   1.190 
    system.time(mclapply(seq(100), function(x) reduced.match(a, b), mc.cores=8))
       user  system elapsed 
      4.323   0.248   0.561 
    system.time(mclapply(seq(100), function(x) f3(a, b), mc.cores=8))
       user  system elapsed 
      4.887   0.128   0.699 
    system.time(mclapply(seq(100), function(x) pmatch(a, b), mc.cores=8))
       user  system elapsed 
      2.890   0.056   0.361

    # length(a)==400 & length(b)==3500
    system.time(mclapply(seq(100), function(x) sapply(a, grep, b), mc.cores=8))
       user  system elapsed 
    118.615   0.362  17.957 
    system.time(mclapply(seq(100), function(x) sapply(paste0("^", a), grep, b), mc.cores=8))
       user  system elapsed 
    114.269   0.264  15.435 
    system.time(mclapply(seq(100), function(x) sapply(a, grep, b, fixed=T), mc.cores=8))
       user  system elapsed 
     29.267   0.239   2.337 
    system.time(mclapply(seq(100), function(x) reduced.match(a, b), mc.cores=8))
       user  system elapsed 
      7.686   0.292   0.912 
    system.time(mclapply(seq(100), function(x) f3(a, b), mc.cores=8))
       user  system elapsed 
      5.482   0.121   0.732 
    system.time(mclapply(seq(100), function(x) pmatch(a, b), mc.cores=8))
       user  system elapsed 
      5.001   0.084   0.760 


    # length(a)==200 & length(b)==10500
    system.time(mclapply(seq(100), function(x) sapply(a, grep, b), mc.cores=8))
       user  system elapsed 
    181.080   0.590  29.206 
    system.time(mclapply(seq(100), function(x) sapply(paste0("^", a), grep, b), mc.cores=8))
       user  system elapsed 
    174.396   0.736  24.845 
    system.time(mclapply(seq(100), function(x) sapply(a, grep, b, fixed=T), mc.cores=8))
       user  system elapsed 
     43.470   0.364   3.726 
    system.time(mclapply(seq(100), function(x) reduced.match(a, b), mc.cores=8))
       user  system elapsed 
     11.147   0.275   1.411 
    system.time(mclapply(seq(100), function(x) f3(a, b), mc.cores=8))
       user  system elapsed 
     17.777   0.166   2.745 
    system.time(mclapply(seq(100), function(x) pmatch(a, b), mc.cores=8))
       user  system elapsed 
      8.651   0.091   1.106 

    # length(a)==400 & length(b)==10500
    system.time(mclapply(seq(100), function(x) sapply(a, grep, b), mc.cores=8))
       user  system elapsed 
    356.354   0.898  57.597 
    system.time(mclapply(seq(100), function(x) sapply(paste0("^", a), grep, b), mc.cores=8))
       user  system elapsed 
    343.498   1.033  46.474 
    system.time(mclapply(seq(100), function(x) sapply(a, grep, b, fixed=T), mc.cores=8))
       user  system elapsed 
     86.378   0.276   6.659 
    system.time(mclapply(seq(100), function(x) reduced.match(a, b), mc.cores=8))
       user  system elapsed 
     20.560   0.232   2.256 
    system.time(mclapply(seq(100), function(x) f3(a, b), mc.cores=8))
       user  system elapsed 
     18.651   0.180   2.537 
    system.time(mclapply(seq(100), function(x) pmatch(a, b), mc.cores=8))
       user  system elapsed 
     15.209   0.075   1.999 

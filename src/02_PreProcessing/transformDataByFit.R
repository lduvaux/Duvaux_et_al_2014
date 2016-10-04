fitPolyn2ToData <- function(x,y,use_intercept = T){

    if(use_intercept){
        form <- y ~ a * x^2 + b * x + c
        return (nls(form,start = list(a = 0, b = 1, c = 0)))
    }
    else{
        form <- y ~ a * x^2 + b * x 
        return (nls(form,start = list(a = 0, b = 1)))
    }
}

fitPolyn3ToData <- function(x,y,use_intercept = T){
    if(use_intercept){
        form <- y ~ a * x^3 + b * x^2 + c * x + d
        return (nls(form, start = list(a = 0, b = 1, c = 0,d =0)))   # fit the polynom
    }
    else{
        form <- y ~ a * x^3 + b * x^2 + c * x  
        return (nls(form, start = list(a = 0, b = 1, c = 0)))   # fit the polynom
    }
}


makeModelForY <- function(x, y, use_intercept = T, polynDeg = 1){
	if (polynDeg==1)
		out <- list(mod = NULL,y =y)    
	else if(polynDeg == 2)
		out  <- list(mod = fitPolyn2ToData (x,y,use_intercept),y =y)
	else if(polynDeg == 3)
		out  <- list(mod = fitPolyn3ToData (x,y,use_intercept),y =y)
	else
		stop("wrong degree")
	
	return (out)	
}


transformY <- function(x,model){#, plot_fit = F, main, col.main){
		#just median transfo
		if(is.null(model$mod))
			rat <- median(model$y/x)
			
		else{	
			pred <- predict(model$mod)
			if(sum(is.na(pred)) > 0)
				stop("na in pred")
			
#~ 		if(use_intercept)
#~ 			inters <- coef(model)[length(coef(model))]
#~ 		else 
#~ 			inters <- 0

		# to have the intersept, we ass the model to predict y for x=0 !!
			inters <- predict(model$mod,newdata=list(x=0))
			if( sum(pred== 0) > 0 ){
				print (min(pred))
				print(c)
				print(sum((pred  - inters) == 0))
		#       stop("denomin == 0")
				stop("numerator == 0")
				}
			
			#~ 	rat <- (x + cte - inters ) / (pred + cte - inters)
			#~     rat <- (x + cte) / (pred + cte)	# add by Ludo
			#    rat <- (x + cte - c ) / (pred + cte - c)
			#~     new_y <- y*rat

			rat <- pred / x
    }
    
    new_y <- model$y/rat
    if( sum(is.na(new_y)) > 0)
        stop("na in returned val")
    
    return(new_y)
}


plotYTransfo <- function(new_y,model,x,main,col.main, thresholds){
		samp=sample(1:length(model$y), 2000)
		plot(new_y[samp] ~ x[samp], pch = 20, ylim = c(0,0.02), xlim = c(0,0.02),cex= 0.3, main=main, col.main=col.main, col="red")
#~        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey94")
#~        points(new_y[samp] ~ x[samp], pch = 20, cex= 0.3, col="red")
		points(model$y[samp]~ x[samp], pch = "+", col="green",cex= 0.5)#, main=names(model$y))
		abline( a = 0,b=1,col = 'black', lwd = 1.5, lty = 1)
#~		sapply(c(0.25, 0.75, 1.25, 1.75), abline, a=0, col="green")

        abline(b=1.5,lty=4, a=0, lwd = 1.5)
		abline(b=0.5,lty=4, a=0, lwd = 1.5)
		abline(b=0.75,lty=3, a=0, lwd = 1.5)
		abline(b=1.25,lty=3, a=0, lwd = 1.5)
		abline(b=0.25,lty=3, a=0, lwd = 1.5)
        
		simx <- seq(from = thresholds[1], to = thresholds[2], length.out = 100);
		if (!is.null(model$mod)) lines(predict(model$mod,list(x = simx)) ~ simx,col ='blue', lwd = 2.25, lty = 2)

}

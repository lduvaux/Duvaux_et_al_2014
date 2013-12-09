## PREPROCESSING FUNCTIONS HERE ##
PrePro_quantiles_ctrl <- function(control_vec, thres){
	control_vec <- sqrt(control_vec)
	quant <- quantile(control_vec, thres)
	lower <- names(which(control_vec<quant[1]))
	higher <- names(which(control_vec>quant[2]))
	res <- list(threshold=thres, quantiles=quant, low_ouliers=lower, high_outliers=higher, sqrt_transformation=T)
	return(res)
}

PrePro_showError <- function(err,col_race,return_median_error = F){
    boxplot(err,outline=F,col=col_race,ylim = c(-.25, .25))
    if(return_median_error)
        return(apply(abs(err),2,median))
}


PrePro_findMeaninglessRows <- function(mat){

    invalid_rows <- apply(mat,1,function(row){all(row == row[1])})
    
    print(sprintf("A total of %i rows had the exact same values for all column",length(invalid_rows[invalid_rows])))
    return(invalid_rows)
    
}

PrePro_showError <- function(err,col_race, return_median_error= FALSE, border=1){
    boxplot(err,outline=F,col=col_race,ylim = c(-0.25, +0.25), border=border)
    return(apply(abs(err),2,median))
}



plotAllFits <- function(models,transformed_data_list,control,plot_mains, thresholds, pdfname){
	mains <- plot_mains$mains
	col.mains <- plot_mains$col.mains

	pdf(pdfname,width = 16, height = 9)
		layout(matrix(1:8,2,4,byrow = T))
		
		for(i in 1:length(models)){
			plotYTransfo(transformed_data_list[[i]],models[[i]],x=control,mains[i],col.mains[i], thresholds);
		}
		
	dev.off()	
}
		
		
		
PrePro_prePrecessing <- function(control, data_mat, polyn_deg = 2, use_intercept = TRUE , weight_error = TRUE, thresholds, plot_mains = NULL, doInval_row=F)
# thresholds -> has to be expressed as the same order as the square root of x
{
    control <- sqrt(control)
    fit <- list(use_intercept=use_intercept, polyn_deg=polyn_deg)
    
    val <- (control  > thresholds[1] & control  < thresholds[2])	# keep baits where control is reliable only
    control <- control[val]
    data_mat0 <- data_mat
    data_mat <- sqrt(data_mat)[val,]
   
	#	list of columns(indivs)
	tmp_list_of_cols <- split(data_mat, rep(1:ncol(data_mat), each = nrow(data_mat)))
	names(tmp_list_of_cols) <- colnames(data_mat)
	#	compute models (list of models)
	models <- lapply(tmp_list_of_cols,makeModelForY,x=control,use_intercept = use_intercept, polynDeg = polyn_deg)
	#	transform all elements of the list
	transformed_data_list <- lapply(models,transformY,x = control)
	
	#~ 	if we want to plot the list plot_mains will not be null
	if(!is.null(plot_mains)) {
		if (polyn_deg==1) {
			polyn <- "median"
			intercpt <- ""
		}
		else {
			polyn <- paste("polyn", polyn_deg, sep="")
			if (use_intercept) 
				intercpt <-  "_Intercept" 
			else
				intercpt <- "_NoIntercept"
		}
		pdfname <- paste("Res_FitLm2Data_", polyn, intercpt , ".pdf", sep="")
		plotAllFits(models,transformed_data_list,control,plot_mains, thresholds, pdfname)
	}
    #merging the data back into a  matrix
    data_mat <- do.call("cbind",transformed_data_list)
    rownames(data_mat) <- rownames(data_mat0)[val]	# add name to new sqrt data
    ratio_mat <- apply(data_mat,2,function(x,ctrl,cte){x/(ctrl + cte)},ctrl = control, cte = 0)
    rounded_ratio_mat <- apply(ratio_mat,2,PrePro_roundToZeroFive)
    
    
    if (doInval_row) 
		inval_rows <- PrePro_findMeaninglessRows(rounded_ratio_mat) 
    
    else 
		inval_rows=NULL	# remove lines where rounded values are equals for all individuals
   
	# matrix of errors
    error_mat <- ratio_mat - rounded_ratio_mat
    weight <- ifelse(0,1,rounded_ratio_mat == 1)
    error_mat_weighted <- weight * error_mat

#	preparation to return data
#~     rownames(rounded_ratio_mat) <- 1:nrow(rounded_ratio_mat)
    
    if(weight_error)
		out_list <- list(raw_y=data_mat0, new_sqrt_y=data_mat, new_ratio = ratio_mat, error_round = error_mat_weighted, good_baits=val, noInfo=inval_rows, fit=fit)
		
	else
		out_list <- list(raw_y=data_mat0, new_sqrt_y=data_mat, new_ratio = ratio_mat, error_round = error_mat, good_baits=val, noInfo=inval_rows, fit=fit)
		
	if (doInval_row)	# if we have removed colums with no variation
	{
		out_list$raw_y <- data_mat0[!inval_rows,]
		out_list$new_sqrt_y <- out_list$new_sqrt_y[!inval_rows,]
		out_list$new_ratio <- out_list$new_ratio[!inval_rows,]
		out_list$error_round <- out_list$error_round[!inval_rows,]
		out_list$good_baits <- out_list$good_baits[!inval_rows]
		out_list$fit <- fit
    }
    return(out_list)
}

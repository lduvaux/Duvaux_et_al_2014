#!/bin/Rscript
source("../utils/functions.R")
source("../utils/globalCtes.R")
source("params.R")
source("./functions.R")

main <- function(argv){

    load(PREVIOUS_DATA)
    
    raw_data_file <- RAW_DATA
    load(raw_data_file)
    data_mat <- PrePro_fetchRawDataMat(raw_data_file, BAD_GUYS); dbgPrint(dim(data_mat))
    control_vec <- PrePro_fetchRawCtrl(raw_data_file)
    
    lib_name <- PrePro_LibName(indiv_details=INDIV_DETAILS,indiv=colnames(data_mat))
    col_lib <- PrePro_LibColo(lib_name, col_lib_raw)
    
    ratio_mat <- apply(data_mat,2,function(x,ctrl,cte){x/(ctrl + cte)},ctrl = control_vec, cte = 0)
	max <- 6
	pdf(PDF_NAME, width=12,height=8)
	par(layout(matrix(1:12,4,3)))
	par(mar=c(1,1,1,1))
	for(i in 1:ncol(ratio_mat)){
		xlab= expression("Depth"["target"]/"Depth"["control"])
# 		plot(1,1, xlab = xlab,ylab="Density",ylim=c(0,10),xlim=c(0.4,max/4),pch="",axes=F)
		plot(1,1, xlab = xlab,ylab="Density",ylim=c(0,10),xlim=c(0,max/4),pch="",axes=F)
#~ 		axis(1,tick = F,label=F,at=(0:8)/4)
		axis(1,tick = F,label=F)

		box()
#~ 		plot(1,1, main = colnames(ratio_mat)[i],ylim=c(0,10),xlim=c(0.4,max/4),pch="",axes=F)
        abline(v=0)
		abline(v=1,lty=2)
		
		abline(v=1.5,lty=4)
		abline(v=0.5,lty=4)
		
		abline(v=0.75,lty=3)
		abline(v=1.25,lty=3)
		abline(v=0.25,lty=3)
		
		d1 <- density(data_list$new_ratio[data_list$new_ratio[,i]<max,i],bw=.02)
		d2 <- density(na.omit(ratio_mat[ratio_mat[,i]<max,i]),bw=.02)
		lines(d1$y ~ d1$x,lwd=2,col="red")
		lines(d2$y ~ d2$x,lwd=2,col="green")
		
		abline(v = median(na.omit(ratio_mat[,i])),col="green")
		abline(v = median(data_list$new_ratio[,i]),col="red")
		
#~ 		text("topleft",)
		legend("topleft", colnames(ratio_mat)[i],cex=1.5,col="blue",bg="white", text.col=col_lib[i]) 
		}
	dev.off()

    #############################
    cat("\n")
    print(" #### save results")
	outFileName <- argv[1]
    ver(sprintf("Saving *DUMMY* data to %s",outFileName))
    dummy <- numeric()
    save (dummy, file=outFileName)
}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);

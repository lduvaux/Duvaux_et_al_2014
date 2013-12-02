#!/bin/Rscript
source("../utils/functions.R")
source("../utils/globalCtes.R")
source("params.R")
source("./ProcessingSegmentationFunctions.R")

main <- function(argv){

    load(PREVIOUS_DATA)
    
    raw_data_file <- RAW_DATA
    print(names(data_list))    
    
    fit <- data_list$fit
    print(fit)
    
    
    # load bait info
    bait_info <- PrePro_fetchBaitInfo(BAIT_GOOD_Q10RATIO); dbgPrint(dim(bait_info))
    bait_categ <- PrePro_BaitCateg(bait_info$name); dbgPrint(length(bait_categ))	# used PrePro_BaitCateg2 instead (use mclapply)?
    
    good <- data_list$good_baits
    control_vec <- sqrt(PrePro_fetchRawCtrl(raw_data_file)[good])

    bait_info <- bait_info[good,]
    bait_categ <- bait_categ[good]
        
    new_sqrt_y <- data_list$new_sqrt_y
    

    contigs_uniq <- unique(bait_info$chrom)
    resSegmentation <- pro_runPerContig(contigs_uniq, contigs=bait_info$chrom, new_sqrt_y, control_vec, nind=ncol(new_sqrt_y))

    alpha_matrix <- resSegmentation$AlphaMatrix
    
    good_baits <- do.call('cbind',strsplit(rownames(alpha_matrix),"_"))
    df <- as.data.frame(t(good_baits)[,1:4])
    colnames(df) <- c("GeneCateg","Gene","Exon","Bait")
    gene_per_categ <- aggregate(Gene ~ GeneCateg, df,function(x)length(unique(x)))

    #############################
    outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    save (fit,alpha_matrix,gene_per_categ, file=outFileName)
}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);

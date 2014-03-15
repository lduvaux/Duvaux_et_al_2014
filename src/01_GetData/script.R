#!/bin/Rscript

# script to detect CNV in my targets
# setup inputs and run PropSeq

source("../utils/functions.R")
source("../utils/globalCtes.R")
source("params.R")


main <- function(argv){
    
	#############END OF TO_CLEANUP###########
		
		
	################### 0) preliminary set up
	# 0.0) list all files to be processed
	targ_fils <- dir(DFIL, TARG_PATT, full.names=T)
	metrics_fils <- dir(DFIL, METRICS_PATT, full.names=T)

    set.seed(0)

	# 0.1) Normalization1: N=C/R where C=mean coverage per bp per bait, R total nuber of aligned read for this individual
		# 0.1.0) define the number of good targets (less than 5% of reads with Q<10)
	tab_good <- read.delim(FIL_GOOD, stringsAsFactors =F)
	temp_tab <- read.delim(targ_fils[1])
	good <- temp_tab$name%in%tab_good$name	# targets to be kept after reads cleaning
	temp_tab <- temp_tab[good,]
	identical (as.character(tab_good$name), as.character(temp_tab$name))
	ntar <- nrow(temp_tab); rm(temp_tab)

		# 0.1.1) set up the matrix
	print("\n############ prepare the matrix of of normalized counts for all individuals #############\n")
# 	tested_nbreads <- numeric(length(metrics_fils))
	tested_tab <- matrix(data=NA, nrow=ntar, ncol=length(metrics_fils))
	for (ite0 in 1:length(metrics_fils))
	{
		# total number of reads
		fil <- metrics_fils[ite0]
		nbread <- read.delim(fil, comment.char = "", skip=6)[,"TOTAL_READS"]
# 		tested_nbreads[ite0] <- nbread
		
		# mean coverage per target
		fil2 <- targ_fils[grep(fil,targ_fils)]
		tab <- read.delim(fil2); tab <- tab[good,]
		tested_tab[,ite0] <- tab[,"mean_coverage"]/nbread
		if (round(ite0%%(length(metrics_fils)/100),0)==0) print(paste(round((ite0*100)/length(metrics_fils),2), "% done", sep=""))
	}
	print("\n############ matrix of tested individuals: done #############\n")

	# 0.2) fetch sequencing information
		# 0.2.1) information from metric files
	ID <- as.character(sapply(metrics_fils, function(x) unlist(strsplit(x, "//|_R1s1a1"))[2]))
	nbtested <- as.character(sapply(ID, function(x) unlist(strsplit(x, "_"))[2]))
# 	nb.reads.nodup <- tested_nbreads
	colnames(tested_tab) <- ID
	vrace <- sapply(ID, function(x) unlist(strsplit(x, "_"))[1])

		# 0.2.2) sort the two spec files to have them in the same order as the metrics files (as thus as the columns order in the matrix)
	tspec_fil <- read.delim(SPEC_FIL)[seq(1, 241, by=2),1:2]
	baba0 <- grep("AD067GACXX", tspec_fil$fastq.file);tspec_fil=tspec_fil[-baba0,]	# remove file "en trop"
	indiv <- as.character(tspec_fil[,1])
	indi <- sapply(nbtested, function(x) which(indiv==x))	# index to have the same order as in metric file list
	tspec_fil <- tspec_fil[indi,]
		# 0.2.3) set up last fields
	run_date <- as.numeric(sapply(as.character(tspec_fil[,2]), function (x) unlist(strsplit(x, "/|_"))[3]))
	run_nb <- as.numeric(sapply(as.character(tspec_fil[,2]), function (x) unlist(strsplit(x, "/|_"))[4]))
	mach_run_nb <- as.character(sapply(as.character(tspec_fil[,2]), function (x) unlist(strsplit(x, "/|_"))[5]))
	lane <- as.numeric(sapply(as.character(tspec_fil[,2]), function (x) unlist(strsplit(x, "/|_"))[6]))
	vrunID <- paste(run_date, lane, mach_run_nb, sep="_")
	runID <- sort(unique(vrunID))

		# 0.2.4)  set up ref vectors (median of the REF individuals) and graph parameters
	x <- sapply(1:nrow(tested_tab), function(x) median(tested_tab[x,CTRL_GUYS]))	# vector of ref data

		# 0.2.5) save data
	tested_tab_raw <- tested_tab
	good_targets <- sapply(tab_good[,"name"], Prepro_fixName)
	ctrl_vec <- x
	ctrl_vec_raw <- ctrl_vec
	
	
	#############END OF TO_CLEANUP###########
 
    ######################
    outFileName <- argv[1]
    ver(sprintf("Saving data to %s",outFileName))
    save(tested_tab_raw, good_targets, ctrl_vec_raw, vrace, vrunID, runID, file=outFileName);
}

argv <- commandArgs(TRUE)[1]
if(DEBUG)
	traceback(main(argv));
if(!DEBUG)
	main(argv);
	














source("../10_Functions_GLMMs/functions.R")
source("../11_GLMTruncFetchData/glm_common.R")

###############
get_SeqDepth <- function(fil_seq_depth)
{

	tab <- read.delim(fil_seq_depth, stringsAsFactors = F)
	bad <- which("Target_Name"==colnames(tab))
    rownames(tab) <- tab[,bad]
	tab <- tab[,-bad]
	tab <- as.matrix(tab)
	return(tab)
}

###############
get_readsNber <- function(fil_fastqc, bad_file)
{
	tab <- read.delim(fil_fastqc, stringsAsFactors = F)
	bad <- grep(bad_file, tab[,"RunInd"])
	tab <- tab[-bad,]
	return(tab)
}

###############
set_ftab <- function (tab_reads, badFcell, badIndiv, Met_indiv)
{
	# 1) define date, flowcell and lane
	SqcingDate <- sapply(tab_reads$RunInd, ext_runID, 1)
	Flowcell <- sapply(tab_reads$RunInd, ext_runID, 3)
	Lane <- sapply(tab_reads$RunInd, ext_runID, 4)
	GenepoolTag <- sapply(tab_reads$RunInd, ext_runID, 5)

	# 2) find correspondance between individuals
	ind <- sapply(tab_reads$dir, new_grep, names(Met_indiv))

	# 3) final table
		# 3.1) set up
	clones <- names(Met_indiv)[ind]
	tab <- data.frame(clones, dir=tab_reads[,"dir"], NberOfReads=tab_reads[,"ReadNber"], ReadLength=tab_reads[,"length"], MedSeqDepth=Met_indiv[ind], Encode=tab_reads[,"Encode"], file1=tab_reads[,"fil1"], file2=tab_reads[,"fil2"], SqcingDate, Flowcell, Lane, GenepoolTag)

		# 3.2) order
	ind2 <- order(SqcingDate, Flowcell, Lane, GenepoolTag)
	tab <- tab[ind2,]
		
		# 3.3) add vector of use for CNV
	CNVestimation <- rep(T, nrow(tab))
	bad1 <- which(tab[,"Flowcell"]==badFcell)
	bad2 <- sapply(badIndiv, function(x) which(x==tab[,"dir"]))
	CNVestimation[c(bad1, bad2)]=F
	tab <- data.frame(tab, CNVestimation)
	return(tab)
}

ext_runID <- function(runID, field)
{
	res <- unlist(strsplit(runID, "_"))[field]
	return(res)
}

new_grep <- function(pattern, x)
{
	res <- grep(paste("_", pattern, "_", sep=""), x)
	return(res)
}

######################
get_insert_size <- function(fil, rpint=F){
    if (rpint) print(fil)
    con <- file(fil, open="r")
    txt <- readLines(con)
    close(con)
    gd_vec <- txt[8]
    ins_size <- as.numeric(unlist(strsplit(gd_vec, "\t"))[[1]])
    return(ins_size)
}

batch_ins_size <- function(path, pattern){
    files <- dir(path, pattern)
    sample_name <- sapply(files, collapse_elements)
    fils <- dir(path, pattern, full.names = T)
    insert_size <- as.numeric(sapply(fils, get_insert_size))
    tab <- data.frame(files, sample_name, insert_size)
    return(tab)
}

######################
get_capture_metrics <- function(fil){
    tab <- read.delim(fil, skip=6, stringsAsFactors = F)
    enrichment <- round(tab [, "FOLD_ENRICHMENT"],2)
    efficiency <- round(tab [, "ON_BAIT_BASES"]/tab [, "PF_UQ_BASES_ALIGNED"],4)
    PTB_30 <- round(tab [, "PCT_TARGET_BASES_30X"], 3)
    res <- c(enrichment, efficiency, PTB_30)
    return(res)
}


batch_capture_metrics <- function(path, pattern){
    files <- dir(path, pattern)
    sample_name <- sapply(files, collapse_elements)
    fils <- dir(path, pattern, full.names = T)
    tab <- t(sapply(fils, get_capture_metrics))
    colnames(tab) <- c("Enrichment", "Efficiency", "PTB30X")
    indiv <- sapply(fils, function(x) unlist(strsplit(x, "\\/|_R1s1a1"))[5])
    rownames (tab) <- indiv
    return(tab)
}






source("./warnings2logfile.R")
source("../utils/getter_functions.R")
source("../utils/globalCtes.R")

library(parallel)
library(optimalCaptureSegmentation)

# setup nk for 'findOptimalSegmentations' (nk= the maximal nb of segments k
	# having different Copy Number in a given chromosome)
Pro_setup_nk <- function(nkmax, long)
{
	temp <- ifelse(
		long<=nkmax, long, ifelse(
			long>=(nkmax*4), floor(long/4), ifelse(
				long<=(nkmax*2), nkmax, floor(long/2)
			)
		)
	)
	return(temp)
}

# change modifiedBIC=Inf to -inf in 'findOptimalSegmentations' results
Pro_changeInf <- function(modifiedBIC, solutions, tagInf)
{
	if (sum(modifiedBIC==Inf, na.rm=T)>0) 
		for(tag in tagInf) solutions[[tag]]$modifiedBIC=-Inf
	return(solutions)
}

# basic function to run the 'findOptimalSegmentations':
	# perform a segmentation analysis across nk values of k
	# check results for NAs and Inf in modifiedBIC
	# create log files
# remark concerning optimalCaptureSegmentation: In some cases, the function can produce warnings if a segmentation solution produces “sumErrorSquared” equal to 0 or very low (e.g. -1.552981e-25) - usually in cases where k is high. The first case produces an estimation of “sigmaEstimated” equal to 0 and a BIC value equal to “inf”. In order to discard these patterns, I changed the respective value of BIC value from “inf” to “-inf”. The second case give NaN for both values. I have discarded these patterns as well.
Pro_runSegmentation <- function(y, x, nk, contig, individual, fil0, fil="UnsuitableIndByContigs.log", fil2="NoMaxDetected.log", fil3="Best-k_is_nk.log", fil4="Zero_in_y.log", allLog_in_one=T, pref)
{
	if (allLog_in_one) {fil=fil0; fil2=fil0; fil3=fil0; fil4=fil0}
	# detect modified BIC with NA values
		# setup data
	long=length(x)
	warning_zero_in_y(y, contig, individual, long, fil0)
		# run findOptimalSegmentations
	solutions <- findOptimalSegmentations(y, x, nk)
	modifiedBIC <- sapply(solutions, function(sol) sol$modifiedBIC)
	
	# detect NAs in solutions
	tagNA <- Pro_setTagNA(modifiedBIC=modifiedBIC, contig=contig, individual=individual, nk=nk, pref="First solutions", fil=fil0, long=long)
	good <- rep(T, length(modifiedBIC))
	if (!is.null(tagNA)) good[tagNA]=F
	warning_only_NA(good, contig, individual, pref, fil0)
	
	# change Inf for -Inf values
	tagInf <- Pro_setTagInf(good=good, modifiedBIC=modifiedBIC, contig=contig, individual=individual, nk=nk, pref=pref, fil=fil0, long=long)
	solutions <- Pro_changeInf(modifiedBIC, solutions, tagInf)
	modifiedBIC <- sapply(solutions, function(sol) sol$modifiedBIC)
	
	# warnings to log files
	warning_suitable_BIC(tagNA, tagInf, modifiedBIC, contig, individual, pref, fil0)
	warning_noMax_BIC(modifiedBIC, contig, individual, pref, fil0)
	
	return(list(solutions=solutions, modifiedBIC=modifiedBIC, tagNA=tagNA, tagInf=tagInf, Nbaits=long))
}

# where needed, run again 'Pro_runSegmentation' by increasing nk
Pro_runSeg_nkExten <- function(y, x, nk, contig, individual, logf)
{
	AllSeg <- Pro_runSegmentation(y, x, nk, contig, individual, pref="First solutions", fil0=logf)
	long <- AllSeg$Nbaits
	# if best n == nkmax, rerun with a higher nkmax till best n != nkmax
	if (length(which.max(AllSeg$modifiedBIC))>0)
	{
		ite=1
		while (which.max(AllSeg$modifiedBIC)==nk & nk<(long-1))
		{
			if (ite==1) 
				warning_Bestk_is_nk(contig, individual, ite, pref="First solutions", fil=logf)
			else 
				warning_Bestk_is_nk(contig, individual, ite, pref=paste("Loop", ite-1, sep=" "), fil=logf)
			
			if (floor(nk*1.5)<long) nk=floor(nk*1.5) else nk=long-1
			AllSeg <- Pro_runSegmentation(y, x, nk, contig, individual, pref=paste("Loop", ite, sep=" "), fil0=logf)
			
			ite=ite+1
		}
	}	
	return(AllSeg)
}

# from 'Pro_runSeg_nkExten' results, select the best segmentation for my 
	# purpose (here the solution with the highest modifiedBIC)
Pro_select_BestSeg <- function(AllSeg, y)
{
	warning("The algorithm of 'Pro_select_BestSeg' is based on the observation that the 
	following patterns are always associated together:
	  - only NaN in modifiedBIC
	  - No max observed in modifiedBIC
	  - No suitable values in modifiedBIC
	  - k=1 (i.e. I've never found k=1 independently of the three above observations).
	They are observed in either of these two cases
	  i) there is only one bait on the chromosome under study
	  ii) the y vector contain only 0 values (i.e. depth of coverage=0 for all baits 
			of the chromosome under study).")
	solutions <- AllSeg$solutions
	vBIC <- AllSeg$modifiedBIC
	# For (i) and (ii), we can apply:
	if (length(solutions)==1|(length(unique(y))==1 & y[1]==0)| (length(unique(vBIC))==1 & is.na(vBIC[1])))
		best <- 1
	else
		best <- which.max(vBIC)
		# for all other cases, we apply:
	return(best)
	
}

# attribute an alpha value to each bait based on the best segementation
	# NB: segmentation is displayed with upper limit of each segment included, 
	# e.g: example$segmenation 
		# 6 17 
	# means that the first segment lies from bait 1 to bait 6, and segment 2 from bait 7 to bait 17.
Pro_alphaPerBait <- function(AllSeg, best)
{
	warning("Bear in mind that segmentation is displayed with upper limit
	of each segment included, e.g.: 
	example$segmenation 
	  6 17
	means that the first segment lies from bait 1 to bait 6, and segment 2 from bait 7 to bait 17")
	k <- AllSeg$solutions[[best]]$k
	seg <- AllSeg$solutions[[best]]$segmentation
	alphas_raw <- AllSeg$solutions[[best]]$alphaEstimated
	
	if (k==1) 
		alphasBestSeg <- rep(alphas_raw, seg)
	else {
		vrep=c(seg[1], sapply(2:length(seg), function(x) seg[x]-seg[x-1])) 	# setup vector of repetition
		alphasBestSeg <- rep(alphas_raw, vrep)
	}
	return(alphasBestSeg)
}

# perform a segmentation analysis of one individual of a given contig
Pro_used_in_lapply <- function(individual, tab, x, nk, contig, logf)
{
    if (DEBUG) print(individual)
    if (is.matrix(tab)) y <- tab[,individual] else y <-  tab[individual]
	AllSeg <- Pro_runSeg_nkExten(y, x, nk, contig, individual, logf=logf)
	best <- Pro_select_BestSeg(AllSeg, y)
	alphasBestSeg <- Pro_alphaPerBait(AllSeg, best)
	
	return(list(ResAllSeg=AllSeg, alphasBestSeg=alphasBestSeg))
}

# perform a segmentation analysis for all individuals of a given contig
Pro_runSegPerContig <- function(ite0, contigs_uniq, contigs=bait_info$chrom, new_sqrt_y, control_vec, nkmax=12, logf="OptimalSegmentations.log")
{
	# 1) prepare data to run the 'findOptimalSegmentations' function
	contig=contigs_uniq[ite0]
	ind=which(contigs==contig)
	tab=new_sqrt_y[ind,]
	print(paste("Contig #", ite0,": contig", contig, sep=" "))
	
	# 2) run the function
	if (length(ind)>1) {
		nk=Pro_setup_nk(nkmax, nrow(tab))
		individual <- colnames(tab)}
	else {
		warning_1bait_only(contig, logf)
		nk=1 ; individual <- names(tab)
	}
#	list_seg_solut <- mclapply(individual, Pro_used_in_lapply, tab, control_vec[ind], nk, contig)
	list_seg_solut <- lapply(individual, Pro_used_in_lapply, tab=tab, x=control_vec[ind], nk=nk, contig=contig, logf=logf)
	printprogress(ite0, length(contigs_uniq))
	names(list_seg_solut) <- individual
	return(list_seg_solut)
}

# setup the matrix of alphas from the segmentation results of all contigs
Pro_setup_AlphaMatrix <- function(all_list_seg_solut, individual, nind)
{
	AlphaMatrix <- lapply(individual, function(x) sapply(1:length(all_list_seg_solut), function(y) all_list_seg_solut[[y]][[x]]$alphasBestSeg))
	vec <- unlist(AlphaMatrix)
	nro <- length(vec)/nind
	AlphaMatrix <- matrix(vec, nro, nind)
	colnames(AlphaMatrix) <- individual
	return(AlphaMatrix)
}

# final function to be used in 'main_prosegment.R': 
	# perform a thorough segmentation analysis for all contigs
pro_runPerContig <- function(contigs_uniq, contigs=bait_info$chrom, new_sqrt_y, control_vec, suf_log="^.*log$", nind)
{

	individual <- colnames(new_sqrt_y)

	if (length(dir(pattern=suf_log))!=0) file.remove(dir(pattern=suf_log))
	all_list_seg_solut <- mclapply(1:length(contigs_uniq), Pro_runSegPerContig, contigs_uniq, contigs, new_sqrt_y, control_vec)
#~ 	all_list_seg_solut <- lapply(1:length(contigs_uniq), Pro_runSegPerContig, contigs_uniq, contigs, new_sqrt_y, control_vec)
#	all_list_seg_solut <- lapply(1:length(contigs_uniq), Pro_runSegPerContig, contigs_uniq, contigs, new_sqrt_y, control_vec)	# best for debbugging purposes
	names(all_list_seg_solut) <- contigs_uniq#[280:285]
	AlphaMatrix <- Pro_setup_AlphaMatrix(all_list_seg_solut, individual, nind)
	rownames(AlphaMatrix) <- rownames(new_sqrt_y)
	return(list(List_RawRes=all_list_seg_solut, AlphaMatrix=AlphaMatrix))
}


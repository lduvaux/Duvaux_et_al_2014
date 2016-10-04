warning_1bait_only <- function(contig, fil)
{
	cat(paste("###### Contig:", contig, " | Only 1 bait available ######\n", sep=""), file=fil, append=T)
}

warning_zero_in_y <- function(y, contig, individual, long, fil)
{
	if (length(unique(y))==1 & y[1]==0)
		cat(paste("Contig:", contig, " ; Individual:", individual, " | Only 0 in y ; long=", long, "\n", sep=""), file=fil, append=T)
	else
		if (0%in%y) 
			cat(paste("Contig:", contig, " ; Individual:", individual, " | Y contains 0 values\n", sep=""), file=fil, append=T)
}

warning_only_NA <- function(good, contig, individual, pref, fil)
{
	if (sum(good)==0) cat(paste("Contig:", contig, " ; Individual:", individual, " | ", pref, " | Only NaN in modifiedBIC\n", sep=""), file=fil, append=T)
}

warning_suitable_BIC <- function(tagNA, tagInf, modifiedBIC, contig, individual, pref, fil)
{
	if (length(tagNA)+length(tagInf)==length(modifiedBIC)) cat(paste("Contig:", contig, " ; Individual:", individual, " | ", pref, " | No suitable modifiedBIC values\n", sep=""), file=fil, append=T)
}

warning_noMax_BIC <- function(modifiedBIC, contig, individual, pref, fil)
{
	if (length(which.max(modifiedBIC))==0) cat(paste("Contig:", contig, " ; Individual:", individual, " | ", pref, " | No max modifiedBIC detected. Only NaN in modifiedBIC?\n", sep=""), file=fil, append=T)
}

warning_Bestk_is_nk <- function(contig, individual, ite, pref, fil)
{
	cat(paste("Contig:", contig, " ; Individual:", individual, " | ", pref, " | The best k is nk -> process to nk extension\n", sep=""), file=fil, append=T)
}

Pro_setTagNA <- function(modifiedBIC, contig, individual, nk, long, fil="Warnings_Pro_setTagNA.txt", pref="")
{
	if (sum(is.na(modifiedBIC))>0) {
		k <- which(is.na(modifiedBIC))
		cat(paste("Contig:", contig, " ; Individual:", individual, " | ", pref, " | NAs in modifiedBIC ; nk=", nk, " ; k=", k, " ; long=", long, "\n", sep=""), file=fil, append=T)
		tagNA=which(is.na(modifiedBIC))
	} else 
		tagNA=NULL
	return(tagNA)
}

Pro_setTagInf <- function(good, modifiedBIC, contig, individual, nk, long, fil="Warnings_Pro_setTagInf.txt", pref="")
{
	if (sum(modifiedBIC[good]==Inf)>0) {
		k <- which(modifiedBIC==Inf)
		cat(paste("Contig:", contig, " ; Individual:", individual, " | ", pref, " | modifiedBIC contains Infs ; nk=", nk, " ; k=", k, " ; long=", long, "\n", sep=""), file=fil, append=T)
		tagInf=which(modifiedBIC==Inf)
	} else 
		tagInf=NULL
}

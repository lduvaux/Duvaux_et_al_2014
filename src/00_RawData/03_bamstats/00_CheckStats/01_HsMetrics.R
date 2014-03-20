####### I) global metrics
rm(list=ls())
setwd("/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/")

fils=dir(pattern="HsMetrics$")

tab=matrix(nrow=length(fils), ncol=40, data=NA)
for (ite0 in 1:length(fils))
{
	ite=fils[ite0]
	IN<-file(ite,open="r")
	tfil<-readLines(IN) #on slurpe le fichier de sortie de SITES
	close(IN)
	
	if (ite0==1) colnames(tab)=unlist(strsplit(tfil[7], "\t"))
	vstat=unlist(strsplit(tfil[8], "\t"))
	tab[ite0, 1:length(vstat)]=vstat
	print(ite)
}

tab=cbind(file=fils, tab)
tab[1:2,]
tab[,1:3]

subtab=tab[,c(1,7, 12:17,19, 21:24, 26:32, 37:38)]
write.table(tab, file="/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/00_CheckStats/HsMetrics.txt", quote=F, sep="\t", row.names=F)

pdf("/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/00_CheckStats/HsMetrics.pdf")
layout(matrix(1:4, nrow=2, ncol=2, byrow = TRUE))
for (ite in 2:ncol(subtab))
{
	vec=as.numeric(subtab[,ite])
	lims=range(vec)
	l1=floor(lims[1]);l2=ceiling(lims[2])
	man=paste("range: ", lims[1], "-", lims[2], sep="") 
#	hist(vec, xlab=colnames(subtab)[ite], xlim=c(l1, l2))
	hist(vec, xlab=colnames(subtab)[ite], main=man)
}
dev.off()


######### II) metrics for targets
	# clean and load library
rm(list=ls())
library(graphics)
library(fields)
	# load table of metrics per individual (no target details)
tab=read.delim("/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/00_CheckStats/HsMetrics.txt", stringsAsFactors = F)

	# change directory and list files
setwd("/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/")
fils=dir(pattern="Targets$")
identical(tab[,1], gsub(".Targets", "", fils))	# check if the order is the same in tab and in the listing of fils
ord=order(sapply(1:nrow(tab), function(x) tab[x,7]))	# vector of individuals (from files) ordered by total sequencing depth (coverage)


# II.1) set up the matrix for 2D grid plot (matrix of mean coverage for individual i at target j)
ttab=as.data.frame(matrix(nrow=357840, ncol=9, data=NA))	# table of metrics per targets per individuals (direct concatenation of all "*_NoDup.HsMetrics.Targets" files)
mato0=matrix(nrow=120, ncol=2982, data=NA)	# matrix of mean coverage for individual i at target j ; targets ordered like in the original list file ; files ordered alpahabetically
inc=1
for (ite0 in 1:length(fils))
{
	ite=fils[ite0]
	tt=read.delim(ite)
	nom=as.character(tt[,"name"])
	if (ite0==1) nom0=nom
	#print(identical(nom0, nom))	# test if target names are always in the same order (so if tables are comparable) ==> yes, always true
	tt=cbind(file=ite, tt)
	if (ite0==1) colnames(ttab)=colnames(tt)
#	print(dim(tt))
	ttab[inc:(inc+nrow(tt)-1),]=tt
	mato0[ite0,]=tt[,8]
	inc=inc+nrow(tt)
	if (ite0%%10==0) print(round(ite0/length(fils)*100, 2))
}


# II.2) distribution of the mean sequencing depth of targets across all individuals
pdf("/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/00_CheckStats/HsMetrics.targets.pdf")
layout(matrix(1:4,2,2, byrow = TRUE))
hist(ttab[,"mean_coverage"])
hist(subset(ttab, ttab[,"mean_coverage"]<=2000)[,"mean_coverage"])
rm(ttab)
dev.off()


# II.3) estimate the number of targets with at least N individuals with a sequencing depth (coverage) of at least coveX for all race
vN=c(10,11,12)
matf=t(mato0)	# matrix of mean coverage for target i for individual j ; targets ordered like in the original list file ; files ordered alpahabetically
colnames(matf)=gsub("_R1s1a1_NoDup.HsMetrics.Targets", "", fils)
vrace=unique(sapply(colnames(matf), function(x) unlist(strsplit(x, "_"))[1]))
nrace=table(sapply(colnames(matf), function(x) unlist(strsplit(x, "_"))[1]))
vcove=c(10, 20, 30, 50)
vgtarg=matrix(nrow=length(vN), ncol=length(vcove), data=NA, dimnames=list(paste("N", vN, sep=""), paste(vcove, "X",sep="")))
for (N0 in 1:length(vN))
{
	N=vN[N0]
	print(paste("number of individuals:", N))
	for (cove0 in 1:length(vcove))
	{
		cove=vcove[cove0]
		mtest=sapply(1:nrow(matf), function(x) ifelse(matf[x,]<cove, 0, 1))	# matrix of coverage test result for individuals i and target j
		mat.targ.count=matrix(nrow=nrow(matf), ncol=length(vrace), data=NA)	# nber of individuals for target i and race j with at least coveX coverage
		for(race0 in 1:length(vrace))
		{
			race=vrace[race0]
			mtest2=mtest[grep(race, rownames(mtest)),]
			targ.count=sapply(1:ncol(mtest2), function(x) sum(mtest2[,x]))	# vector of the number of individuals for target i for race race0
			mat.targ.count[, race0]=targ.count
			# print(dim(mtest2))	# number of individuals per race ok
		}
		rownames(mat.targ.count)=nom
		colnames(mat.targ.count)=vrace
		# vector of the nb of targets with at least 10 individuals per race with coverage X
		vgood=sapply(1:nrow(mat.targ.count), function(x) ifelse(length(which(mat.targ.count[x,]<N))==0, T, F))
		vgtarg[N0, cove0]=length(vgood[vgood])
		print(length(vgood[vgood]))
	}

}
	# write matf table in a file
matf=matf[,ord]	# mean coverage for target i for individual j ordered by increasing sequencing coverage 
matf=cbind(Target_Name=as.character(tt[,6]), matf)
write.table(matf, file="/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/00_CheckStats/CoveragePerTargetPerInd.txt", quote=F, sep="\t", row.names=F)


# II.4) 2D colored grid for coverage per individual per target (matrix of mean coverage for individual i at target j)
# see "/windows/Documents and Settings/Ludovic/Mes documents/01_professionnel/06_these/export_cca/20_rédaction/papier1/soumis/06_ABC_biais/10_check_Tm/01_10.000simuls/06_Tm/CV_distrib.R"

colnames(mato0)=1:ncol(mato0)	# mean coverage for individual i for target j
rownames(mato0)=1:nrow(mato0)
vthres=c(50, 100, 200, 500)
for (thres in vthres)
{
	mato=mato0
	mato[which(mato>thres)]=thres	# fix the upper limit in order to display what is relevant only 
	mato=mato[ord,]	# vector of individuals (from files) ordered by total sequencing depth (coverage)
	pdfnom=paste("/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/00_CheckStats/HsMetrics.targets.grid-Thres", thres,".pdf", sep="")
	pdf(pdfnom, width=16, height=16)

	cocos=c("darkgreen", "yellow", "orange","red", "darkred")
	jet.colors <- colorRampPalette(cocos[length(cocos):1])
	nbcol=10
	color=jet.colors(nbcol)
	maxi=1
	brk=round(seq(0,maxi, length.out=(nbcol+1)),2)

	# draw
	image.plot(mato, col=color, ylab="y", xlab="x", xaxt="n")
	title("test", line=1)
	axis(1, at=seq(from=0, to=1,length=12), labels=seq(from=10,to=120, length=12))
	dev.off()
	print(mean(mato))
}








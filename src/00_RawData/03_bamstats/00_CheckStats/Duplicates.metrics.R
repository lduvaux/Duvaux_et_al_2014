rm(list=ls())
	# change directory and list files
setwd("/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/")
fils=dir(pattern="_NoDup.DupMetrics$")

tab=matrix(nrow=length(fils), ncol=9, data=NA)
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

range(as.numeric(tab[,"PERCENT_DUPLICATION"]))

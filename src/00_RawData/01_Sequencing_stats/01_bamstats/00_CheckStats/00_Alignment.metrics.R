rm(list=ls())
	# change directory and list files
setwd("/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/")
fils=dir(pattern="_NoDup.AliSumMetrics$")

tab=matrix(nrow=length(fils), ncol=25, data=NA)
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

setwd("/data/bo1ld/06_NERC_capture/02_TrueStampyRuns/01_run1/03_bamstats/00_CheckStats")
mean(as.numeric(tab[,"PCT_PF_READS_ALIGNED"]))
range(as.numeric(tab[,"PCT_PF_READS_ALIGNED"]))
hist(as.numeric(tab[,"PCT_PF_READS_ALIGNED"]))
dev.off()

# 1) looking % aligned reads per race
races=unique(sapply(tab[,"file"], function(x) unlist(strsplit(x, "_"))[1]))
tabali=NULL
for (ite0 in 1:length(races))
{
	race=races[ite0]
	tab0=tab[grep(race, tab[,"file"]),]
	moy=mean(as.numeric(tab0[,"PCT_PF_READS_ALIGNED"]))
	ic=range(as.numeric(tab0[,"PCT_PF_READS_ALIGNED"]))
	tabali=rbind(tabali, c(nrow(tab0), moy, ic))
}
colnames(tabali)=c("N", "PCAligned_Mean", "PCAligned_Min", "PCAligned_Max")
tabali=data.frame(races, tabali)
write.table(tabali, file="AlignementStatsPerRace.txt", row.names=F, quote=F, sep=" \t")

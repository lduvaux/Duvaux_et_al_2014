rm(list=ls())

# I) check alignment summary metrics
    # I.1) collect stats
        # list files
fils <- dir(pattern="_NoDup.AliSumMetrics$", path="../", full.names = T)

        # read files
tab <- matrix(nrow=length(fils), ncol=25, data=NA)
for (ite0 in 1:length(fils))
{
	ite <- fils[ite0]
	IN <- file(ite,open="r")
	tfil <- readLines(IN) #on slurpe le fichier de sortie de SITES
	close(IN)
	
	if (ite0==1) colnames(tab)=unlist(strsplit(tfil[7], "\t"))
	vstat <- unlist(strsplit(tfil[8], "\t"))
	tab[ite0, 1:length(vstat)]=vstat
	print(ite)
}
tab <- cbind(file=fils, tab)
tab[1:2,]
tab[1:10,1:3]

        # histogram
pdf()
mean(as.numeric(tab[,"PCT_PF_READS_ALIGNED"]))
range(as.numeric(tab[,"PCT_PF_READS_ALIGNED"]))
hist(as.numeric(tab[,"PCT_PF_READS_ALIGNED"]), breaks=20)
dev.off()

    # I.2) looking % aligned reads per race
races <- unique(sapply(tab[,"file"], function(x) unlist(strsplit(x, "_"))[1]))
tabali <- NULL
for (ite0 in 1:length(races))
{
	race <- races[ite0]
	tab0 <- tab[grep(race, tab[,"file"]),]
	moy <- mean(as.numeric(tab0[,"PCT_PF_READS_ALIGNED"]))
	ic <- range(as.numeric(tab0[,"PCT_PF_READS_ALIGNED"]))
	tabali <- rbind(tabali, c(nrow(tab0), moy, ic))
}
colnames(tabali) <- c("N", "PCAligned_Mean", "PCAligned_Min", "PCAligned_Max")
tabali <- data.frame(races, tabali)
write.table(tabali, file="AlignementStatsPerRace.txt", row.names=F, quote=F, sep=" \t")


# II) compute the number of reads (or bases) on baits
fhs <- dir(pattern="HsMetrics$", path="../", full.names = T)
ths <- as.data.frame(matrix(nrow=length(fhs), ncol=4, data=NA))
colnames(ths) <- c("File", "BasesON", "TotalBases", "PCT_BASES_ON")
filnam <- gsub(".AliSumMetrics", "", fils)
ths[,1] <- filnam
if (identical(filnam, gsub(".HsMetrics", "", fhs))){
    for (ite0 in 1:length(fhs))
    {
        ite <- fhs[ite0]
        IN <- file(ite,open="r")
        tfil <- readLines(IN) #on slurpe le fichier de sortie de SITES
        close(IN)
        titles <- unlist(strsplit(tfil[7], "\t"))[-1]
        values <- unlist(strsplit(tfil[8], "\t"))[-1]
        if (ite0==1) {
            ind_ON <- which(titles=="ON_BAIT_BASES")
            ind_TotReads <- which(titles=="TOTAL_READS")
            }
            
        BasesON <- as.numeric(values[ind_ON])
        TotalBases <- as.numeric(values[ind_TotReads]) * as.numeric(tab[ite0, "MEAN_READ_LENGTH"])
        PCT_BASES_ON <- round(BasesON*100/TotalBases, 2)
        ths[ite0, "BasesON"] <- BasesON
        ths[ite0, "TotalBases"] <- TotalBases
        ths[ite0, "PCT_BASES_ON"] <- PCT_BASES_ON
        print(ite)
    }
} else {
    print("WARNING: files not in the same order")
}














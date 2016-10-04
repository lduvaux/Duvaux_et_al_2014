Prepro_fixTargetName <- function(name){
	# fix pipes

	is_promoter <- ifelse(length(grep('^PMT_',name)) > 0, TRUE,FALSE)
	if( !is_promoter)
		out <- sub('\\.','_',name)
	else
		out <- sprintf("%s_1",name)
	return(out)
}

tab <- read.delim("./tablemappingOnV2.LociCategory.Aliases3_20130508.txt")
NewTargetName <- sapply(tab$Aliases, Prepro_fixTargetName)
tab$NewTargetName <- NewTargetName
write.table(tab, file="../MappingOnV2.NewExonName_20131009.txt", sep="\t", row.names=F, quote=T)

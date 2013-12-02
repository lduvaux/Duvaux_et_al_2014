## JUST CONVENIANCE GETTER FUNCTIONS ##
PrePro_findIndex <- function (request, dbb){
	index <- sapply(request, function (x) which(x==dbb))
	return(index)
}

PrePro_fetchRaces <- function(dataFile, ctrl_guys, bad_guys=NULL){
    load(dataFile)
    ind <- PrePro_findIndex(ctrl_guys,names(vrace))
    
    print(vrace)
    print(ind)
    vrace[ind] <- "Medicago_ctrl"
    if (!is.null(bad_guys)[1]){
		ind <- PrePro_findIndex(bad_guys,names(vrace))
		vrace <- vrace[-ind]
	}
    return(vrace)
}

PrePro_fetchUniqRaces <- function(r){
	races <- r
    med <- which(races=="Medicago_ctrl")
    if (length(med) > 0)
		races[med] <- "Medicago"
    return(unique(races))
}

PrePro_fetchRacesPrior <- function(tabrace){
    tabrace <- tabrace[-which(is.na(tabrace[,1])),]
    vrace <- tabrace$Final
    names(vrace) <- tabrace$SampleName
    return(vrace)
}


PrePro_fetchRawDataMat <- function(dataFile, bad_guys=NULL){
    load(dataFile)
    rownames(tested_tab_raw) <- good_targets
    if (!is.null(bad_guys)[1]){
		ind <- PrePro_findIndex(bad_guys, colnames(tested_tab_raw))
		tested_tab_raw <- tested_tab_raw[, -ind]
	}
    return(tested_tab_raw)
}

PrePro_fetchColours <- function(races, races_uniq, race_colo_raw){
    ind <- PrePro_findIndex(races, races_uniq)
    colos <- race_colo_raw[ind]
    names(colos) <- names(races)
    return(colos)
}

PrePro_fetchBaitInfo <- function(datafile){
	temp <- read.delim(datafile, stringsAsFactors =F )
	return(temp)
}

PrePro_fetchRawCtrl <- function(dataFile){
    load(dataFile)
    names(ctrl_vec_raw) <- good_targets
    return(ctrl_vec_raw)
}

PrePro_BaitCateg <- function(vecto){
	temp <- sapply(vecto, function(x) unlist(strsplit(x, "_"))[1])
	return(temp)
}

PrePro_BaitCateg2 <- function(vecto){
	temp <- unlist(mclapply(vecto, function(x) unlist(strsplit(x, "_"))[1]))
	names(temp) <- vecto
	return(temp)
}

PrePro_LibName <- function (indiv_details, indiv){
	lib_name=sapply(1:nrow(indiv_details), function (x) paste(indiv_details[x,"run.date"], indiv_details[x,"flow.cell"], indiv_details[x,"mach.run.nb"], sep="_"))
	
	index_lib <- PrePro_findIndex(indiv, indiv_details[,"ID"])
	lib_name <- lib_name[index_lib]
	
	return(lib_name)
}

PrePro_LibColo <- function (lib_name, col_lib){
	lib_name_uniq <- sort(unique(lib_name))

	index_col <- PrePro_findIndex(lib_name, lib_name_uniq)
	return(col_lib[index_col])
}

PrePro_fetchInfo_TarGene <- function(datafile){
	temp <- read.delim(datafile, stringsAsFactors =F )
	return(temp)
}

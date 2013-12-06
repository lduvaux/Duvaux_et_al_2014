DEBUG <- TRUE
VERBOSE <- 1

PrePro_roundToZeroFive <- function(dat){
    return(round(dat*2)/2)
}


ver <- function(str){
try(
	if(VERBOSE > 0)
		print(str)
	)
}

printprogress <- function(nb, long, step=2.5)
{
	jalon <- seq(0,100, step)
	born <- jalon*long/100
	b1 <- born[born<=nb]
	b1 <- b1[length(b1)]
	
	if (nb==ceiling(b1)) 
		print(paste(round(nb*100/290,2), "% done", sep=""))
}

########### 0) utilities bait names
Pro_ExonName <- function(nom, char_split="_")
{
	elements <- unlist(strsplit(nom, char_split))
	gene <- paste(elements[1:3], collapse="_")
	return(gene)
}

Pro_geneName <- function(nom, char_split="_")
{
	elements <- unlist(strsplit(nom, char_split))
	gene <- paste(elements[1:2], collapse="_")
	return(gene)
}

Prepro_fixName <- function(str){
	# fix pipes
	tmp_strs  <- strsplit(str,'\\|')[[1]]
	name <- tmp_strs[[1]][1]
	
	is_promoter <- ifelse(length(grep('^PMT_',name)) > 0, TRUE,FALSE)
	if( !is_promoter){
		out <- ifelse(length(strsplit(name,'_')[[1]]) == 2, sprintf("%s_1",name),name )	# your forgot the [[1]] here
		out <- sub('\\.','_',out)
		}
	else{
		tmp_strs  <- strsplit(name,'_')[[1]]
		if(length(tmp_strs) == 2){
			out <- sprintf("%s_%s_1_1",tmp_strs[1], tmp_strs[2])
		}
		else{
			out <- sprintf("%s_%s_1_%s",tmp_strs[1], tmp_strs[2], tmp_strs[3])
		}
		
	}
	test <- length(unlist(strsplit(out, "_")))
	if (test!=4)
		stop(paste("There should be 4 elements in: ", out, "(initital name: ", str, ")", sep=""))
	else
		return(out)
}

Draw_pdf <- function(fonction, nom_pdf)
{
	pdf(nom_pdf, onefile = TRUE)
	fonction
	dev.off()
}

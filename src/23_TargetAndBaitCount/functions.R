source("../10_GeneralStatsA/functions.R")
source("../11_GLMTruncFetchData/glm_common.R")

###############
fill_count_table <- function(Count_table, all_TARG_INFO, gene_categ, pro_SEGMENT, baits_BE4_CLEAN)
{
	# 1) initial dataset
		# 1.1) fill targets
	tab0 <- read.delim(all_TARG_INFO, stringsAsFactors=F)
	NTargSet0_Categ <- table(tab0$GeneCateg)
	ind <- PrePro_findIndex(gene_categ[1:length(NTargSet0_Categ)], names(NTargSet0_Categ))
	Count_table[,2] <- NTargSet0_Categ[ind]

		# 1.2) fill genes
	targets_gene0 <- unique(sapply(tab0$NewTargetName, Pro_geneName))
	Count_table[,3] <- table(sapply(targets_gene0, get_element))[ind]

	# 2) after GMAP cleaning
		# 2.1) fill targets
	bad_target_GMAP <- which(is.na(tab0$contigV2)==T)
	tab1 <- tab0[-bad_target_GMAP, ]
	Count_table[,4] <- table(tab1$GeneCateg)[ind]

		# 2.2) fill genes
	targets_gene1 <- unique(sapply(tab1$NewTargetName, Pro_geneName))
	Count_table[,5] <- table(sapply(targets_gene1, get_element))[ind]

		# 2.3) fill complete genes
	IncGenes_GMAP <- unique(sapply(tab0$NewTargetName[bad_target_GMAP], Pro_geneName))
	IncGenes_ind <- targets_gene1%in%IncGenes_GMAP
	NonTrimGenes_GMAP <- targets_gene1[!IncGenes_ind]
	Count_table[,6] <- table(sapply(NonTrimGenes_GMAP, get_element))[ind]

	# 3) final dataset
	load(pro_SEGMENT)	# reference for remaining baits
		# 3.1)  final remaining targets
	ListBaits <- rownames(alpha_matrix)
	exons_ListBaits <- unique(sapply(ListBaits, Pro_ExonName))	# needed to count the exons
	Count_table[,7] <- table(sapply(exons_ListBaits, get_element))[ind]

		# 3.2) final remaining genes
	Tab2gene <- unique(sapply(ListBaits, Pro_geneName))
	Count_table[,9] <- table(sapply(Tab2gene, get_element))[ind]

		# 3.3) the number of complete targets and genes
	refBaits <- read.delim(baits_BE4_CLEAN, stringsAsFactors=F, header=F)[,5]	# reference for baits remaining after GMAP (warning: targets with no overlap!) but before anything else
	refBaits <- sapply(refBaits, Prepro_fixName)
	rmbaits <- which(refBaits%in%ListBaits==F)
				# E) remove incomplete targets
	IncTarg <- unique(sapply(refBaits[rmbaits], Pro_ExonName))	# can partially present or completely absent
	AllTarg <- unique(sapply(refBaits, Pro_ExonName))
	CompTarg <- AllTarg[!AllTarg%in%IncTarg]	# include only targets that are not in IncTarg, i.e. in the list of incomplete targets
	Count_table[,8] <- table(sapply(CompTarg, get_element))[ind]
				# F) remove incomplete genes
	IncGenes_BaitCleaning <- unique(sapply(refBaits[rmbaits], Pro_geneName))
	AllGenes <- unique(sapply(refBaits, Pro_geneName))
	NonTrimGenes_BaitCleaning <- AllGenes[!AllGenes%in%IncGenes_BaitCleaning]
	NonTrimGenes <- NonTrimGenes_BaitCleaning[!NonTrimGenes_BaitCleaning%in%IncGenes_GMAP]
	Count_table[,10] <- table(sapply(NonTrimGenes, get_element))[ind]

	res <- list(Count_table=Count_table, NonTrimGenes=NonTrimGenes)
	
	return(res)
}

get_element <- function(char_string, sep="_", element=1)
{
	return(unlist(strsplit(char_string, sep))[element])
}

###############



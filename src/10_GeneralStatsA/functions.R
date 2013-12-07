############## I) frequency of partially and fullly duplicated GENES per gene category
{
	# I.1) set up the list of result
		# it can be per unit/race or per GeneCateg/race
Pro_test_partial <- function(IndivRace, BaitsNames, alpha_matrix, compute_prop=T) {
	units <- unique(BaitsNames)
	if (compute_prop) {
		mat_ext_dup <- mclapply(units, sub_test_partialGen2, IndivRace, alpha_matrix, BaitsNames)
		names(mat_ext_dup) <- units}
	else
		mat_ext_dup <- sapply(seq(units), function (x) sub_test_partialGen1(units[x], IndivRace, alpha_matrix, BaitsNames))
	return(mat_ext_dup) }

sub_test_partialGen2 <- function(unit0, IndivRace, alpha_matrix, BaitsNames) {
	ind <- which(BaitsNames==unit0)
	mat_1gen <- alpha_matrix[ind,]
	races <- unique(IndivRace)
	vect_gen <- lapply(races, sub_test_partialGenRace, mat_1gen, IndivRace)
	names(vect_gen) <- races
	return(vect_gen)}

sub_test_partialGenRace <- function(race, mat_1gen, IndivRace) {
	ind <- which(IndivRace==race)
	if (is.vector(mat_1gen)) {
		mat_1gen_1race <- mat_1gen[ind]
		vectest <- ifelse(mat_1gen_1race==1, "1_NoDup", "3_CpDup")}
	else {
		mat_1gen_1race <- mat_1gen[,ind]
		vectest <- apply(mat_1gen_1race, 2, sub_test_partialVect)}
	fqcy <- table(vectest)
	test_prop <- c("1_NoDup","2_PtDup","3_CpDup")%in%names(fqcy)
	if (F%in%test_prop) {
		prop_raw <- fqcy
		fqcy <- rep(0,3)
		names(fqcy) <- c("1_NoDup","2_PtDup","3_CpDup")
		fqcy[test_prop] <- prop_raw
	}
	
	prop <- table(vectest)/length(vectest)
	test_prop <- c("1_NoDup","2_PtDup","3_CpDup")%in%names(prop)
	
	if (F%in%test_prop) {
		prop_raw <- prop
		prop <- rep(0,3)
		names(prop) <- c("1_NoDup","2_PtDup","3_CpDup")
		prop[test_prop] <- prop_raw
	}
	ll <- list(vectest=vectest, fqcy=fqcy, prop=prop)
	return(ll)
	}

		# I.1.1) test
sub_test_partialVect <- function(vect) {	# test one individual
	unic <- unique(vect)
	if (length(unic)!=1)
		test <- "2_PtDup"
	else if (unic==1)
		test <- "1_NoDup"
	else
		test <- "3_CpDup"   # not necessarily the same number of copies all along!
	return(test)}

		# I.1.2) coding duplications and computing proportion simultaneously
sub_test_partialGen1 <- function(unit0, IndivRace, alpha_matrix, BaitsNames) {
	ind <- which(BaitsNames==unit0)
	mat_1gen <- alpha_matrix[ind,]
	races <- unique(IndivRace)
	if (is.vector(mat_1gen))
		vect_gen <- ifelse(mat_1gen==1, "1_NoDup","3_CpDup")
	else 
		vect_gen <- sub_test_partialGentTest(mat_1gen)
	return(vect_gen)}

sub_test_partialGentTest <- function(mat_1gen, IndivRace) {
	vectest <- apply(mat_1gen, 2, sub_test_partialVect)
	return(vectest)}


# 1.3) set up summary matrix per GenCateg
Pro_ExtDup <- function(liste, Genes_GeneCateg){
	GeneCategs <- unique(Genes_GeneCateg[,2])
	list_extDup_GeneCateg <- lapply(GeneCategs, sub_ExtDup_1GenCateg_Allrace, liste, Genes_GeneCateg)
	names(list_extDup_GeneCateg) <- GeneCategs
	return(list_extDup_GeneCateg)}

sub_ExtDup_1Gen_Allrace <- function(subliste) {
	mat_sub_ExtDup_1Gen_Allrace <- t(sapply(subliste, function(sol) sol$fqcy))
	return(mat_sub_ExtDup_1Gen_Allrace)}

sub_ExtDup_1GenCateg_Allrace <- function(GeneCateg, liste, Genes_GeneCateg) {
	ind <- which(Genes_GeneCateg[,2]==GeneCateg)
	sublist <- liste[ind]
	for(ite in seq(sublist)) {
		if (ite==1) 
			mat_extDup_AllGen_Allrace <- sub_ExtDup_1Gen_Allrace(sublist[[ite]])
		else
			mat_extDup_AllGen_Allrace <- mat_extDup_AllGen_Allrace+sub_ExtDup_1Gen_Allrace(sublist[[ite]])}
	return(mat_extDup_AllGen_Allrace)}
Pro_ExtDup_AllCateg <- function(liste){
	for(ite in seq(liste)) {
		if (ite==1) 
			fin_list <- liste[[ite]]
		else
			fin_list <- fin_list+liste[[ite]]}
	return(fin_list)}

	####
Pro_setup_list_extDup <- function(mat_extDup_AllCateg, list_extDup_GeneCateg) {
	list_extDup <- vector("list", length(list_extDup_GeneCateg)+1)
	list_extDup[[1]] <- mat_extDup_AllCateg
	list_extDup[2:(length(list_extDup_GeneCateg)+1)] <- list_extDup_GeneCateg
	names(list_extDup)[1] <- "AllCateg"
	names(list_extDup)[2:(length(list_extDup_GeneCateg)+1)] <- names(list_extDup_GeneCateg)
	return(list_extDup)}

	####
Pro_setup_mat_extDup <- function(list_extDup) {
	mat <- matrix(NA, 3, length(list_extDup))
	for(ite in seq(list_extDup)) {
		vect <- apply(list_extDup[[ite]], 2, sum)
		mat[,ite] <- vect}
	colnames(mat) <- names(list_extDup)
	rownames(mat) <- c("1_NoDup","2_PtDup","3_CpDup")
	return(mat)}
}

############## II) frequency of duplicated BAITS/genes per gene category

#~pro_detectCNV <- function(categ, BaitsGeneCateg, alpha_matrix, fqcy=F, raw_scor=F) {
#~	index_baits <- which(BaitsGeneCateg==categ)
#~	tab <- alpha_matrix[index_baits,]
#~	score <- sum(apply(tab, 1, sub_get_cnv_score, raw_scor))	# look each bait independentlty
#~	res <- ifelse(fqcy, score/length(index_baits), score)
#~	return(res)}

#~sub_get_cnv_score <- function(vect, raw_scor){	# if raw_score ==T, detect polymorphic units (binary results)
#~	res <- ifelse(raw_scor, !(length(unique(vect))==1 & unique(vect)[1]==1), sum(vect!=1))	# if the number of possible states is >1 or if the number of states is 1 but with a value !=1, thus there is CNV
#~	return(res)}

	####
#~pro_detectCNV_RaceCateg <- function(categ, BaitsGeneCateg, race, IndivRace, alpha_matrix, fqcy=F, raw_scor=F) {
#~	index_baits <- which(BaitsGeneCateg==categ)
#~	index_race <- which(race==IndivRace)
#~	tab <- alpha_matrix[index_baits,index_race]
#~	score <- sum(apply(tab, 1, sub_get_cnv_score, raw_scor))
#~	res <- ifelse(fqcy, score/length(index_baits), score)
#~	return(res)}

	####
pro_polymGene <- function(gene, BaitsGeneNames, alpha_matrix, fqcy=F, blocks=F, groups, CpDup_Only) {
#~    print(gene)
	index_baits <- which(BaitsGeneNames==gene)
	tab <- alpha_matrix[index_baits,]

	if (blocks) {
		if (blocks && missing(groups))
			stop("Vector of individual's groups is absent")
        if (fqcy)
            res <- numeric(length(groups))
        else
            res <- character(length(groups))

		for (ite in 1:length(groups)){
			gp <- groups[[ite]]
			grep_patt <- paste(paste("^", gp, sep=""), collapse="|")
#~			print(gp)
			
			if (is.vector(tab)) {
                indiv <- grep(grep_patt, names(tab))
				tab2 <- tab[indiv]}
			else {
                indiv <- grep(grep_patt, colnames(tab))
                tab2 <- tab[,indiv]}
            if (fqcy)
                res[ite] <- test_pro_fqcy(tab2, gp, gene, CpDup_Only, blk=blocks)
            else
                res[ite] <- test_pro_polymGene(tab2)}
	} else {
        if (fqcy)
            res <- test_pro_fqcy(tab, gp, gene, CpDup_Only)
        else
            res <- test_pro_polymGene(tab)}
	return(res)}

test_pro_polymGene <- function(tab_2)
{
	if (is.matrix(tab_2)) {
		vec_test <- apply(tab_2, 2, sub_test_partialVect)	# test the polymorphism of each indiv
		res <- sub_get_polymGene(vec_test)}
	else {
		test <- unique(tab_2)
		res <- ifelse(length(test)==1 & test[1]==1, "1_NoDup", "3_CpDup")}	# if length==1 and the only value is 1, sothere is no CNV (value=0), if not because there is only one bait so the gene is completely duplicated for at least one indiv (value=1)
	return(res)}

sub_get_polymGene <- function(vec_test){
	if ("3_CpDup"%in%vec_test)
		res <- "3_CpDup"
	else if ("2_PtDup"%in%vec_test)
		res <- "2_PtDup"
	else
		res <- "1_NoDup"
	return(res)}
    
test_pro_fqcy <- function(tab, gp, gene, CpDup_Only, blk)
{
    # compute fqcy of all kind of duplications (partial & complete)
    if (!CpDup_Only) {
#~        print("All_Dup")
        if (is.matrix(tab)) {
            vec_test <- apply(tab, 2, sub_test_partialVect)	# test the polymorphism of each indiv
            lg_CN1 <- length(grep("^1_NoDup$", vec_test))}
        else {
            vec_test <- tab
            lg_CN1 <- length(which(tab==1))}
            
        sub_res <- round(1-(lg_CN1/length(vec_test)),3)  # result for all kind of duplications
        if (sub_res==1) {
            if (blk)
                print(paste("No individual of the group ", gp, " with CN=1 for the gene ", gene, sep=""))
            else
                print(paste("No individual with CN=1 for the gene ", gene, sep=""))}
        }

    # compute fqcy of complete duplication only
    else {
#~        print("Complete Dup only")
        if (is.matrix(tab)) {
            vec_test <- apply(tab, 2, sub_test_partialVect)
            lg_CpDup <- length(grep("^3_CpDup$", vec_test))
            }
        else {
            vec_test <- tab
            lg_CpDup <- length(which(tab!=1))}

        sub_res <- round(lg_CpDup/length(vec_test),3)    # result for all complete duplications
        }
	return(sub_res)}




############## III) polymorphism within and between races
{
	############# III.1) subsample the main matrix
	{
pro_distCN_per_Race <- function(race, IndivRace, alpha_matrix){
	tab <- sub_subset_tab(race, IndivRace, alpha_matrix, byrow=F)
	distr <- table(tab)
	return(distr)}

sub_subset_tab <- function(request, database, matrice, byrow=T){
	ind <- which(database==request)
	if (byrow)
		mat <- matrice[ind,]
	else
		mat <- matrice[,ind]
	return(mat)}

	####
pro_distCN_per_GeneCateg <- function(categ, BaitsGeneCateg, alpha_matrix){
	tab <- sub_subset_tab(categ, BaitsGeneCateg, alpha_matrix, byrow=T)
	distr <- table(tab)
	return(distr)}
	}

	############# III.2) autapomorphies
	{
pro_mat_autapo <- function(IndivRace, categ, BaitsGeneCateg, alpha_matrix){
	print(categ)
	tab <- sub_subset_tab(categ, BaitsGeneCateg, alpha_matrix, byrow=T)
	spec_categ <- apply(tab, 1, sub_count_spec, IndivRace)
	mat_spec <- sub_SetSpec_mat(tab, spec_categ)
	return(mat_spec)}
	
sub_count_spec <- function(bait, IndivRace, detect_fixed=F)
{
	distrib_matrix <- table(bait, IndivRace)
	raw_spec_bait <- apply(distrib_matrix, 1, sub_detect_spec)
	
	if (detect_fixed) {
		autapo <- apply(distrib_matrix, 1, sub_detect_alone)
		if (T%in%autapo){
			synapo <- rep("None", length(autapo)); names(synapo) <- as.character(names(raw_spec_bait))
			ind <- which(autapo==T)
			for (ite in ind)
			{
				vect <- distrib_matrix[ite,]
				ind_indiv <- which(vect!=0)
				test <- sub_detect_alone(distrib_matrix[, ind_indiv])
				if (test) {
					synapo[ite] <- colnames(distrib_matrix)[ind_indiv]
					names(synapo)[ite] <- as.character(names(raw_spec_bait))[ite]}
			}
		}
		else
			{synapo <- "None"; names(synapo)="1"}
		return(synapo)}
	else
		return(raw_spec_bait)
}

sub_detect_spec <- function(vect, IndivRace)
{
	auta <- sub_detect_alone(vect)
	if (auta) spec=names(vect)[which(vect!=0)] else spec="None"
	return(spec)
}

sub_detect_alone <- function(vect)
{
	zeros <- which(vect==0)
	lon <- length(zeros)
	if (lon==(length(vect)-1)) res <- T else res <- F
	return(res)
}

sub_SetSpec_mat <- function(tab, spec_categ)
{
	seqq <- seq(0,5, 0.5)
	mat <- matrix(ncol=length(seqq), nrow=nrow(tab), data=NA)
	colnames(mat)=seqq; rownames(mat)=rownames(tab)
	
	if (class(spec_categ)=="list"){
		for(i in 1:length(spec_categ))
		{
			# print(i)
			ind <- PrePro_findIndex(names(spec_categ[[i]]), seqq)
			mat[i,ind]=as.character(spec_categ[[i]])
		}
	}
	else if (class(spec_categ)=="character")
		mat[,"1"]="None"
	else 
		stop("unexpected features: fix the function")
	return(mat)
}

	####
pro_spec_perRace <- function(race, mat_spec) {
	ind <- which(mat_spec==race, arr.ind=T)
	res <- length(unique(ind[,1]))
	return(res)}

	####
pro_div <- function(x, vect) {x/vect}
	}

	############# III.3) autapomorphies
	{
pro_mat_synapo <- function(IndivRace, categ, BaitsGeneCateg, alpha_matrix){
	print(categ)
	tab <- sub_subset_tab(categ, BaitsGeneCateg, alpha_matrix, byrow=T)
	spec_categ <- apply(tab, 1, sub_count_spec, IndivRace, detect_fixed=T)
	mat_spec <- sub_SetSpec_mat(tab, spec_categ)
	return(mat_spec)}
}
	}

	############# III.4) Compute 'genotype' frequencies
Pro_geno_allRaces <- function(races, alpha_matrix, IndivRace)
{
	MostFqAltGeno <- matrix(ncol=length(races), nrow=nrow(alpha_matrix), data=NA)	# most frequent Alternative genotype
	FqAltGeno <- matrix(ncol=length(races), nrow=nrow(alpha_matrix), data=NA)	# fqce alternative genotype
	rownames(MostFqAltGeno) <- rownames(alpha_matrix); colnames(MostFqAltGeno) <- races
	rownames(FqAltGeno) <- rownames(alpha_matrix); colnames(FqAltGeno) <- races
	
	for (ite in 1:length(races))
	{
		mat <- sub_fq_allBaits(ite, races, alpha_matrix, IndivRace, MostFqAltGeno, FqAltGeno)
		MostFqAltGeno[,ite] <- mat[,2]
		FqAltGeno[,ite] <- mat[,1]
	}
	res <- list(MostFqAltGeno=MostFqAltGeno, FqAltGeno=FqAltGeno)
	return(res)
}

sub_fq_allBaits <- function(ite, races, alpha_matrix, IndivRace, MostFqAltGeno, FqAltGeno)
{
	race <- races[ite]
	print(race)
	mat <- t(sapply(1:nrow(alpha_matrix), function(x) sub_geno_fq(alpha_matrix[x,], race, IndivRace)))
	return(mat)
}

sub_geno_fq <- function(bait, race, IndivRace)#, x)
{
	ind <- grep(race, IndivRace)
	distrib_matrix <- table(bait[ind], IndivRace[ind])
	maxfq_no1 <- sub_maxfq_no1(distrib_matrix[,1])
#	if (length(maxfq_no1)==0) stop(x)
	return(maxfq_no1)
}

sub_maxfq_no1 <- function(vec_fq)
{
	if (length(vec_fq)>1) {
		bad <- which(names(vec_fq)==1)
		if (length(bad)>0) 
			vec_fq_gd <- vec_fq[-bad]
		else 
			vec_fq_gd <- vec_fq
		
		good <- vec_fq_gd[which.max(vec_fq_gd)]
		fq <- good/sum(vec_fq)
		genotype <- as.numeric(names(good))
		}
	else
		{fq <- 0; genotype <-  (-1)}
	ll <- c(round(fq,2), round(genotype,1))
	return(ll)
}












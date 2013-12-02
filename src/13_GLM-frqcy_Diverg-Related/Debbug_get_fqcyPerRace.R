source("./functions.R")

# test set_Pre_GLMtable function (folder 11_)
CLUSTERS <- list(Cytisus="Cytisus", Lathyrus="Lathyrus", L.corn.="L.corn.", L.ped.="L.ped.", Ononis="Ononis", Medicago="Medicago", Pisum="Pisum", Trifolium="Trifolium")

list_gen <- ListGenes
Baits_Gn_Name <- BaitsGeneNames
alpha_mat <- alpha_matrix
Tab_Genes_Info <- Genes_Info
Categ_4GLM <- CATEG_FOR_GLM
CompGenes <- CompGenes
test_blocks <- T
list_groups <- CLUSTERS


# test pro_polymGene function (folder 10_)
gene <- "Control_g209"
fqcy <- T
blocks <- T
groups <- list_groups

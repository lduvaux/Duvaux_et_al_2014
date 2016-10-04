    inf_PMT <- names(gini_gene_rf)[grep("^PMT_", names(gini_gene_rf))]
    N_inf_PMT <- length(inf_PMT)

    N_pairs_inf_PMT <- (length(inf_PMT)^2-length(inf_PMT))/2
    N_pairs_same_cont_PMT <- round(P_PMT*N_pairs_inf_PMT)
    N_pairs_dif_cont_PMT <- N_pairs_inf_PMT - N_pairs_same_cont_PMT

    ind_same_contig_PMT <- sample((1:length(vec_test_cont_PMT))[vec_test_cont_PMT], N_pairs_same_cont_PMT)
    ind_dif_contig_PMT <- sample((1:length(vec_test_cont_PMT))[!vec_test_cont_PMT], N_pairs_dif_cont_PMT)

    vec <- vec_cont_val_PMT[c(ind_same_contig_PMT, ind_dif_contig_PMT)]
    chosen_PMT <- unique(unlist(strsplit(vec, ":")))









      P_PMT <- as.numeric(P_PMT)
    info_TargGene <- read.delim(info_TargGene_fil)
    # look for contigs
    v <- sapply(bait_names, collapse_elements, sep="_", what=1:3, colla="_")
    vcontigs <- sapply(v, get_contig, info_TargGene)

    # 1) draw PMT
        # test same contig for all possible pairs of baits (bait_names)
    PMTs_ind <- grep("^PMT_", bait_names)
    vcontigs_PMT <- vcontigs[PMTs_ind]
    bait_names_PMT <- bait_names[PMTs_ind]

        # test pair same contig
    mat_test_cont_PMT <- outer(vcontigs_PMT, vcontigs_PMT, "==")
    vec_test_cont_PMT <- mat_test_cont_PMT[upper.tri(mat_test_cont_PMT, diag = FALSE)]

        # all pair compositions
    mat_cont_val_PMT <- outer(bait_names_PMT, bait_names_PMT, paste, sep=":")
    vec_cont_val_PMT <- mat_cont_val_PMT[upper.tri(mat_cont_val_PMT, diag = FALSE)]

        # set up the number of pairs from the same contig to draw
    inf_PMT <- names(gini_gene_rf)[grep("^PMT_", names(gini_gene_rf))]
    v_PMT <- rep(NA, length(inf_PMT))
    while (NA%in%v_PMT)
    {
        n_same <- round(P_PMT*100)
        n_dif <- round((1-P_PMT)*100)
        ind_same_contig_PMT <- sample((1:length(vec_test_cont_PMT))[vec_test_cont_PMT], n_same)
        ind_dif_contig_PMT <- sample((1:length(vec_test_cont_PMT))[!vec_test_cont_PMT], n_dif)
        vec <- vec_cont_val_PMT[c(ind_same_contig_PMT, ind_dif_contig_PMT)]
        chosen_PMT <- unique(unlist(strsplit(vec, ":")))

    }


system.time(
    get_baits_per_pairs(bait_names=bait_nam, P_PMT=P_PMT_cont_PMT_obs, P_Gn=P_Gn_cont_Gn_obs, P_Gn_PMT=P_Gn_cont_PMT_obs, info_TargGene_fil=INFO_TARGENE_FILE, gini_gene_rf=gini_gene_rf, ds_pmt=ds_PMT, ds_gns=ds_Gn, inc=10, inc2=5, inc3=5, verbose=1)
)

system.time(
    sims <- sapply(1:100, function (x)
        get_baits_per_pairs(x, bait_names=bait_nam, P_PMT=P_PMT_cont_PMT_obs, P_Gn=P_Gn_cont_Gn_obs, P_Gn_PMT=P_Gn_cont_PMT_obs, info_TargGene_fil=INFO_TARGENE_FILE, gini_gene_rf=gini_gene_rf, ds_pmt=ds_PMT, ds_gns=ds_Gn, inc=10, inc2=5, inc3=5, verbose=0)$Gns
    )
)

print(aa <- sapply(1:ncol(sims), function(x) get_P_same_contig(sims[,x], sims[,x], INFO_TARGENE_FILE, verbose=F)))

sapply(2:ncol(sims), function(x,y) identical(sort(sims[,1]), sort(sims[,x])))

iden <- function(x, y, tab=sims)
{
    res <- sapply(seq(y), function(a) identical(sort(sims[,x]), sort(sims[,a])))
    return(res)
}

res <- sapply(1:ncol(sims), iden, 1:ncol(sims), tab=sims)

system.time(
    sims2 <- sapply(1:100, function (x)
        get_baits_per_pairs(x, bait_names=bait_nam, P_PMT=P_PMT_cont_PMT_obs, P_Gn=P_Gn_cont_Gn_obs, P_Gn_PMT=P_Gn_cont_PMT_obs, info_TargGene_fil=INFO_TARGENE_FILE, gini_gene_rf=gini_gene_rf, ds_pmt=ds_PMT, ds_gns=ds_Gn, inc=10, inc2=5, inc3=5, verbose=0)$PMTs
    )
)




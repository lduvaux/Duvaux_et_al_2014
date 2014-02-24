makeTransparent<-function(someColor, alpha=50)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

draw_star_bar0 <- function(tab, i, yl, colo, alph=100)
{
    polygon(x=c(tab$startV2[i],tab$stopV2[i],tab$stopV2[i],tab$startV2[i]), y=c(yl[1], yl[1], yl[2], yl[2]),col=makeTransparent(colo, alpha=alph), border = NA)
}

draw_star_bar <- function(tab, i, yl, colo, alph=100)
{
    polygon(x=c(tab$Start[i],tab$End[i],tab$End[i],tab$Start[i]), y=c(yl[1], yl[1], yl[2], yl[2]),col=makeTransparent(colo, alpha=alph), border = NA)
}

get_link_xcoor <- function(cnv_targ, ftab_cnv)
{
    # subset table
    x_coord <- sort(unlist(lapply(cnv_targ, function(x) range(ftab_cnv[grep(paste(x, "_", sep=""), rownames(ftab_cnv)),"bary"]))))    # define range of cnv_targets
    x_coord <- x_coord[-c(1, length(x_coord))]# remove first and last points as not used for linking
    return(x_coord)
}

draw_lines <- function(v_sta, v_end, l_wd=0.5, l_ty=3, colo="red")
{
    points(x=v_sta, y=v_end, lwd=l_wd, lty=l_ty, type="l", col=colo)
}

plot_CNV_chr <- function(tab_star, tab_tar, tab_cnv, yli=c(-0.5, 2.5), centz=c(0.75,1.25), c_ex=.9, l_wd=.9, races, transp=F, alphaa=80)
#zcent: central zone, use the delimit the rounding area where alpha wil lbe rounded to 1
{
    # 1) define preliminary parameters
        # 1.0) misc
    sub_targ <- tab_star$Name
    contig <- unique(tab_star$Contig)
    cnv_starg <- rownames(tab_cnv)
    gene <- collapse_elements(sub_targ[1], what=1:2)
    print(gene)
    targ <- tab_tar$NewTargetName
    cnv_targ <- unique(sapply(cnv_starg, collapse_elements))

    races_uniq <- unique(races)
    col_races <- PrePro_fetchColours(races, races_uniq, race_colo_raw)
        # 1.1) table of subtarg present and absent from the final dataset
    conc <- tab_star[,"Name"]%in%rownames(tab_cnv)
    tab_star_gd <- tab_star[conc, ]
    tab_star_gd[,"bary"] <- apply(tab_star_gd[, c("Start", "End")], 1, median)

        # 1.2) set up final tab of CNV
#~    ind <- sapply(rownames(tab_cnv), function(x) which(x==tab_star_gd$Name))
    ind <- match(rownames(tab_cnv), tab_star_gd$Name)
    ftab_cnv <- cbind(bary=tab_star_gd[ind,"bary"], tab_cnv)
    ftab_cnv <- ftab_cnv[order(ftab_cnv[, "bary"]), ]

    # 2) plot an empty graph
    rg_coo <- c(min(c(tab_tar$startV2,tab_tar$stopV2), na.rm=T), max(c(tab_tar$startV2,tab_tar$stopV2), na.rm=T))# define range of initial targets
    print("Range of displayed marker:")
    print(range(rg_coo))
    plot(rg_coo, rep(1, 2), xlim=rg_coo, ylim=yli, type="n", main=paste(gene, " (", contig, ")", sep=""), xlab="Contig coordinate (bp)", ylab="Alpha (CNV relative to standard)")
    points(x=c(-10, 10000000), y=c(0,0), type="l")	# add a lower box (x range is over large on purpose) 

    # 3) draw the rounding area
    polygon(x=c(-10, 10000000, 10000000, -10), y=c(centz[1], centz[1], centz[2], centz[2]), col=makeTransparent("red", alpha=60), border = NA)

    # 4) draw the spanning of subtargets as originally designed
    for (sta in 1:nrow(tab_tar))
    {draw_star_bar0(tab=tab_tar, i=sta, yl=yli+c(-0.5, +0.5), colo="gray", alph=80)}

    # 5) draw
        # the initial sets of subtargets...
    for (sta in 1:nrow(tab_star))
    {draw_star_bar(tab=tab_star, i=sta, yl=c(yli[1]-0.5, 0), colo="black", alph=90)}
        # ... and the final sets of subtargets
    if(nrow(tab_star_gd)>0){
        for (sta in 1:nrow(tab_star_gd))
        {draw_star_bar(tab=tab_star_gd, i=sta, yl=c(yli[1]-0.5, 0), colo="darkgoldenrod2", alph=100)}
    }

    # 6) set up the vector of linking lines between targets
    if (length(cnv_targ)>1) 
        x_coord <- get_link_xcoor(cnv_targ, ftab_cnv)

    # 7) plot alpha segments for each individual
    if (is.vector(ftab_cnv)){
        ftab_cnv <- t(as.matrix(ftab_cnv))}

    med_ctrl <- names(races)[races%in%"Medicago_ctrl"]
    ctrl_in_ftab_cnv <- sapply(med_ctrl, function(x) which(colnames(ftab_cnv)==x))
    v_rdom_draw <- (2:ncol(ftab_cnv))[-(ctrl_in_ftab_cnv-1)]
    rdom_draw <- c(sample(v_rdom_draw, length(v_rdom_draw)), ctrl_in_ftab_cnv)
    for (i in rdom_draw)
    {
        # 7.1) set up individual parmater
        indiv <- colnames(ftab_cnv)[i]
        if (transp)
            colo <- makeTransparent(col_races[which(indiv==names(col_races))], alphaa)
        else
            colo <- col_races[which(indiv==names(col_races))]

        # 7.2) draw linking lines between targets
        if (length(cnv_targ)>1) {
            x_coord_ind <- sapply(x_coord, grep, ftab_cnv[,"bary"])
            y_coord <- ftab_cnv[x_coord_ind,i]
            sapply(seq(1, length(x_coord), by=2), function(x) draw_lines(v_sta=x_coord[x:(x+1)], v_end=y_coord[x:(x+1)], colo=colo))
        }

        # 7.3) draw lines between subtargets of the same targets
        if (length(cnv_targ)>1) {
            x_coord_sta0 <- c(1, x_coord_ind)
            x_coord_sta <- x_coord_sta0[seq(1, length(x_coord_sta0), by=2)]
            x_coord_end0 <- c(x_coord_ind, nrow(ftab_cnv))
            x_coord_end <- x_coord_end0[seq(1, length(x_coord_end0), by=2)]
            sapply(seq(length(x_coord_sta)), function(x)
                points(x=ftab_cnv[x_coord_sta[x]:x_coord_end[x],1], y=ftab_cnv[x_coord_sta[x]:x_coord_end[x],i], col=colo, type="o", pch=20, lwd=l_wd, cex=c_ex)
            )
        }
        else {
            points(x=ftab_cnv[,1], y=ftab_cnv[,i], col=colo, type="o", pch=20, lwd=l_wd, cex=c_ex)
        }
    }
    cat("\n")

    return(col_races)
}

get_data4plot <- function(gn, t_targ, subtarg, alpha_matrix, adapt=T)
{
            ind <- grep (paste(gn, "_", sep=""), t_targ[,"NewTargetName"])
            lt_targ <- t_targ[ind, ]

            ind <- grep (paste(gn, "_", sep=""), subtarg[,"Name"])
            lt_star <- subtarg[ind, ]
            n_ctig <- sort(table(lt_star$Contig), decreasing=T)
            
            ind <- grep (paste(gn, "_", sep=""), rownames(alpha_matrix))
            lcnv <- alpha_matrix[ind,]
            if (is.vector(lcnv)) {
                lcnv <- t(as.matrix(lcnv))
                rownames(lcnv) <- rownames(alpha_matrix)[ind]
            }
            
            # check contig
            if (length(n_ctig)>1) {
                print("WARNING: gene on several contigs")
                
                # detect good subtarget in initial subtargets
                gd_ctig <- names(n_ctig)[1]
                ind <- lt_star$Contig%in%gd_ctig
                lt_star <- lt_star[ind,]    # subset initial table of subtargets
                
                # remove subtargets from bad contigs in lcnv
                    # fetch contigs of star in lt_cnv
#~                ind <- match(rownames(lcnv), subtarg[, "Name"])
                ind <- sapply(rownames(lcnv), function(x) which(x==subtarg[, "Name"]))
                cnv_contigs <- subtarg[ind,"Contig"]
                ind <- cnv_contigs%in%gd_ctig
                lcnv <- lcnv[ind, ]

                # same with lt_targ table
                ind <- lt_targ$contigV2%in%gd_ctig
                lt_targ <- lt_targ[ind,]
            }
            LengthV1 <- max(c(lt_targ$startV1,lt_targ$stopV1), na.rm=T) - min(c(lt_targ$startV1,lt_targ$stopV1), na.rm=T) + 1
            LengthV2 <- max(c(lt_targ$startV2,lt_targ$stopV2), na.rm=T) - min(c(lt_targ$startV2,lt_targ$stopV2), na.rm=T) + 1
            Fold <- round(LengthV1/LengthV2, 2)
            stats <- c(gn, LengthV1, LengthV2, Fold)

            mmax <- max(lcnv)
            if (adapt)
                rgg <- c(-0.5, ceiling(mmax))
            else {
                if (mmax< 2.5)
                    rgg <- c(-0.5, 2.5)
                else
                    rgg <- c(-0.5, ceiling(mmax))
            }
        res <- list(lcnv=lcnv, lt_star=lt_star, lt_targ=lt_targ, rgg=rgg, stats=stats)
}

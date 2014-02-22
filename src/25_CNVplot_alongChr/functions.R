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

plot_CNV_chr <- function(tab_cnv, tab_star, tab_tar, yli=c(-0.5, 2.5), centz=c(0.75,1.25))
#zcent: central zone, use the delimit the rounding area where alpha wil lbe rounded to 1
{
    # 1) define sub-tables
        # 1.1) table of subtarg present and absent from the final dataset
    conc <- rownames(tab_cnv)%in%tab_star[,"Name"]
    tab_star_gd <- tab_star[conc, ]
    tab_star_gd[,"bary"] <- apply(tab_star_gd[, c("Start", "End")], 1, median)
    ind <- sapply(rownames(tab_cnv), function(x) which(x==tab_star_gd$Name))

        # 1.2) set up final tab of CNV
    ind <- sapply(rownames(tab_cnv), function(x) which(x==tab_star_gd$Name))
    ftab_cnv <- cbind(bary=tab_star_gd[ind,"bary"], tab_cnv)
    ftab_cnv <- ftab_cnv[order(ftab_cnv[, "bary"]), ]

    # 2) plot an empty graph
    rg_coo <- c(min(tab_tar$startV2), max(tab_tar$stopV2))# define range of initial targets

    plot(rg_coo, rep(1, 2), xlim=rg_coo, ylim=yli, type="n", main=gn, xlab="Chromosome coordinate (bp)", ylab="Alpha (CNV relative to standard)")
    points(x=c(-10, 10000000), y=c(0,0), type="l")	# add a lower box (x range is over large on purpose) 

    # 3) draw the rounding area
    polygon(x=c(-10, 10000000, 10000000, -10), y=c(centz[1], centz[1], centz[2], centz[2]), col=makeTransparent("red", alpha=60), border = NA)

    # 4) draw the spanning of subtargets as originally designed
    for (sta in 1:nrow(tab_tar))
    {draw_star_bar0(tab=tab_tar, i=sta, yl=yli+c(-0.5, +0.5), colo="gray", alph=80)}

    # 5) draw
        # the initial 
    for (sta in 1:nrow(tab_star))
    {draw_star_bar(tab=tab_star, i=sta, yl=c(yli[1]-0.5, 0), colo="black", alph=90)}
        # and final sets of subtargets
    if(nrow(tab_star_gd)>0){
        for (sta in 1:nrow(tab_star_gd))
        {draw_star_bar(tab=tab_star_gd, i=sta, yl=c(yli[1]-0.5, 0), colo="darkgoldenrod2", alph=100)}
    }

    # 6) plot alpha segments for each individual
    for (i in 2:ncol(ftab_cnv))
    {
        
    }






#######################################
#######################################
#######################################
#######################################






    for (ilist in 1:length(newxsta))
    {
        # set up vectors for plotting
        dista <- sapply(2:length(xsta), function(z) xsta[z]-xend[z-1]-1)
        
        # set up variables and tables
            # all initially captured targets on this contig
        ista1=which.max(stnoverl0$Start[stnoverl0$Start<=xsta[1]])
        if (xsta[length(xsta)]>stnoverl0$Start[length(stnoverl0$Start)]) iend1=length(stnoverl0$Start) else iend1=which(stnoverl0$Start==min(stnoverl0$Start[stnoverl0$Start>=xsta[length(xsta)]]))
        stnoverl=stnoverl0[ista1:iend1,]	# all initially captured targets on this contig
        svgen=unique(sapply(stnoverl$targetID, function(x) unlist(strsplit(x, ".", fixed=T))[1]))
        
            # sub-table of sub-targets remaining after the first sift on this contig
        ista3=which.max(sgootarg0$start[sgootarg0$start<=xsta[1]])
        if(xsta[length(xsta)]>sgootarg0$start[length(sgootarg0$start)]) iend3=length(sgootarg0$start) else iend3=which(sgootarg0$start==min(sgootarg0$start[sgootarg0$start>=xsta[length(xsta)]]))
        sgootarg=sgootarg0[ista3:iend3,]	# sub-table of sub-targets remaining after the first sift on this contig
            # index of good targets on this contig
        indgd=indgd0[ista3:iend3]
        
            # number of different targets (no sub-targets)
        itargs=which(dista>1)
        n.targ=length(itargs)+1	# different from nrow(stnoverl)
            
            # xlim for graphs
        xli=range(c(sgootarg$start, sgootarg$end))	# xlim of plot
            # set up individual graphic options
        pp=rep(20, 120)#pp=rep(1:10, 12)
        
        for (indiv in 1:ncol(talpha))
        {
            # set up reference individuals
            race0=vrace[indiv]
            icolo=grep(race0, race)
            nomindgd=vindiv[indiv]
            if (nomindgd%in%ref==T) icolo=9
            for (n in 1:n.targ)
            {
                # set up subsets of sub-targets per target
                if (n==1) targ1=1 else targ1=itargs[n-1]+1
                targ2=itargs[n]
                if (n==n.targ) targ2=length(xsta)
                x=sapply(targ1:targ2, function (x) round((sgootarg$end[x]+sgootarg$start[x])/2, 0))
                y=talpha[indgd[targ1:targ2],indiv]
                
                # draw une linking line between targets
                if (n!=1) points(x=c(lastx, x[1]), y=c(lasty, y[1]), lwd=0.5, lty=3, type="l", col=colos0[icolo])
                
                lastx=x[length(x)]
                lasty=y[length(y)]
                
                if (n==1 & indiv==1) 
                {
                    # plot an empty graph
                    if (length(svgen)>6) 
                    {
                        lsvgen=length(svgen)
                        midsvgen=ceiling((3*lsvgen)/8)
                        intiti=paste(paste(svgen[1:midsvgen], collapse=";"), "\n", paste(svgen[(midsvgen+1):lsvgen], collapse=";"), sep="")
                        titi=paste(contig, " (", intiti, ")", sep="")
                    } else titi=paste(contig, " (", paste(svgen, collapse=";"), ")", sep="")
                    plot(x, y, xlim=xli, ylim=yli, col=colos[icolo], type="n", pch=pp[indiv], main=titi, xlab="position (pb)", ylab="alpha (nb of copies)")
                    points(c(-10, 10000000), c(0,0), type="l")	# add a lower box
                    # draw the confidence area
                    polygon(x=c(-10, 10000000, 10000000, -10), y=c(bornes[1], bornes[1], bornes[2], bornes[2]), col=makeTransparent("red", alpha=60), border = NA)
                    # draw the spanning of targets as originally designed
                    for (tata in 1:nrow(stnoverl))
                    {
                        polygon(x=c(stnoverl$Start[tata],stnoverl$End[tata],stnoverl$End[tata],stnoverl$Start[tata]), y=c(yli[1]-0.5,yli[1]-0.5, yli[2]+0.5, yli[2]+0.5), col=makeTransparent("gray", alpha=80), border = NA)
                    }
                    # draw the sub-targets remaining after first sift
                    for (tata in 1:nrow(sgootarg))
                    {
                        polygon(x=c(sgootarg$start[tata],sgootarg$end[tata],sgootarg$end[tata],sgootarg$start[tata]), y=c(yli[1]-0.5,yli[1]-0.5, 0, 0),col=makeTransparent("black", alpha=100), border = NA)
                    }
                    # draw the final sample of loci
                    if (contig%in%names(sampling))
                    {
                        i.sam=which(contig==names(sampling))
                        sam.sta=sampling[[i.sam]]$start
                        sam.end=sampling[[i.sam]]$end
                        for (tata in 1:length(sam.sta))
                        {
                            polygon(x=c(sam.sta[tata],sam.end[tata],sam.end[tata],sam.sta[tata]), y=c(yli[1]-0.5,yli[1]-0.5, 0, 0),col=makeTransparent("darkgoldenrod2", alpha=100), border = NA)
                        }
                    }
                    
                }
                # plot the actual representation of alphas
                points(x, y, col=colos[icolo], type=typ, pch=pp[indiv], lwd=0.9)
            }
        }
        print(ite)
        print(contig)
        inc=inc+1
        
        # draw the legend if last box
        if (inc==(nr*nc)) 
        {
            plot(1,1, xaxt="n", yaxt="n") 
            legend("center", legend=c(race, "Medicago: reference"), col = colos, lty = 1, bg = 'gray75', cex=2.5)
            inc=1
        } 
        
    }
}

    library(VennDiagram)
    library(grid)
    library(gridBase)
    library(lattice)

    # data
    l1 <- list(Deletion=1:1420, Insertion=967:2042)
    l2 <- list(Deletion=1:502, Insertion=324:660)
    l3 <- list(Deletion=1:142, Insertion=85:184)
    l4 <- list(Deletion=1:161, Insertion=22:217)
    venns <- list(Subtargets=l1, Targets=l2, Genes=l3, Promoters=l4)

    
    gl <- grid.layout(nrow=2, ncol=2)

    # setup viewports
    vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
    vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 
    vp.3 <- viewport(layout.pos.col=1, layout.pos.row=2) 
    vp.4 <- viewport(layout.pos.col=2, layout.pos.row=2) 

    # init layout
    pushViewport(viewport(layout=gl))
    
    
    for (i in 1:4){
        # access the relevant viewport
        vp <- paste("vp.", i, sep="")
        pushViewport(get(vp))

        # draw the venn diagram
        temp <- venn.diagram(venns[[i]], fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 1,cat.fontface = 2, lty =2, filename = NULL, sub=names(venns)[i], margin = 0.2, sub.pos = c(0.5, 0.78), sub.col="blue")

        # start new base graphics in first viewport
        grid.draw(temp)

        # done with this viewport
        popViewport()
    }

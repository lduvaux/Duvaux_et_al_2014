get_elements <- function(x, sep="_", what=c(1)){
    return(strsplit(x, sep)[[1]][what])
}

collapse_elements <- function(nom, sep="_", what=1:3, colla="_"){
    res <- paste(get_elements(nom, sep, what), collapse=colla)
    return(res)
}

caract_bait <- function(vect_bait){
    tt <- table(vect_bait)
    vtest <- as.numeric(names(tt))
    
    if (length(which(vtest<1))>0) res1 <- 1
    else res1 <- 0
    
    if (length(which(vtest>1))>0) res2 <- 1
    else res2 <- 0

    if (res1==1 & res2==1) res3 <- 1
    else res3 <- 0
    
    return(c(res1, res2, res3))}

caract_bait_race <- function(race, noms, vect_bait, x){
    ind <- grep(race, noms)
    res <- caract_bait(vect_bait[ind])
    if (x%%200==0) print(x)
    return(res)}

caract_targ <- function(tab_targ){

    tt <- table(as.vector(tab_targ))
    vtest <- as.numeric(names(tt))
    
    if (length(which(vtest<1))>0) res1 <- 1
    else res1 <- 0
    
    if (length(which(vtest>1))>0) res2 <- 1
    else res2 <- 0

    if (res1==1 & res2==1) res3 <- 1
    else res3 <- 0
    
    return(c(res1, res2, res3))}

get_caract_targ <- function(target, subtargets, alpha_matrix){
    ind <- grep(paste(target, "_", sep=""), subtargets, fixed=T)
    tab_targ <- alpha_matrix[ind, ]
    res <- caract_targ(tab_targ)
    return(res)}




draw_venn <- function(){
    # libraries
    library(VennDiagram)
    library(grid)
    library(gridBase)
    library(lattice)

    # create the diagrams
    temp1 <- venn.diagram(list(B = 1:1800, A = 1571:2020),
        fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 1,cat.fontface = 2,
        lty =2, filename = NULL)
    temp2 <- venn.diagram(list(A = 1:1800, B = 1571:2020),
        fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 1,cat.fontface = 2,
        lty =2, filename = NULL)    




    plo <- venn.diagram(list4venn,fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3, filename = "Pairwise_Venn_diagram.tiff")
    grid.draw(plo)



    


    # start new page
    plot.new() 

    pdf("testpdf", width = 14, height = 7)
    # setup layout
    gl <- grid.layout(nrow=1, ncol=2)
    # grid.show.layout(gl)

    # setup viewports
    vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
    vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 

    # init layout
    pushViewport(viewport(layout=gl))
    # access the first position
    pushViewport(vp.1)

    # start new base graphics in first viewport
    par(new=TRUE, fig=gridFIG())

    grid.draw(temp2)

    # done with the first viewport
    popViewport()

    # move to the next viewport
    pushViewport(vp.2)

      grid.draw(temp2)

    # done with this viewport
    popViewport(1)

    dev.off()
}

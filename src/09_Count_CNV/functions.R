get_elements <- function(x, sep="_", what=c(1)){
    return(strsplit(x, sep)[[1]][what])
}

collapse_elements <- function(nom, sep="_", what=1:3, colla="_"){
    res <- paste(get_elements(nom, sep, what), collapse=colla)
    return(res)
}

caract_bait <- function(vect_bait){
# vect_bait can be a vector (bait) or a table (target or gene)
    tt <- table(as.vector(vect_bait))
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
    res <- caract_bait(tab_targ)
    return(res)}



caract_bait_race <- function(race, noms, vect_bait, x){
    ind <- grep(race, noms)
    res <- caract_bait(vect_bait[ind])
    if (x%%200==0) print(x)
    return(res)}


draw_venn <- function(venns){

#~    plot.new() 
    # setup layout
    gl <- grid.layout(nrow=2, ncol=2)
    # grid.show.layout(gl)

    # setup viewports
    vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
    vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 
    vp.3 <- viewport(layout.pos.col=1, layout.pos.row=2) 
    vp.4 <- viewport(layout.pos.col=2, layout.pos.row=2) 

    # init layout
    pushViewport(viewport(layout=gl))
    
    # access the first position
    for (i in 1:4){
        vp <- paste("vp.", i, sep="")
        pushViewport(get(vp))

        coco <- c("red", "blue")
        temp <- venn.diagram(venns[[i]], fill = coco, alpha = c(0.5, 0.5), cex = 1,cat.fontface = 2, lty =2, filename = NULL, sub=names(venns)[i], margin = 0.2, sub.pos = c(0.5, 0.78), sub.col="black", sub.fontface = "bold", sub.cex = 1.5, cat.col =coco)

        # start new base graphics in first viewport
#~        if (i==1) par(new=TRUE, fig=gridFIG())
        grid.draw(temp)

        # done with the first viewport
        popViewport()
    }
}

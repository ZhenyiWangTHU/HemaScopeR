#' Run copykat
#'
#' @param st_obj The Seurat object
#' @param save_path The path to save results
#' @param assay The assay to use
#' @param ... Parameters of `copykat::copykat`
#'
#' @import copykat
#' @import ggplot2
#' @import Seurat
#'
#' @export
#'
st_CNV <- function(
    st_obj,
    save_path,
    assay = 'Spatial',
    LOW.DR = 0.05,
    UP.DR = 0.1,
    win.size = 25,
    distance = "euclidean",
    genome = "hg20",
    n.cores = 1,
    species = 'human'
){

    if(is.null(genome)){
        if(species == 'human'){
            genome = 'hg20'
        }else if(species == 'mouse'){
            genome = 'mm10'
        }else{
            stop(paste0(species), ' is not supported in copykat.')
        }
    }

    if(!dir.exists(save_path)){
        dir.create(save_path)
    }
    setwd(save_path)

    if(packageVersion('SeuratObject') >= '5.0.0'){
        rawmat <- GetAssayData(st_obj, layer = 'counts',
                               assay = assay)
    }else{
        rawmat <- GetAssayData(st_obj, slot = 'counts',
                               assay = assay)
    }

    copykat.test <- copykat::copykat(
        rawmat = as.matrix(rawmat),
        LOW.DR = LOW.DR,
        UP.DR = UP.DR,
        win.size = win.size,
        distance = distance,
        genome = genome,
        n.cores = n.cores,
        plot.genes = 'FALSE')

    saveRDS(copykat.test, file.path(save_path, 'copykat_result.rds'))

    if(!dir.exists(file.path(save_path, 'pdf'))){
        dir.create(file.path(save_path, 'pdf'))
    }
    if(!dir.exists(file.path(save_path, 'png'))){
        dir.create(file.path(save_path, 'png'))
    }
    copykatPlot(copykat.test = copykat.test,
                save_path = save_path)

    pred.test <- data.frame(copykat.test$prediction)
    pred.test <- pred.test[make.names(colnames(st_obj)), 'copykat.pred']
    st_obj <- AddMetaData(st_obj, pred.test, col.name = 'CNV_state')
    st_obj$CNV_state[st_obj$CNV_state == 'aneuploid'] <- 'Aneuploid'
    st_obj$CNV_state[st_obj$CNV_state == 'diploid'] <- 'Diploid'
    suppressMessages(suppressWarnings(
        p <- SpatialDimPlot(st_obj, 'CNV_state', stroke = NA) +
            scale_fill_manual(name = 'CNV state',
                              values = list('Aneuploid' = rgb(254,129,125,maxColorValue = 255),
                                            'Diploid' = rgb(129,184,223,maxColorValue = 255),
                                            'not.defined' = 'grey50')) +
            theme(legend.key = element_blank()) +
            labs(name = 'CNV state') +
            guides(fill=guide_legend(override.aes = list(size=4)))
    ))
    saveImage(save_path,
              p,
              'CNV_Spatial',
              height = 4,
              width = 4)

    return(st_obj)
}


#' copykat plot
#'
#' @param copykat.test The results of `copykat::copykat`
#' @param save_path The path to save figures
#'
#' @import copykat
#' @import RColorBrewer
#' @import parallelDist
#'
copykatPlot <- function(
    copykat.test,
    save_path
){

    pred.test <- data.frame(copykat.test$prediction)
    rownames(pred.test) <- make.names(rownames(pred.test))
    if(length(which(pred.test$copykat.pred=="not.defined")) > 0){
        pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]
    }
    CNA.test <- data.frame(copykat.test$CNAmat)

    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

    chr <- as.numeric(CNA.test$chrom) %% 2+1
    rbPal1 <- colorRampPalette(c('black','grey'))
    CHR <- rbPal1(2)[as.numeric(chr)]
    chr1 <- cbind(CHR,CHR)

    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    com.preN <- pred.test$copykat.pred
    pred <- rbPal5(2)[as.numeric(factor(com.preN))]

    cells <- rbind(pred,pred)
    col_breaks = c(seq(-1,-0.4,length=50),
                   seq(-0.4,-0.2,length=150),
                   seq(-0.2,0.2,length=600),
                   seq(0.2,0.4,length=150),
                   seq(0.4, 1,length=50))

    pdf(file.path(save_path, 'pdf', 'CNV_heatmap.pdf'),
        width = 10, height = 10)
    heatmap.3(t(CNA.test[, rownames(pred.test)]),dendrogram="r",
              distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"),
              hclustfun = function(x) hclust(x, method="ward.D2"),
              ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              keysize=1, density.info="none", trace="none",
              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(5,5))
    # heatmap.3(t(CNA.test[,8:ncol(CNA.test)]),dendrogram="r",
    #           distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"),
    #           hclustfun = function(x) hclust(x, method="ward.D2"),
    #           ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
    #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
    #           keysize=1, density.info="none", trace="none",
    #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
    #           symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(5,5))
    legend("topright", paste("pred.",names(table(com.preN)),sep=""),
           pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
    dev.off()

    png(file.path(save_path, 'png', 'CNV_heatmap.png'),
        width = 1000, height = 1000)
    heatmap.3(t(CNA.test[, rownames(pred.test)]),dendrogram="r",
              distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"),
              hclustfun = function(x) hclust(x, method="ward.D2"),
              ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              keysize=1, density.info="none", trace="none",
              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(5,5))
    # heatmap.3(t(CNA.test[,8:ncol(CNA.test)]),dendrogram="r",
    #           distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"),
    #           hclustfun = function(x) hclust(x, method="ward.D2"),
    #           ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
    #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
    #           keysize=1, density.info="none", trace="none",
    #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
    #           symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(5,5))
    legend("topright", paste("pred.",names(table(com.preN)),sep=""),
           pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
    dev.off()


    tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
    tumor.cells <- make.names(tumor.cells)
    tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
    hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
    hc.umap <- cutree(hcc,2)

    rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
    subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
    cells <- rbind(subpop,subpop)

    pdf(file.path(save_path, 'pdf', 'Tumor_subpopulations.pdf'),
        width = 10, height = 10)
    heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
              ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              keysize=1, density.info="none", trace="none",
              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(5,5))

    legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')
    dev.off()

    png(file.path(save_path, 'png', 'Tumor_subpopulations.png'),
        width = 1000, height = 1000)
    heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
              ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              keysize=1, density.info="none", trace="none",
              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(5,5))

    legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')
    dev.off()
}

#' Run copykat for scRNA-seq data
#'
#' @param sc_object The Seurat object
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
sc_CNV = function(
    sc_object=NULL,
    save_path=NULL,
    assay = 'RNA',
    LOW.DR = 0.05,
    UP.DR = 0.1,
    win.size = 25,
    distance = "euclidean",
    genome = "hg20",
    n.cores = 1,
    species = 'hsa'
){
    require(copykat)

    if(is.null(genome)){
        if(species == 'hsa'){
            genome = 'hg20'
        }else if(species == 'mmu'){
            genome = 'mm10'
        }else{
            stop(paste0(species), ' is not supported in copykat.')
        }
    }

    if(!dir.exists(save_path)){
        dir.create(save_path)
    }
    
    temp.wd <- getwd()
    setwd(save_path)

    copykat.test <- copykat::copykat(
        rawmat = as.matrix(GetAssayData(sc_object,
                                        slot = 'counts',
                                        assay = assay)),
        LOW.DR = LOW.DR,
        UP.DR = UP.DR,
        win.size = win.size,
        distance = distance,
        genome = genome,
        n.cores = n.cores)

    saveRDS(copykat.test, file.path(save_path, 'copykat_result.rds'))
    pred.test <- data.frame(copykat.test$prediction)
    sc_object@meta.data$CNV_state <- mapvalues(rownames(sc_object@meta.data),
                                               from=pred.test$cell.names,
                                               to=pred.test$copykat.pred)
    
    pdf('./CNV_state.pdf', width = 7, height = 6)
      print(DimPlot(sc_object, group.by = "CNV_state", label = T, pt.size = 1,label.size = 8) + scale_color_manual(values = c("red", "grey","blue")))
    dev.off()

    png('./CNV_state.png', width = 700, height = 600)
      print(DimPlot(sc_object, group.by = "CNV_state", label = T, pt.size = 1,label.size = 8) + scale_color_manual(values = c("red", "grey","blue")))
    dev.off()

    setwd(temp.wd)
    return(sc_object)
}
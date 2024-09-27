#' Loading data
#'
#' @param input.data.dirs A character of data path, where there are filtered_feature_bc_matrix.h5
#' and spatial folder
#' @param output.dir A character of path to store the h5ad file
#' @param sampleName A character of the name of the sample
#' @param rds.file description
#'
#'
#' @details
#' This function loads data and generates corresponding h5ad file
#'
#' @import Seurat
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
st_Loading_Data <- function(
        input.data.dir,
        output.dir = '.',
        sampleName = 'Hema_ST',

        rds.file = FALSE,
        filename = "filtered_feature_bc_matrix.h5",
        assay = "Spatial",
        slice = "slice1",
        filter.matrix = TRUE,
        to.upper = FALSE
){
    if(!dir.exists(output.dir)){
        dir.create(output.dir)
    }

    if(rds.file){
        st_obj <- readRDS(input.data.dir)

        if(!(assay %in% names(st_obj@assays))){
            stop(paste0(assay, ' is not in the assays of the object.'))
        }
        if(!(slice %in% names(st_obj@images))){
            stop(paste0(slice, ' is not in the images of the object.'))
        }

        st_obj@active.assay <- assay
    }else{
        st_obj <- Load10X_Spatial(input.data.dir,
                                  filename = filename,
                                  assay = assay,
                                  slice = slice,
                                  filter.matrix = filter.matrix,
                                  to.upper = to.upper,
                                  image = NULL)
    }
    st_obj@project.name <- sampleName
    st_obj$orig.ident <- factor(sampleName)
    st_obj@active.ident <- st_obj$orig.ident

    # suppressMessages(suppressWarnings(
    #     SaveH5Seurat(st_obj,
    #                  filename = file.path(output.dir, 'st_object.h5Seurat'),
    #                  overwrite = TRUE)
    # ))
    # suppressMessages(suppressWarnings(
    #     Convert(file.path(output.dir, 'st_object.h5Seurat'), dest = "h5ad",
    #             overwrite = TRUE)
    # ))

    if(rds.file){
        Rds2H5(object = st_obj,
               file.dir = output.dir,
               assay = assay,
               slice = slice)
    }else{
        tmp <- file.copy(from = file.path(input.data.dir, filename),
                         to = file.path(output.dir, filename),
                         overwrite = TRUE)
        tmp <- file.copy(from = file.path(input.data.dir, 'spatial'),
                         to = output.dir,
                         recursive = TRUE,
                         overwrite = TRUE)
    }

    return(st_obj)
}

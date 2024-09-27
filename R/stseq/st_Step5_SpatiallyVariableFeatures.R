#' Find spatially variable features
#'
#' @param st_obj The SeuratObject of 10X visium data
#' @param output.dir A character of path to store the results and figures
#' @param assay The assay used
#' @param selection.method Method for selecting spatially variable features,
#' `markvariogram` or `moransi`
#' @param n.top.show The number of top genes shown in the figure
#' @param n.col The number of columns in the output figure
#' @param verbose verbose as `Seurat`
#'
#' @details
#' This function finds spatially variable features and plots the distribution of
#' these genes
#'
#'
#' @import Seurat
#' @import Rfast2
#'
#' @export
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'

st_SpatiallyVariableFeatures <- function(
    st_obj,
    output.dir = '.',
    assay = 'SCT',
    selection.method = 'moransi',
    n.top.show = 10,
    n.col = 5,
    verbose = FALSE
){
    if(!dir.exists(output.dir)){
        dir.create(output.dir)
    }

    suppressWarnings(
        st_obj <- FindSpatiallyVariableFeatures(object = st_obj,
                                                assay = assay,
                                                selection.method = selection.method,
                                                verbose = verbose)
    )
    meta.features <- st_obj@assays[[assay]]@meta.features
    meta.features <- meta.features[!is.na(meta.features), ]
    write.csv(meta.features,
              file.path(output.dir, 'SVF.csv'),
              row.names = T)

    svfs <- mySpatiallyVariableFeatures(object = st_obj,
                                        assay = assay,
                                        selection.method = selection.method)
    p.svfs <- mySpatialFeaturePlot(st_obj = st_obj,
                                   features = svfs[1:n.top.show],
                                   n.col = n.col)
    saveImage(output.dir,
              p.svfs,
              'SpatiallyVariableFeatures',
              height = 2.5*ceiling(length(p.svfs) / n.col),
              width = 3*min(length(p.svfs), n.col))

    return(st_obj)
}


#' SpatiallyVariableFeatures revision
#'
#'
mySpatiallyVariableFeatures <- function(
    object,
    assay = "SCT",
    selection.method = "moransi"
){
    meta.features <- object@assays[[assay]]@meta.features
    moransi_cols <- grep(paste0("^", selection.method), colnames(meta.features), value = TRUE)
    filtered_data <- meta.features[meta.features[[paste0(selection.method, ".spatially.variable")]], moransi_cols]
    sorted_data <- filtered_data[order(filtered_data[[paste0(selection.method, ".spatially.variable.rank")]]), ]
    return(rownames(sorted_data))
}



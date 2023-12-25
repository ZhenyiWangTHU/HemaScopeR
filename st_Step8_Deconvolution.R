#' The pipeline for integrated analysis of 10x visium data and scRNA-seq data
#'
#' @param st.data.dir A character of data path, where there are filtered_feature_bc_matrix.h5
#' and spatial folder
#' @param sc.h5ad.dir A character of path of h5ad file of scRNA-seq reference data
#' @param save.path A character of path to save data
#' @param sc.labels.key A character indicating the key of scRNA-seq data labels
#' used during running cell2location
#' @param species A character, `human` or `mouse`
#' @param condaenv A character of the name of conda environment used here
#' @param use.gpu A bool value indicating whether to use GPU
#' @param sc.max.epoch A integer representing the maximum epochs of training scRNA-seq data
#' @param st.max.epoch A integer representing the maximum epochs of training ST data
#' @param use.Dataset 'HematoMap' or 'LymphNode'
#'
#' @import reticulate
#' @import Seurat
#' @import SeuratDisk
#' @import RColorBrewer
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
st_Deconvolution <- function(
    st.data.dir,
    sc.h5ad.dir = NULL,
    library_id = 'Hema_ST',
    st_obj = NULL,
    save_path = NULL,
    sc.labels.key = 'seurat_clusters',
    species = 'mouse',
    sc.max.epoch = 1000,
    st.max.epoch = 10000,
    use.gpu = TRUE,
    condaenv = 'r-reticulate',
    use.Dataset = 'LymphNode'
){
    use_condaenv(condaenv)
    source_python(file.path(system.file(package = "HemaScopeR"),
                            "python/cell2loc.py"))
    if(is.null(save_path)){
        save_path = getwd()
    }
    if(!dir.exists(save_path)){
        dir.create(save_path)
    }
    if(!dir.exists(file.path(save_path, 'pdf'))){
        dir.create(file.path(save_path, 'pdf'))
    }
    if(!dir.exists(file.path(save_path, 'png'))){
        dir.create(file.path(save_path, 'png'))
    }

    ### To do
    if(is.null(sc.h5ad.dir)){
        if(use.Dataset == 'HematoMap'){
            default_path <- file.path(system.file(package = 'HemaScopeR'),
                                      ifelse(species == 'mouse',
                                             yes = 'extdata/HematoMap_mouse.csv',
                                             no = 'extdata/HematoMap_human.csv'))
        }else if(use.Dataset == 'LymphNode'){
            default_path <- file.path(system.file(package = 'HemaScopeR'),
                                      ifelse(species == 'mouse',
                                             yes = 'extdata/LymphNode_mouse.csv',
                                             no = 'extdata/LymphNode_human.csv'))
        }else{
            stop(paste0('No dataset called ', use.Dataset))
        }
    }

    run_cell2loc(st_data_path = st.data.dir,
                 sc_h5ad_path = sc.h5ad.dir,
                 library_id = library_id,
                 labels_key = sc.labels.key,
                 save_path = save_path,
                 sc_max_epochs = as.integer(sc.max.epoch),
                 st_max_epochs = as.integer(st.max.epoch),
                 use_gpu = as.logical(use.gpu),
                 species = species,
                 default_path = default_path)

    # suppressMessages(suppressWarnings(
    #     Convert(file.path(save_path, 'st.h5ad'),
    #             dest = 'h5seurat',
    #             assay = 'Spatial',
    #             overwrite = TRUE)
    # ))
    # suppressMessages(suppressWarnings(
    #     st_obj.cell2loc <- LoadH5Seurat(file.path(save_path, 'st.h5seurat'),
    #                                     assay = 'Spatial')
    # ))
    # # saveRDS(st_obj.cell2loc,
    # #         file.path(save_path, 'st_cell2loc.rds'))
    #
    # labels.name <- colnames(st_obj.cell2loc@meta.data)
    # labels.name <- labels.name[7:length(labels.name)]
    # cell2loc.meta <- st_obj.cell2loc[[labels.name]]
    # write.csv(cell2loc.meta,
    #           file.path(save_path, 'cell2loc_res.csv'),
    #           row.names = TRUE)
    cell2loc.meta <- read.csv(file.path(save_path, 'cell2loc_res.csv'),
                              row.names = 1)
    labels.name <- colnames(cell2loc.meta)

    if(is.null(st_obj)){
        st_obj.cell2loc <- Load10X_Spatial(st.data.dir)
    }else{
        st_obj.cell2loc <- st_obj
    }
    st_obj.cell2loc <- st_obj.cell2loc[, rownames(cell2loc.meta)]
    st_obj.cell2loc <- AddMetaData(st_obj.cell2loc, cell2loc.meta)

    nCount.key <- grep('^nCount_', colnames(st_obj.cell2loc@meta.data), value = TRUE)[1]
    if(is.na(nCount.key)){
        st_obj.cell2loc <- AddMetaData(st_obj.cell2loc,
                                       colSums(st_obj.cell2loc@assays[[st_obj.cell2loc@active.assay]]@counts),
                                       col.name = 'nCount')
        nCount.key <- 'nCount'
    }

    # absolute abundance
    FeatureColors.one <- colorRampPalette(colors = brewer.pal(n = 9, name = 'YlOrRd'))
    for(labels in labels.name){
        p.label.spatial <- mySpatialFeaturePlot(st_obj = st_obj.cell2loc,
                                                features = labels,
                                                legend.name = labels,
                                                legend.color = FeatureColors.one)
        saveImage(save_path,
                  p.label.spatial,
                  labels,
                  height = 4,
                  width = 4)
    }

    p.labels.spatial <- mySpatialFeaturePlot(st_obj = st_obj.cell2loc,
                                             features = labels.name,
                                             n.col = 4)
    saveImage(save_path,
              p.labels.spatial,
              'Cell2loc_total',
              height = 2.5*ceiling(length(p.labels.spatial) / 4),
              width = 3*min(length(p.labels.spatial), 4))

    # relative abundance
    # for(labels in labels.name){
    #     st_obj.cell2loc <- AddMetaData(st_obj.cell2loc,
    #                                    metadata = st_obj.cell2loc[[labels]] /
    #                                        st_obj.cell2loc[[nCount.key]] *
    #                                        median(st_obj.cell2loc[[nCount.key]][, 1]),
    #                                    col.name = paste0(labels, '_normalized'))
    #     p.label.spatial <- mySpatialFeaturePlot(st_obj = st_obj.cell2loc,
    #                                             features = paste0(labels, '_normalized'),
    #                                             legend.name = paste0(labels, '_norm'),
    #                                             legend.color = FeatureColors.one)
    #     saveImage(save_path,
    #               p.label.spatial,
    #               paste0(labels, '_normalized'),
    #               height = 4,
    #               width = 4)
    # }
    # p.labels.spatial <- mySpatialFeaturePlot(st_obj = st_obj.cell2loc,
    #                                          features = paste0(labels.name, '_normalized'),
    #                                          legend.name = paste0(labels.name, '_norm'),
    #                                          n.col = 4)
    # saveImage(save_path,
    #           p.labels.spatial,
    #           'Cell2loc_total_norm',
    #           height = 2.5*ceiling(length(p.labels.spatial) / 4),
    #           width = 3.2*min(length(p.labels.spatial), 4))

    saveRDS(st_obj.cell2loc,
            file.path(save_path, 'st_cell2loc.rds'))

    if(!is.null(st_obj)){
        return(st_obj.cell2loc)
    }
}

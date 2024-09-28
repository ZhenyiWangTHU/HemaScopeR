#' Normalization, PCA, clustering and visualization
#'
#' @param st_obj The SeuratObject of 10X visium data
#' @param output.dir A character of path to store the results and figures
#' @param normalization.method `SCTransform` or other choice of `normalization.method` in `NormalizaData`
#' @param npcs The number of dimensions for analysis
#' @param pcs.used The principal components (PCs) to use for analysis
#' @param resolution An float parameter of `FindClusters` of `Seurat`
#' @param verbose verbose as `Seurat`
#'
#' @details
#' This function performs normalization, PCA, clustering and visualization on the
#' SeuratObject of 10X visium data and returns the processed SeuratObject
#'
#' @import Seurat
#' @import ggplot2
#' @import RColorBrewer
#' @import dplyr
#'
#' @export
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
st_Clustering <- function(
    st_obj,
    output.dir = '.',
    assay = 'Spatial',
    normalization.method = 'SCTransform',
    npcs = 50,
    pcs.used = 1:10,
    resolution = 0.8,
    max.n.cluster = 30,
    verbose = FALSE
){
    if(!dir.exists(output.dir)){
        dir.create(output.dir)
    }

    if(normalization.method == 'SCTransform'){
        st_obj <- st_obj %>% SCTransform(assay = assay, verbose = verbose) %>%
            RunPCA(assay = 'SCT', npcs = npcs, verbose = verbose)
    }else{
        st_obj <- st_obj %>% NormalizeData(normalization.method = normalization.method, verbose = verbose) %>%
            FindVariableFeatures(verbose = verbose) %>%
            ScaleData(verbose = verbose) %>%
            RunPCA(assay = assay, npcs = npcs, verbose = verbose)
    }

    suppressMessages(suppressWarnings(
        st_obj <- st_obj %>% FindNeighbors(reduction = "pca", dims = pcs.used, verbose = verbose) %>%
            FindClusters(resolution = resolution, verbose = verbose) %>%
            RunUMAP(reduction = "pca", dims = pcs.used, verbose = verbose) %>%
            RunTSNE(reduction = "pca", dims = pcs.used, verbose = verbose)
    ))
    n.cluster <- length(unique(st_obj$seurat_clusters))

    while (n.cluster > max.n.cluster) {
        print(paste0('Now the number of clusters is ', n.cluster, ', which exceeds max.n.cluster. Trying less resolution...'))
        resolution = resolution / 2
        st_obj <- st_obj %>% FindClusters(resolution = resolution, verbose = verbose)
        if(n.cluster == length(unique(st_obj$seurat_clusters))){
            print('The number of clusters will not change.')
            break
        }
        n.cluster <- length(unique(st_obj$seurat_clusters))
    }

    ##### Visualization #####
    suppressMessages(suppressWarnings(
        p.cluster.spatial <- SpatialDimPlot(st_obj, group.by = 'seurat_clusters', stroke = NA) +
            scale_fill_manual(name = 'Cluster',
                              values = getDefaultClusterColor(n.cluster)) +
            theme(legend.position = 'right',
                  legend.key = element_blank()) +
            guides(fill=guide_legend(override.aes = list(size=4)))
    ))
    saveImage(output.dir,
              p.cluster.spatial,
              'cluster_spatial',
              height = 4,
              width = 4.5)

    p.cluster.umap <- DimPlot(st_obj, reduction = 'umap') +
        scale_color_manual(name = 'Cluster',
                           values = getDefaultClusterColor(n.cluster)) +
        theme(legend.position = 'right',
              legend.key = element_blank()) +
        guides(color=guide_legend(override.aes = list(size=4)))
    saveImage(output.dir,
              p.cluster.umap,
              'cluster_UMAP',
              height = 4,
              width = 5)

    p.cluster.tsne <- DimPlot(st_obj, reduction = 'tsne') +
        scale_color_manual(name = 'Cluster',
                           values = getDefaultClusterColor(n.cluster)) +
        theme(legend.position = 'right',
              legend.key = element_blank()) +
        guides(color=guide_legend(override.aes = list(size=4)))
    saveImage(output.dir,
              p.cluster.tsne,
              'cluster_tSNE',
              height = 4,
              width = 5)

    ##### Save the data #####
    saveRDS(st_obj, file.path(output.dir, 'object_Clustering.rds'))

    return(st_obj)
}


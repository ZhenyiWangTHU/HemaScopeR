#' Find DEGs of clusters
#'
#' @param st_obj The SeuratObject of 10X visium data
#' @param output.dir A character of path to store the results and figures
#' @param ident.label A character of the label set to `active.ident`
#' @param only.pos A bool value to indicate whether only return positive markers
#' @param min.pct The parameter of `FindAllMarkers`
#' @param logfc.threshold The parameter of `FindAllMarkers`
#' @param test.use The test used in `FindAllMarkers`
#' @param verbose verbose as `Seurat`
#'
#' @details
#' This function finds and visualizes markers for each of the identity classes and
#' returns the marker dataframe.
#'
#' @import Seurat
#' @import dplyr
#' @import ggplot2
#' @import SeuratObject
#'
#' @export
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
st_Find_DEGs <- function(
        st_obj,
        output.dir = '.',
        ident.label = 'seurat_clusters',
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25,
        test.use = 'wilcox',
        verbose = FALSE
){
    if(!dir.exists(output.dir)){
        dir.create(output.dir)
    }

    st_obj <- SeuratObject::SetIdent(st_obj, value = st_obj@meta.data[, ident.label])
    n.cluster <- length(unique(st_obj@active.ident))
    st_obj.markers <- FindAllMarkers(st_obj,
                                     only.pos = only.pos,
                                     min.pct = min.pct,
                                     logfc.threshold = logfc.threshold,
                                     test.use = test.use,
                                     verbose = verbose)
    st_obj.markers.top5 <- st_obj.markers %>%
        group_by(cluster) %>%
        top_n(n = 5, wt = .data[[grep('avg_log', colnames(st_obj.markers), value = T)]])
    st_obj.markers.top5 <- st_obj.markers.top5[!duplicated(st_obj.markers.top5$gene), ]
    p.degs.dot <- DotPlot(st_obj,
                          features = st_obj.markers.top5$gene,
                          cols=c("lightgrey",'red'),
                          group.by = ident.label) +
        theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5)) +
        scale_y_discrete(limits = rev(levels(st_obj@meta.data[, ident.label])))
    saveImage(output.dir,
              p.degs.dot,
              'DEGs_dot',
              height = 3+0.2*n.cluster,
              width = 2+0.25*nrow(st_obj.markers.top5))

    write.csv(st_obj.markers,
              file.path(output.dir, 'markers.csv'),
              row.names = T)

    return(st_obj.markers)
}

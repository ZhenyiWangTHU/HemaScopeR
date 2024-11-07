#' Dimension reduction by HematoMap.
#'
#' Function \code{run_HematoMap} performs dimension reduction using HematoMap.
#'
#' @param sc_object A Seurat object containing single-cell RNA-seq data.
#' @param output.dir The directory where cluster tree plots will be saved.
#'
#' @details
#' This function performs dimension reduction using HematoMap on a Seurat object.
#' It generates cluster tree plots and saves them as PDF files and PNG files in the specified output directory.
#' 
#' @return run_HematoMap returns cluster tree plots in both PDF and PNG formats.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#'
#' @export
run_HematoMap = function(sc_object = NULL,
                         output.dir = NULL){
    # Construct HematoMap object
    hmap <- HematoMap::CreateSubclusterObject(sc_object)
    saveRDS(hmap, file = paste0(gsub("/Step5.Visualization/$", "/RDSfiles", output.dir),'/hmap.RDS'))

    # Visualization for input data (hierarchy tree)
    # pdf
    pdf(paste0(output.dir, '/hierarchy tree of input data for mapping cell types.pdf'), width=8, height=8)
    print(plotClusterTree(hemato.subc = hmap, group.subc = "external",
                    color.mapping = "cell.type", title = paste('The input', nrow(sc_object@meta.data), 'cells'),
                    point.size = 20, label.cell = TRUE))
    dev.off()

    pdf(paste0(output.dir, '/hierarchy tree of input data for mapping proportions.pdf'), width=8, height=8)
    print(plotClusterTree(hemato.subc = hmap, group.subc = "external",
                    color.mapping = "cell.percentage", title = paste('The input', nrow(sc_object@meta.data), 'cells'),
                    point.size = 20, label.cell = TRUE))
    dev.off()

    # Visualization for the normal BMMCs (hierarchy tree)
    pdf(paste0(output.dir, '/hierarchy tree of normal BMMCs for mapping cell types.pdf'), width=8, height=8)
    print(plotClusterTree(hemato.subc = hmap, group.subc = "reference",
                    color.mapping = "cell.type", title = "Normal BMMCs",
                    point.size = 20, label.cell = TRUE))
    dev.off()

    pdf(paste0(output.dir, '/hierarchy tree of normal BMMCs for mapping proportions.pdf'), width=8, height=8)
    print(plotClusterTree(hemato.subc = hmap, group.subc = "reference",
                    color.mapping = "cell.percentage", title = "Normal BMMCs",
                    point.size = 20, label.cell = TRUE))
    dev.off()

    # png
    png(paste0(output.dir, '/hierarchy tree of input data for mapping cell types.png'), width=800, height=800)
    print(plotClusterTree(hemato.subc = hmap, group.subc = "external",
                    color.mapping = "cell.type", title = paste('The input', nrow(sc_object@meta.data), 'cells'),
                    point.size = 20, label.cell = TRUE))
    dev.off()

    png(paste0(output.dir, '/hierarchy tree of input data for mapping proportions.png'), width=800, height=800)
    print(plotClusterTree(hemato.subc = hmap, group.subc = "external",
                    color.mapping = "cell.percentage", title = paste('The input', nrow(sc_object@meta.data), 'cells'),
                    point.size = 20, label.cell = TRUE))
    dev.off()

    # Visualization for the normal BMMCs (hierarchy tree)
    png(paste0(output.dir, '/hierarchy tree of normal BMMCs for mapping cell types.png'), width=800, height=800)
    print(plotClusterTree(hemato.subc = hmap, group.subc = "reference",
                    color.mapping = "cell.type", title = "Normal BMMCs",
                    point.size = 20, label.cell = TRUE))
    dev.off()

    png(paste0(output.dir, '/hierarchy tree of normal BMMCs for mapping proportions.png'), width=800, height=800)
    print(plotClusterTree(hemato.subc = hmap, group.subc = "reference",
                    color.mapping = "cell.percentage", title = "Normal BMMCs",
                    point.size = 20, label.cell = TRUE))
    dev.off()

}
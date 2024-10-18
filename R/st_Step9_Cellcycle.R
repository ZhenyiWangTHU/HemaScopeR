#' Calculate cell cycle scores
#'
#' @param st_obj The Seurat object
#' @param save_path The path to save results
#' @param s.features A vector of features associated with S phase
#' @param g2m.features A vector of features associated with G2M phase
#' @param species A character, `human` or `mouse`
#' @param FeatureColors.bi A function that interpolates a set of given colors to create new color palettes and color ramps
#'
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @import RColorBrewer
#'
#' @export
#'
st_Cell_cycle <- function(
        st_obj,
        save_path = '.',
        s.features = NULL,
        g2m.features = NULL,
        species = 'human',
        FeatureColors.bi = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu')))
){
    if(!dir.exists(save_path)){
        dir.create(save_path)
    }
    s.features = unlist(
        ifelse(is.null(s.features),
               ifelse(species == 'human',
                      yes = list(Seurat::cc.genes.updated.2019$s.genes),
                      no = list(convertHumanGene(Seurat::cc.genes.updated.2019$s.genes))),
               list(s.features))
    )
    g2m.features = unlist(
        ifelse(is.null(g2m.features),
               ifelse(species == 'human',
                      yes = list(Seurat::cc.genes.updated.2019$g2m.genes),
                      no = list(convertHumanGene(Seurat::cc.genes.updated.2019$g2m.genes))),
               list(g2m.features))
    )

    st_obj <- CellCycleScoring(
        st_obj,
        s.features = s.features,
        g2m.features = g2m.features
    )

    p.S.Score.spatial <- mySpatialFeaturePlot(st_obj = st_obj,
                                              features = 'S.Score',
                                              legend.name = 'S.Score',
                                              legend.color = FeatureColors.bi)
    saveImage(save_path,
              p.S.Score.spatial,
              'S_Score_spatial',
              height = 4,
              width = 4)

    p.G2M.Score.spatial <- mySpatialFeaturePlot(
        st_obj = st_obj,
        features = 'G2M.Score',
        legend.name = 'G2M.Score',
        legend.color = FeatureColors.bi
    )
    saveImage(save_path,
              p.G2M.Score.spatial,
              'G2M_Score_spatial',
              height = 4,
              width = 4)

    if('umap' %in% names(st_obj@reductions)){
        suppressMessages(suppressWarnings(
            p.S.Score.UMAP <- FeaturePlot(
                st_obj, features = 'S.Score',
                reduction = 'umap'
            ) +
                theme(legend.position = 'right') +
                scale_color_gradientn(
                    name = 'S.Score',
                    colors = FeatureColors.bi(n = 100)
                )
        ))
        saveImage(save_path,
                  p.S.Score.UMAP,
                  'S_Score_UMAP',
                  height = 4,
                  width = 4)

        suppressMessages(suppressWarnings(
            p.G2M.Score.UMAP <- FeaturePlot(
                st_obj, features = 'G2M.Score',
                reduction = 'umap'
            ) +
                theme(legend.position = 'right') +
                scale_color_gradientn(
                    name = 'G2M.Score',
                    colors = FeatureColors.bi(n = 100)
                )
        ))
        saveImage(save_path,
                  p.G2M.Score.UMAP,
                  'G2M_Score_UMAP',
                  height = 4,
                  width = 4)
    }

    if('tsne' %in% names(st_obj@reductions)){
        suppressMessages(suppressWarnings(
            p.S.Score.TSNE <- FeaturePlot(
                st_obj, features = 'S.Score',
                reduction = 'tsne'
            ) +
                theme(legend.position = 'right') +
                scale_color_gradientn(
                    name = 'S.Score',
                    colors = FeatureColors.bi(n = 100)
                )
        ))
        saveImage(save_path,
                  p.S.Score.TSNE,
                  'S_Score_TSNE',
                  height = 4,
                  width = 4)

        suppressMessages(suppressWarnings(
            p.G2M.Score.TSNE <- FeaturePlot(
                st_obj, features = 'G2M.Score',
                reduction = 'tsne'
            ) +
                theme(legend.position = 'right') +
                scale_color_gradientn(
                    name = 'G2M.Score',
                    colors = FeatureColors.bi(n = 100)
                )
        ))
        saveImage(save_path,
                  p.G2M.Score.TSNE,
                  'G2M_Score_TSNE',
                  height = 4,
                  width = 4)
    }

    return(st_obj)
}

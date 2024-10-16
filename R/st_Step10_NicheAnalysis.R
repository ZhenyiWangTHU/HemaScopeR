#' Co occurence analysis
#'
#' @param st_obj SeuratObject
#' @param features A vector of features
#' @param save_path Save path
#' @param slice slice
#' @param p The power
#' @param method correlation or Wasserstein
#'
#' @import transport
#' @import Seurat
#' @import ggplot2
#' @import pheatmap
#'
#' @export
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
Co_occurenceScore <- function(
        st_obj,
        features,
        save_path = NULL,
        slice = 'slice1',
        method = 'Wasserstein',
        p = 2
){
    res <- data.frame(matrix(0,
                             nrow = length(features),
                             ncol = length(features)))
    colnames(res) <- features
    rownames(res) <- features

    for(i in features){
        if(i %in% names(st_obj@meta.data)){
            next
        }else if(i %in% rownames(GetAssayData(st_obj))){
            st_obj[[i]] <- GetAssayData(st_obj, slot = 'data')[i, ]
        }else{
            stop(paste0('Feature ', i, ' is not in the meta.data or genes of object.'))
        }
    }

    if(method == 'Wasserstein'){
        pos <- GetTissueCoordinates(st_obj[[slice]])

        for(i in 1:length(features)){
            a <- wpp(coordinates = pos,
                     mass = st_obj[[features[i]]][, 1] / sum(st_obj[[features[i]]]))
            for(j in i:length(features)){
                b <- wpp(coordinates = pos,
                         mass = st_obj[[features[j]]][, 1] / sum(st_obj[[features[j]]]))
                res[i, j] <- CalculateWasserstein(a, b, p)
            }
        }

        res <- res + t(res)

        res <- 1 / (1 + res)
        diag(res) <- NA

    }else if(method == 'correlation'){
        for(i in 1:length(features)){
            for(j in i:length(features)){
                res[i, j] <- cor(st_obj[[features[i]]][, 1],
                                 st_obj[[features[j]]][, 1],
                                 method = 'pearson')
            }
        }

        res <- res + t(res)
        diag(res) <- NA

    }else{
        stop(paste0('No method called ', method))
    }

    write.csv(res,
              file.path(save_path, 'Coexistence.csv'))

    p <- pheatmap(res,
                  main = 'Coexistence score')

    # # res <- res / rowSums(res)
    #
    # df <- data.frame(feature1 = matrix(0, nrow = ncol(res) * nrow(res)),
    #                  feature2 = matrix(0, nrow = ncol(res) * nrow(res)),
    #                  coIndex = matrix(0, nrow = ncol(res) * nrow(res)))
    #
    # ii = 1
    # for(i in 1:nrow(res)){
    #     for(j in 1:ncol(res)){
    #         df$feature1[ii] <- rownames(res)[i]
    #         df$feature2[ii] <- colnames(res)[j]
    #         df$coIndex[ii] <- res[[i, j]]
    #
    #         ii = ii + 1
    #     }
    # }
    #
    # boundary.levels <- sort(unique(df$feature1))
    # df$feature1 <- factor(df$feature1, levels = boundary.levels)
    # df$feature2 <- factor(df$feature2, levels = boundary.levels)
    #
    # p <- ggplot(df, aes(x = feature1, y = feature2, fill = coIndex)) +
    #     geom_raster() +
    #     theme_linedraw() + theme(panel.grid.major = element_blank()) +
    #     scale_y_discrete(limits = rev(levels(df$feature2))) +
    #     scale_fill_gradientn(colors = c('white', '#D73027'),
    #                          name = 'Score',
    #                          na.value = 'grey50') +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1,
    #                                      vjust = 0.5), axis.title.x = element_blank(),
    #           axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")

    if(!is.null(save_path)){
        saveImage(save_path,
                  p,
                  'Coexistence',
                  height = 1.5+0.25*length(features),
                  width = 1.5+0.25*length(features))
    }else{
        p
    }

    return(df)
}

#' Calculate the wasserstein distance
#'
#' @param a A wpp object
#' @param b A wpp object
#' @param p The power
#'
#' @import transport
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
CalculateWasserstein <- function(
        a,
        b,
        p = 2
){
    require(transport)
    a$mass <- a$mass / sum(a$mass)
    b$mass <- b$mass / sum(b$mass)
    return(wasserstein(a, b, p = p))
}

#' Cluster niches
#'
#' @param st_obj SeuratObject
#' @param features A vector of features
#' @param save_path Save path
#' @param kmeans.n The number of clusters
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
#' @export
#'
#' @import ggplot2
#' @import RColorBrewer
#'
NicheCluster <- function(
        st_obj,
        features,
        save_path = NULL,
        kmeans.n = 4
){
    mat <- data.frame(matrix(NA,
                             nrow = ncol(st_obj),
                             ncol = length(features)))
    rownames(mat) <- colnames(st_obj)
    colnames(mat) <- features
    for(i in features){
        if(i %in% names(st_obj@meta.data)){
            mat[, i] <- st_obj@meta.data[, i]
        }else if(i %in% rownames(GetAssayData(st_obj))){
            mat[, i] <- GetAssayData(st_obj, slot = 'data')[i, ]
        }else{
            stop(paste0('Feature ', i, ' is not in the meta.data or genes of object.'))
        }
    }
    res <- kmeans(mat, centers = kmeans.n)
    st_obj <- AddMetaData(st_obj,
                          res$cluster,
                          col.name = 'NicheCluster')
    suppressMessages(suppressWarnings(
        p.nich <- SpatialDimPlot(st_obj, group.by = 'NicheCluster', stroke = NA) +
            scale_fill_manual(name = 'NicheCluster',
                              values = getDefaultClusterColor(kmeans.n)) +
            theme(legend.position = 'right',
                  legend.key = element_blank()) +
            guides(fill=guide_legend(override.aes = list(size=4)))
    ))
    if(!is.null(save_path)){
        saveImage(save_path,
                  p.nich,
                  'niche_spatial',
                  height = 4,
                  width = 4)
    }else{
        p.nich
    }

    # nich cell type enrichment
    tmp <- data.frame(matrix(NA,
                             nrow = kmeans.n,
                             ncol = length(features)))
    rownames(tmp) <- c(1:kmeans.n)
    colnames(tmp) <- features
    for(i in features){
        i.tmp <- aggregate(mat[, i],
                           by = list(nich.cluster = res$cluster),
                           FUN = mean)
        rownames(i.tmp) <- i.tmp$nich.cluster
        i.tmp <- i.tmp[rownames(tmp), ]
        tmp[, i] <- i.tmp$x
    }

    tmp <- tmp / rowSums(tmp)

    df <- data.frame(nich.cluster = matrix(0, nrow = ncol(tmp) * nrow(tmp)),
                     features = matrix(0, nrow = ncol(tmp) * nrow(tmp)),
                     rate = matrix(0, nrow = ncol(tmp) * nrow(tmp)))

    ii = 1
    for(i in 1:nrow(tmp)){
        for(j in 1:ncol(tmp)){
            df$nich.cluster[ii] <- rownames(tmp)[i]
            df$features[ii] <- colnames(tmp)[j]
            df$rate[ii] <- tmp[[i, j]]

            ii = ii + 1
        }
    }

    cluster.levels <- sort(unique(df$nich.cluster))
    df$nich.cluster <- factor(df$nich.cluster, levels = cluster.levels)

    p <- ggplot(df, aes(x = nich.cluster, y = features, fill = rate)) +
        geom_raster() +
        theme_linedraw() + theme(panel.grid.major = element_blank()) +
        scale_y_discrete(limits = rev(levels(df$features))) +
        scale_fill_gradientn(
            name = 'Prop.',
            colours = colorRampPalette(
                colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu'))
            )(n = 100)
        ) +
        # scale_fill_gradientn(colors = c('white', '#D73027'),
        #                      name = 'Prop.') +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5,
                                         vjust = 1), axis.title.x = element_blank(),
              axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
    if(!is.null(save_path)){
        saveImage(save_path,
                  p,
                  'niche_features',
                  height = 0.5+0.11*length(features),
                  width = 2.0+0.15*kmeans.n)
    }else{
        p
    }

    if('CNV_state' %in% names(st_obj@meta.data)){
        tmp <- aggregate(st_obj$CNV_state, by=list(Niche=st_obj$NicheCluster,
                                                   group=st_obj$CNV_state), length)
        p <- ggplot(tmp, aes(x=Niche, y=x)) +
            geom_col(aes(fill=group), position = 'fill') +
            scale_fill_manual(name = 'CNV_state',
                              values = c("#E41A1C", "grey"),
                              labels = c('aneuploid', 'diploid'),
            ) +
            theme_classic() +
            # xlab('Cluster') +
            xlab(NULL) +
            ylab('Proportion') +
            theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = 'black'),
                  axis.text.y = element_text(color = 'black'))
        if(!is.null(save_path)){
            saveImage(save_path,
                      p,
                      'niche_CNV',
                      height = 1.0,
                      width = 1.0+0.25*kmeans.n)
        }else{
            p
        }
    }

    return(st_obj)
}

#' Niche analysis
#'
#' @param st_obj SeuratObject
#' @param features A vector of features
#' @param save_path Save path
#' @param kmeans.n The number of clusters to cluster
#' @param st_data_path A path containing `spatial` file and `filtered_feature_bc_matrix.h5` file
#' @param slice The slice to use
#' @param species The species this sample belongs to
#' @param pythonPath The path to the Python environment to use for the analysis.
#'
#' @export
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
st_NicheAnalysis <- function(
        st_obj,
        features,
        save_path = '.',
        coexistence.method = 'correlation',
        kmeans.n = 4,
        st_data_path = NULL,
        slice = 'slice1',
        species = 'mouse',
        pythonPath = NULL
        # condaenv = 'r-reticulate'
){
    if(!dir.exists(save_path)){
        dir.create(save_path)
    }

    # Global Nich
    # print('Coexistence analysis...')
    df <- Co_occurenceScore(
        st_obj = st_obj,
        features = features,
        save_path = save_path,
        slice = slice,
        method = coexistence.method
    )

    # Local Nich
    st_obj <- NicheCluster(
        st_obj = st_obj,
        features = features,
        save_path = save_path,
        kmeans.n = kmeans.n
    )

    # Nich interaction
    if(!is.null(st_data_path)){
        write.csv(st_obj@meta.data,
                  file.path(save_path, 'metadata.csv'),
                  row.names = TRUE)
        st_Interaction(
            st_data_path = st_data_path,
            metadata_path = file.path(save_path, 'metadata.csv'),
            label_key = 'NicheCluster',
            save_path = save_path,
            species = species,
            pythonPath = pythonPath
            # condaenv = condaenv
        )
    }

    return(st_obj)
}

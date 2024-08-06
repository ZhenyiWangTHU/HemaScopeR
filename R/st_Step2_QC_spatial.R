#' Quality control on ST data
#'
#' @param st_obj The SeuratObject of 10X visium data
#' @param output.dir A character of path to store the results and figures
#' @param min.gene An integer representing the minimum number of genes detected in a spot
#' @param max.gene An integer representing the maximum number of genes detected in a spot
#' @param min.nUMI An integer representing the minimum number of nUMI detected in a spot
#' @param max.nUMI An integer representing the maximum number of nUMI detected in a spot
#' @param min.spot An integer representing the minimum number of spots expressing each gene
#' @param species A character representing the species of sample, `human` or `mouse`
#' @param bool.remove.mito A bool value indicating whether removing mitochondrial gene
#' @param SpatialColors A function that interpolates a set of given colors to create new color palettes and color ramps
#'
#' @details
#' This function performs basic quality control on the SeuratObject of 10X visium data and
#' returns the filtered SeuratObject
#'
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
QC_Spatial <- function(
    st_obj,
    output.dir = '.',
    min.gene = 200,
    min.nUMI = 500,
    max.gene = Inf,
    max.nUMI = Inf,
    min.spot = 3,
    species = 'human', # 'human','mouse'
    bool.remove.mito = TRUE,
    SpatialColors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
){
    if(!dir.exists(output.dir)){
        dir.create(output.dir)
    }

    active.assay <- st_obj@active.assay

    ##### Plot nUMI #####
    p.nUMI.spatial <- mySpatialFeaturePlot(st_obj = st_obj,
                                           features = paste0('nCount_', active.assay),
                                           legend.name = 'nUMI',
                                           legend.color = SpatialColors)
    saveImage(output.dir,
              p.nUMI.spatial,
              'nUMI_spatial',
              height = 4,
              width = 4)

    p.nUMI.violin <- VlnPlot(st_obj, features = paste0('nCount_', active.assay)) +
        xlab(NULL) +
        ggtitle('nUMI') +
        NoLegend() +
        suppressWarnings(scale_fill_manual(values = '#7FC97F')) +
        theme(axis.text.x = element_text(angle = 0,
                                         vjust = 0,
                                         hjust = 0.5))
    if(min.nUMI > 0){
        p.nUMI.violin <- p.nUMI.violin +
            geom_hline(yintercept = min.nUMI, linetype="dashed")
    }
    if(!is.infinite(max.nUMI)){
        p.nUMI.violin <- p.nUMI.violin +
            geom_hline(yintercept = max.nUMI, linetype="dashed")
    }
    saveImage(output.dir,
              p.nUMI.violin,
              'nUMI_violin',
              height = 4,
              width = 4)

    ##### Plot nGene #####
    p.nGene.spatial <- mySpatialFeaturePlot(st_obj = st_obj,
                                            features = paste0('nFeature_', active.assay),
                                            legend.name = 'nGene',
                                            legend.color = SpatialColors)
    saveImage(output.dir,
              p.nGene.spatial,
              'nGene_spatial',
              height = 4,
              width = 4)

    p.nGene.violin <- VlnPlot(st_obj, features = paste0('nFeature_', active.assay)) +
        xlab(NULL) +
        ggtitle('nGene') +
        NoLegend() +
        suppressWarnings(scale_fill_manual(values = '#BEAED4')) +
        theme(axis.text.x = element_text(angle = 0,
                                         vjust = 0,
                                         hjust = 0.5))
    if(min.gene > 0){
        p.nGene.violin <- p.nGene.violin +
            geom_hline(yintercept = min.gene, linetype="dashed")
    }
    if(!is.infinite(max.gene)){
        p.nGene.violin <- p.nGene.violin +
            geom_hline(yintercept = max.gene, linetype="dashed")
    }
    saveImage(output.dir,
              p.nGene.violin,
              'nGene_violin',
              height = 4,
              width = 4)

    ##### Show areas with low quality #####
    st_obj$Quality <- ifelse(st_obj$nCount_Spatial > max.nUMI | st_obj$nCount_Spatial < min.nUMI |
                                 st_obj$nFeature_Spatial > max.gene | st_obj$nFeature_Spatial < min.gene,
                             'Low',
                             'High')
    suppressMessages(suppressWarnings(
        p.quality.spatial <- SpatialDimPlot(st_obj, group.by = 'Quality',
                                            pt.size.factor = 1.6, stroke = NA) +
            theme(legend.position = 'right',
                  legend.key = element_blank()) +
            scale_fill_manual(
                name = 'Quality',
                values = c('Low'='#377EB8','High'='#E41A1C')
            ) +
            guides(fill = guide_legend(override.aes = list(size = 4)))
    ))
    saveImage(output.dir,
              p.quality.spatial,
              'quality_spatial',
              height = 4,
              width = 4)

    ##### Plot nSpot #####
    if(packageVersion('SeuratObject') >= '5.0.0'){
        counts <- GetAssayData(st_obj, layer = 'data',
                               assay = 'Spatial')
    }else{
        counts <- GetAssayData(st_obj, slot = 'data',
                               assay = 'Spatial')
    }
    nSpot <- data.frame(nSpot = rowSums(as.matrix(counts > 0)))
    p.nSpot.hist <- ggplot(data = nSpot, aes(x=nSpot)) +
        geom_histogram(bins = 30, fill = '#5D69B1') +
        theme_classic()
    if(min.spot > 0){
        p.nSpot.hist <- p.nSpot.hist + geom_vline(xintercept = min.spot, linetype="dashed")
    }
    saveImage(output.dir,
              p.nSpot.hist,
              'nSpot_hist',
              height = 4,
              width = 4)

    ##### Plot mito genes #####
    if(species %in% c('human', 'mouse')){
        if(species == 'human'){
            st_obj[["percent.mt"]] <- PercentageFeatureSet(st_obj, pattern = "^MT-")
        }else if(species == 'mouse'){
            st_obj[["percent.mt"]] <- PercentageFeatureSet(st_obj, pattern = "^mt-")
        }else{
            st_obj[["percent.mt"]] <- 0
        }
        p.mito.spatial <- mySpatialFeaturePlot(st_obj = st_obj,
                                               features = 'percent.mt',
                                               legend.name = 'percent.mt',
                                               legend.color = SpatialColors)
        saveImage(output.dir,
                  p.mito.spatial,
                  'mito_spatial',
                  height = 4,
                  width = 4)
    }

    ##### Save the raw data #####
    saveRDS(st_obj, file.path(output.dir, 'raw_object.rds'))

    ##### Remove spots with low quality #####
    # st_obj <- st_obj %>% subset(nCount_Spatial < max.nUMI & nCount_Spatial > min.nUMI) %>%
    #     subset(nFeature_Spatial < max.gene & nFeature_Spatial > min.gene)
    st_obj <- st_obj %>% subset(Quality == 'High')

    ##### Remove genes with fewer occurrences #####
    filtered_genes <- rownames(nSpot)[nSpot$nSpot > min.spot]
    st_obj <- st_obj[filtered_genes, ]

    ##### Remove mito genes if bool.remove.mito #####
    if(bool.remove.mito){
        print('Removing mitochondrial genes...')
        if(species == 'human'){
            mito.genes <- grep('^MT-', rownames(st_obj), value = T)
        }else if(species == 'mouse'){
            mito.genes <- grep('^mt-', rownames(st_obj), value = T)
        }else{
            mito.genes = c()
        }
        st_obj <- st_obj[setdiff(rownames(st_obj), mito.genes), ]
    }

    ##### Save the filtered data #####
    saveRDS(st_obj, file.path(output.dir, 'object_filtered_QC.rds'))

    return(st_obj)
}

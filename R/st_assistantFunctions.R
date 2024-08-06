#### Functions ####
getDefaultClusterColor <- function(n.cluster){
    brewer.set3 <- c('#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462',
                     '#B3DE69','#FCCDE5','#BC80BD','#CCEBC5','#FFED6F')
    if(n.cluster <= 11){
        return(brewer.set3[1:n.cluster])
    }else{
        return(colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n.cluster))
    }
}

plotSingle <- function(st_obj,
                       feature,
                       legend.name = NULL,
                       legend.color = colorRampPalette(
                           colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu'))
                       ),
                       ...){
    suppressMessages(suppressWarnings(
        p <- SpatialFeaturePlot(st_obj,
                                features = feature,
                                stroke = NA,
                                ...) +
            theme(legend.position = 'right') +
            scale_fill_gradientn(
                name = ifelse(is.null(legend.name),
                              feature,
                              legend.name),
                colours = legend.color(n = 100)
            )
    ))
    return(p)
}

#' Plot features on the image
#'
#' @import Seurat
#' @import patchwork
#' @import RColorBrewer
#'
#' @export
#'
mySpatialFeaturePlot <- function(st_obj,
                                 features,
                                 legend.name = NULL,
                                 legend.color = colorRampPalette(
                                     colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu'))
                                 ),
                                 n.col = 4,
                                 ...){
    if(length(features) == 1){
        return(plotSingle(st_obj,
                          features,
                          legend.name,
                          legend.color,
                          ...))
    }else{
        plots <- vector(mode = 'list', length = length(features))
        for(i in 1:length(features)){
            plots[[i]] <- plotSingle(st_obj,
                                     features[i],
                                     ifelse(is.null(legend.name),
                                            features[i],
                                            legend.name[i]),
                                     legend.color,
                                     ...)
        }
        return(wrap_plots(plots = plots,
                          ncol = min(n.col, length(plots))))
    }
}

plotSingle_merfish <- function(st_obj,
                               feature,
                               legend.name = NULL,
                               legend.color = colorRampPalette(
                                   colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu'))
                               ),
                               ...){
    suppressMessages(suppressWarnings(
        p <- ImageFeaturePlot(st_obj,
                              features = feature,
                              ...) +
            scale_fill_gradientn(
                name = ifelse(is.null(legend.name),
                              feature,
                              legend.name),
                colours = legend.color(n = 100)
            )
    ))
    return(p)
}

#' Plot features on the image
#'
#' @import Seurat
#' @import patchwork
#' @import RColorBrewer
#'
#' @export
#'
myImageFeaturePlot <- function(st_obj,
                               features,
                               legend.name = NULL,
                               legend.color = colorRampPalette(
                                   colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu'))
                               ),
                               n.col = 4,
                               ...){
    if(length(features) == 1){
        return(plotSingle_merfish(st_obj,
                                  features,
                                  legend.name,
                                  legend.color,
                                  ...))
    }else{
        plots <- vector(mode = 'list', length = length(features))
        for(i in 1:length(features)){
            plots[[i]] <- plotSingle_merfish(st_obj,
                                             features[i],
                                             ifelse(is.null(legend.name),
                                                    features[i],
                                                    legend.name[i]),
                                             legend.color,
                                             ...)
        }
        return(wrap_plots(plots = plots,
                          ncol = min(n.col, length(plots))))
    }
}

saveImage <- function(output.dir,
                      p,
                      saveName,
                      height = 4,
                      width = 4){
    if(!dir.exists(file.path(output.dir, 'pdf'))){
        dir.create(file.path(output.dir, 'pdf'))
    }
    if(!dir.exists(file.path(output.dir, 'png'))){
        dir.create(file.path(output.dir, 'png'))
    }

    ggsave(file.path(output.dir, 'pdf', paste0(saveName, '.pdf')),
           p,
           height = height,
           width = width)
    ggsave(file.path(output.dir, 'png', paste0(saveName, '.png')),
           p,
           height = height,
           width = width)
}


#' @param envname The name of the conda environment
#' @param python_version The version of python, 3.9 by default
#' @import reticulate
#' @export
init_condaenv <- function(
        envname = 'r-reticulate',
        python_version = 3.9
){
    reticulate::conda_create(envname = envname,
                             python_version = python_version)
    reticulate::use_condaenv(condaenv = envname)
    reticulate::py_install('commot[tradeSeq]',
                           pip = TRUE)
    reticulate::py_install('cell2location[tutorials]',
                           pip = TRUE)
    # print(paste0('The conda environment named ', envname,
    #              ' is now already available for running cell2location and COMMOT'))
}

#' Convert mouse genes to human genes
#'
#' @param mouse.genes A vector of mouse genes
#' @param host Host to connect to in `useMart` function
#'
#' @import biomaRt
#'
#' @export
#'
convertMouseGene <- function(
        mouse.genes,
        host = "https://dec2021.archive.ensembl.org/"
){
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                     host = host)
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                     host = host)
    m2h.g <- getLDS(attributes = "mgi_symbol",
                    filters = "mgi_symbol",
                    mart = mouse,
                    values = mouse.genes,
                    attributesL = "hgnc_symbol",
                    martL = human,
                    uniqueRows = TRUE)

    m2h.g$Upper.symbol <- toupper(m2h.g$MGI.symbol)
    simple.genes <- m2h.g[m2h.g$HGNC.symbol == m2h.g$Upper.symbol, ]

    trans.genes <- m2h.g[!(m2h.g$MGI.symbol %in% simple.genes$MGI.symbol), ]
    trans.one.genes <- trans.genes[!duplicated(trans.genes$MGI.symbol), ]

    res.genes <- rbind(simple.genes, trans.one.genes)
    res.genes <- subset(res.genes, select = -Upper.symbol)

    return(res.genes)
}

#' Convert human genes to mouse genes
#'
#' @param human.genes A vector of human genes
#' @param host Host to connect to in `useMart` function
#'
#' @import biomaRt
#' @import stringr
#'
#' @export
#'
convertHumanGene <- function(
        human.genes,
        host = "https://dec2021.archive.ensembl.org/"
){
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                     host = host)
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                     host = host)
    h2m.g <- getLDS(attributes = "hgnc_symbol",
                    filters = "hgnc_symbol",
                    mart = human,
                    values = human.genes,
                    attributesL = "mgi_symbol",
                    martL = mouse,
                    uniqueRows = TRUE)

    h2m.g$Lower.symbol <- str_to_title(h2m.g$HGNC.symbol)
    simple.genes <- h2m.g[h2m.g$MGI.symbol == h2m.g$Lower.symbol, ]

    trans.genes <- h2m.g[!(h2m.g$HGNC.symbol %in% simple.genes$HGNC.symbol), ]
    trans.one.genes <- trans.genes[!duplicated(trans.genes$HGNC.symbol), ]

    res.genes <- rbind(simple.genes, trans.one.genes)
    res.genes <- subset(res.genes, select = -Lower.symbol)

    return(res.genes)
}

#' convertMatrixMouseGene
#'
#' Convert mouse genes to human genes in a matrix with mouse genes as rownames
#'
#' @param counts A matrix whose rownames are mouse genes
#' @param host Host to connect to in `useMart` function
convertMatrixMouseGene <- function(
        counts,
        host = "https://dec2021.archive.ensembl.org/"
){
    mouse.genes <- rownames(counts)
    res.genes <- convertMouseGene(mouse.genes = mouse.genes,
                                  host = host)

    rownames(res.genes) <- res.genes$MGI.symbol
    res.genes <- res.genes[intersect(res.genes$MGI.symbol, rownames(counts)), ]

    counts <- counts[res.genes$MGI.symbol, ]
    counts <- aggregate(counts, by = list(res.genes$HGNC.symbol), FUN = mean)
    rownames(counts) <- counts$Group.1
    counts <- subset(counts, select = -Group.1)

    counts <- as(as.matrix(counts), Class = 'dgCMatrix')

    return(counts)
}


#' convertMatrixHumanGene
#'
#' Convert human genes to mouse genes in a matrix with human genes as rownames
#'
#' @param counts A matrix whose rownames are human genes
#' @param host Host to connect to in `useMart` function
convertMatrixHumanGene <- function(
        counts,
        host = "https://dec2021.archive.ensembl.org/"
){
    human.genes <- rownames(counts)
    res.genes <- convertHumanGene(human.genes = human.genes,
                                  host = host)

    rownames(res.genes) <- res.genes$HGNC.symbol
    res.genes <- res.genes[intersect(res.genes$HGNC.symbol, rownames(counts)), ]

    counts <- counts[res.genes$HGNC.symbol, ]
    counts <- aggregate(counts, by = list(res.genes$MGI.symbol), FUN = mean)
    rownames(counts) <- counts$Group.1
    counts <- subset(counts, select = -Group.1)

    counts <- as(as.matrix(counts), Class = 'dgCMatrix')

    return(counts)
}


#' Save object to h5 file and spatial file
#'
#' @param object
#' @param file.dir
#' @param assay
#' @param slice
#'
#' @import rjson
#' @import png
#'
Rds2H5 <- function(
        object,
        file.dir = '.',
        assay = 'Spatial',
        slice = 'slice1'
){
    # H5 file
    outfile <- hdf5r::H5File$new(file.path(file.dir,
                                           'filtered_feature_bc_matrix.h5'),
                                 mode = 'w')
    matrix.grp <- outfile$create_group('matrix')
    matrix.grp[['data']] <- as.integer(object@assays[[assay]]@counts@x)
    matrix.grp[['indices']] <- object@assays[[assay]]@counts@i
    matrix.grp[['indptr']] <- object@assays[[assay]]@counts@p
    matrix.grp[['shape']] <- object@assays[[assay]]@counts@Dim
    matrix.grp[['barcodes']] <- colnames(object)
    features.grp <- outfile$create_group('matrix/features')
    features.grp[['name']] <- rownames(object)
    features.grp[['feature_type']] <- rep("Gene Expression", nrow(object))
    features.grp[['id']] <- rep("Unknown", nrow(object))
    features.grp[['genome']] <- rep("Unknown", nrow(object))
    outfile$close_all()

    # spatial file
    if(!dir.exists(file.path(file.dir, 'spatial'))){
        dir.create(file.path(file.dir, 'spatial'))
    }

    # scalefactors_json.json
    scale.factors <- object@images[[slice]]@scale.factors
    scalefactors <- vector()
    scalefactors[['tissue_hires_scalef']] <- scale.factors$hires
    scalefactors[['tissue_lowres_scalef']] <- scale.factors$lowres
    scalefactors[['fiducial_diameter_fullres']] <- scale.factors$fiducial
    scalefactors[['spot_diameter_fullres']] <- scale.factors$spot
    jsonData <- rjson::toJSON(scalefactors)
    write(jsonData, file.path(file.dir, 'spatial', 'scalefactors_json.json'))

    # tissue_positions_list.csv
    coordinates <- object@images[[slice]]@coordinates
    coordinates$barcode <- rownames(coordinates)
    coordinates <- coordinates[c("barcode","tissue","row","col",
                                 "imagerow","imagecol")]
    write.table(coordinates,
                file.path(file.dir, 'spatial', 'tissue_positions_list.csv'),
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE,
                sep = ',')

    # tissue_lowres_image.png
    png::writePNG(image = object@images[[slice]]@image,
                  target = file.path(file.dir, 'spatial',
                                     'tissue_lowres_image.png'))
    png::writePNG(image = object@images[[slice]]@image,
                  target = file.path(file.dir, 'spatial',
                                     'tissue_hires_image.png'))
}

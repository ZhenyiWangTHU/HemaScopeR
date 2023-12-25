#' run_cell_annotation
#'
#' This function annotates cells using [abcCellmap](http://scrna.sklehabc.com/) and
#' returns the Seurat object with annotation labels. To save memory, the function
#' can perform batch predictions on the data
#'
#' @param object A Seurat object containing count matrix
#' @param split A bool value indicating whether to predict cells in batches
#' @param frag_num A integer indicating the number of cells in one batch
#' @param species A character indicating the species of the sample
#' @param host A character indicating the host to connect to in `useMart` function of `biomaRt`
#' @param assay A character indicating the source `assay` of count matrix
#'
#' @value The following documents are from 'abcCellmap' package: Users can predict the cell types of hematopoietic cells by implementing two approaches (Scmap and Seurat). 
#' Cells in our ABC are labeled by 43 different RNA clusters according to unsupervised clustering of single-cell transcriptional profiles, 
#' and also labeled by 32 immunophenotypic cell types, involving HSPC, B cell, T cell, NK cell, Neutrophil, Monocyte and Erythrocyte population. The format of result are as follows:
#' queryCell, Seurat.RNACluster, Seurat.RNACluster.score, Seurat.Immunophenotype, Seurat.Immunophenotype.score, 
#' scmap.RNACluster, scmap.RNACluster.score, scmap.Immunophenotype, scmap.Immunophenotype.score, scmap.Cell, scmap.Cell.score, pertype.
#' Of which, "queryCell" is the cell information in the query data,
#' "Seurat.RNACluster" is the RNA cluster predicted by Seurat,
#' "Seurat.RNACluster.score" is the prediction score of RNA cluster by Seurat,
#' "Seurat.Immunophenotype" is the immunophenotypic cell type predicted by Seurat,
#' "Seurat.Immunophenotype.score" is the prediction score of immunophenotypic cell type by Seurat,
#' "scmap.RNACluster" is the RNA cluster predicted by scmap,
#' "scmap.RNACluster.score" is the prediction score of RNA cluster by scmap,
#' "scmap.Immunophenotype" is the immunophenotypic cell type predicted by scmap,
#' "scmap.Immunophenotype.score" is the prediction score of immunophenotypic cell type by scmap,
#' "scmap.Cell" is the nearest single cell in our ABC reference predicted by scmap,
#' "scmap.Cell.score" is the prediction score of the nearest cell by scmap,
#' "pertype" means the percentages of top 2 immunophenotypic cell types in corresponding Seurat.RNACluster result.
#'
#' @import Seurat
#' @export
run_cell_annotation <- function(
        object,
        split = TRUE,
        frag_num = 5000,
        species = 'mmu',
        host = "https://dec2021.archive.ensembl.org/",
        assay = 'RNA',
        output.dir = NULL
){
    annotation.labels <- data.frame()
    if(!split){
        res <- cell_annotation_abc(object,
                                   species = species,
                                   host = host,
                                   assay = assay)
        annotation.labels <- rbind(annotation.labels, res)
    }else{
        breaks <- seq(from = 0, to = ncol(object), by = frag_num)
        breaks <- c(breaks, ncol(object))
        for(i in 1:(length(breaks) - 1)){
            if(breaks[i+1] - (breaks[i]+1) < 1){
                next
            }
            sub.object <- object[, (breaks[i]+1):breaks[i+1]]
            res <- cell_annotation_abc(sub.object,
                                       species = species,
                                       host = host,
                                       assay = assay)
            annotation.labels <- rbind(annotation.labels, res)
        }
    }

    annotation.labels <- annotation.labels[colnames(object), ]
    object <- AddMetaData(object,
                          metadata = annotation.labels)

    pdf(paste0(output.dir,'/predicted.Seurat.RNACluster.pdf'), width = 9, height = 6)
    print(DimPlot(object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='Seurat.RNACluster', raster = FALSE))
    dev.off()
    
    png(paste0(output.dir,'/predicted.Seurat.RNACluster.png'), width = 900, height = 600)
    print(DimPlot(object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='Seurat.RNACluster', raster = FALSE))
    dev.off()

    pdf(paste0(output.dir,'/predicted.Seurat.Immunophenotype.pdf'), width = 9, height = 6)
    print(DimPlot(object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='Seurat.Immunophenotype', raster = FALSE))
    dev.off()
    
    png(paste0(output.dir,'/predicted.Seurat.Immunophenotype.png'), width = 900, height = 600)
    print(DimPlot(object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='Seurat.Immunophenotype', raster = FALSE))
    dev.off()

    pdf(paste0(output.dir,'/predicted.scmap.RNACluster.pdf'), width = 11, height = 6)
    print(DimPlot(object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='scmap.RNACluster', raster = FALSE))
    dev.off()
    
    png(paste0(output.dir,'/predicted.scmap.RNACluster.png'), width = 1100, height = 600)
    print(DimPlot(object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='scmap.RNACluster', raster = FALSE))
    dev.off()

    pdf(paste0(output.dir,'/predicted.scmap.Immunophenotype.pdf'), width = 9, height = 6)
    print(DimPlot(object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='scmap.Immunophenotype', raster = FALSE))
    dev.off()
    
    png(paste0(output.dir,'/predicted.scmap.Immunophenotype.png'), width = 900, height = 600)
    print(DimPlot(object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='scmap.Immunophenotype', raster = FALSE))
    dev.off()
    
    return(object)
}

#' cell_annotation_abc
#'
#' This function annotates cells using `abcCellmap`
#'
#' @param object A Seurat object containing count matrix
#' @param species A character indicating the species of the sample
#' @param host A character indicating the host to connect to in `useMart` function of
#' `biomaRt`
#' @param assay A character indicating the source `assay` of count matrix
#'
#' @import Seurat
#' @import abcCellmap
cell_annotation_abc <- function(
    object,
    species = 'mmu',
    host = "https://dec2021.archive.ensembl.org/",
    assay = 'RNA'
){
    require(abcCellmap)
    
    if(packageVersion('SeuratObject') >= '5.0.0'){
        counts <- GetAssayData(object, 
                               layer = 'counts',
                               assay = assay)
    }else{
        counts <- GetAssayData(object, 
                               slot = 'counts',
                               assay = assay)
    }
    
    if(species == 'mmu'){
        counts <- convertMatrixMouseGene(counts = counts,
                                         host = host)
    }else if(species == 'hsa'){
        counts <- counts
    }else{
        print('The species parameter should be "mmu" or "hsa".')
    }
        
    suppressMessages(suppressWarnings(
        res <- abcCellmap(queryData = counts)
    ))

    if(length(which(duplicated(res$queryCell)))>0){
        res <- res[-which(duplicated(res$queryCell)),]
    }

    rownames(res) <- res$queryCell
    res$label <- res$scmap.Immunophenotype
    res <- subset(res, select = -queryCell)

    return(res)
}
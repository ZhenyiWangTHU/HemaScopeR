#' Create Violin Plots for Marker Genes
#' 
#' This function generates violin plots for marker genes across different cell types using the provided expression matrix and cell type annotations.
#' 
#' @param dataMatrix A data frame or matrix representing the expression data, where rows are cells and columns are genes.
#' @param features A character vector specifying the marker genes to plot in the violin plots.
#' @param CellTypes A factor vector containing cell type annotations for each cell.
#' @param cellTypeOrders A character vector specifying the order of cell types for plotting. Defaults to unique values in `CellTypes`.
#' @param cellTypeColors A character vector specifying the colors to use for cell type groups. Defaults to a color palette.
#' @param output.dir The path to the directory where the resulting violin plots will be saved.
#' 
#' @details
#' This function creates violin plots for marker genes across different cell types. You can provide the expression data in the form of a data frame or matrix. The `features` parameter should contain the names of marker genes to be plotted. The `CellTypes` parameter should contain cell type annotations, and `cellTypeOrders` can be used to specify the order in which cell types are displayed in the plots. You can also customize the colors of cell type groups using `cellTypeColors`. The resulting plots are saved to the specified output directory.
#' 
#' @return None. This function generates and saves violin plots as PDF files.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#' 
#' @export

combinedViolinPlot = function(dataMatrix = NULL,
                              features = NULL,
                              CellTypes = NULL,
                              cellTypeOrders = NULL,
                              cellTypeColors = NULL,
                              Org = NULL,
                              output.dir = NULL){
  # data("genecode_geneSymbolandEnsembleID")
  if(Org == 'mmu'){
     load("../data/mouseGeneSymbolandEnsembleID.rdata")
     genecode_geneSymbolandEnsembleID <- mouseGeneSymbolandEnsembleID
  }else if(Org == 'hsa'){
     load("../data/humanGeneSymbolandEnsembleID.rdata")
     genecode_geneSymbolandEnsembleID <- humanGeneSymbolandEnsembleID
  }
        
  lineages_signatures <- features
  if(all(grepl('^ENSMUSG|^ENSG', rownames(dataMatrix)))){
    # if the features are ensemble id, we need convert the rownames from ensembleID to gene name, e.g. 'ENSG00000102145' -> 'Gata1'.
    lineages_signatures <- unique(lineages_signatures)
    lineages_signatures_ensembleId <- mapvalues(x = lineages_signatures,
                                                from = genecode_geneSymbolandEnsembleID$geneName,
                                                to= genecode_geneSymbolandEnsembleID$ensemblIDNoDot,
                                                warn_missing = FALSE)

    dataMatrix <- subset(dataMatrix, rownames(dataMatrix) %in% lineages_signatures_ensembleId)
    rownames(dataMatrix) <- mapvalues(x = rownames(dataMatrix),
                                      from = genecode_geneSymbolandEnsembleID$ensemblIDNoDot,
                                      to= genecode_geneSymbolandEnsembleID$geneName,
                                      warn_missing = FALSE)
  }else{
    lineages_signatures <- unique(lineages_signatures)
    dataMatrix <- subset(dataMatrix, rownames(dataMatrix) %in% lineages_signatures)
  }

  dataMatrix <- as.data.frame(t(dataMatrix))
  dataMatrix$CellTypes <- CellTypes
  dataMatrix<- reshape2::melt(dataMatrix,
                              id.vars="CellTypes",
                              variable.name="Genes",
                              value.name="Scaled_expression_value")
  dataMatrix$Genes <- factor(dataMatrix$Genes, levels = unique(lineages_signatures))
  dataMatrix <- dataMatrix[order(dataMatrix$Genes, decreasing = FALSE),]

  if(is.null(cellTypeOrders)){cellTypeOrders <- unique(CellTypes)}
  dataMatrix$CellTypes <- factor(dataMatrix$CellTypes, levels = cellTypeOrders)

  if(is.null(cellTypeColors)){cellTypeColors <- rbPal5(length(unique(CellTypes)))[as.numeric(factor(unique(CellTypes)))]}
  pdf(paste0(output.dir,'/combinedViolin.pdf'),width=(0.8*length(unique(CellTypes))),height=(0.5*length(unique(features))))
      print(ggplot(data=dataMatrix,
                   aes(x=CellTypes,
                       y=Scaled_expression_value))+
            geom_violin(aes(fill=CellTypes),size=0.1,scale="width")+
                        facet_grid(Genes~.,scales="free")+
                        scale_fill_manual(values = cellTypeColors,
                                          breaks = cellTypeOrders)+
                        theme(panel.grid.major =element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              axis.line = element_line(size = 1,
                                                       colour = "black"),
                              axis.title.x =element_text(size = 30,
                                                         family = "sans",
                                                         color = "black",
                                                         face = "bold"),
                              axis.text.x = element_text(size = 30,
                                                         family = "sans",
                                                         color = "black",
                                                         face = "bold",
                                                         angle= 315,
                                                         vjust = 1,
                                                         hjust = 0),
                              axis.text.y = element_text(size = 8,
                                                         family = "sans",
                                                         color = "black",
                                                         face = "bold",
                                                         vjust = 0,
                                                         hjust = 1),
                              axis.title.y=element_text(size=30,
                                                        family = "sans",
                                                        color = "black",
                                                        face = "bold"),
                              legend.text = element_text(size=30,
                                                         family = "sans",
                                                         color = "black",
                                                         face = "bold"),
                              legend.title = element_blank(),
                              axis.ticks = element_blank(),
                              strip.text.y = element_text(size = 10)))
  dev.off()
   
  png(paste0(output.dir,'/combinedViolin.png'),width=(80*length(unique(CellTypes))),height=(50*length(unique(features))))
      print(ggplot(data=dataMatrix,
                   aes(x=CellTypes,
                       y=Scaled_expression_value))+
            geom_violin(aes(fill=CellTypes),size=0.1,scale="width")+
                        facet_grid(Genes~.,scales="free")+
                        scale_fill_manual(values = cellTypeColors,
                                          breaks = cellTypeOrders)+
                        theme(panel.grid.major =element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              axis.line = element_line(size = 1,
                                                       colour = "black"),
                              axis.title.x =element_text(size = 30,
                                                         family = "sans",
                                                         color = "black",
                                                         face = "bold"),
                              axis.text.x = element_text(size = 30,
                                                         family = "sans",
                                                         color = "black",
                                                         face = "bold",
                                                         angle= 315,
                                                         vjust = 1,
                                                         hjust = 0),
                              axis.text.y = element_text(size = 8,
                                                         family = "sans",
                                                         color = "black",
                                                         face = "bold",
                                                         vjust = 0,
                                                         hjust = 1),
                              axis.title.y=element_text(size=30,
                                                        family = "sans",
                                                        color = "black",
                                                        face = "bold"),
                              legend.text = element_text(size=30,
                                                         family = "sans",
                                                         color = "black",
                                                         face = "bold"),
                              legend.title = element_blank(),
                              axis.ticks = element_blank(),
                              strip.text.y = element_text(size = 15, angle = 0)))
  dev.off()  
}

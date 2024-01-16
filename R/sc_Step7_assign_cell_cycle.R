#' Assign Cell Cycle
#' 
#' This function is used to assign cell cycle phases based on input single-cell RNA sequencing data.
#' 
#' @param sc_object A Seurat object containing single-cell RNA sequencing data.
#' @param counts_matrix The 'counts' slot in the Seurat object.
#' @param data_matrix The 'data' slot in the Seurat object.
#' @param cellcycleCutoff The cutoff value for distinguishing between cycling and quiescent cells. Cells with a G1G2Score below this cutoff are considered quiescent.
#' @param cellTypeOrders The order of cell types for visualization. If not provided, the function will use the unique cell types in the input Seurat object.
#' @param output.dir The path to the directory where the resulting plots will be stored.
#' @param Org A character vector specifying the species of cell cycle genes, can be 'mmu' (mouse) or 'hsa' (human).
#' 
#' @details
#' This function assigns cell cycle phases based on the input single-cell RNA sequencing data by analyzing cell cycle-related genes and generates plots of the cell cycle analysis results. 
#' It uses external data, including lists of cell cycle genes and cell cycle marker genes, which should be prepared and loaded before using the function. 
#' Note that scran was refered to https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html. Cell cycle genes are from Dysfunctional CD8 T Cells Form a Proliferative, 
#' Dynamically Regulated Compartment within Human Melanoma(Table S6)first sheet, supplementary table 1, https://www.nature.com/articles/nature20123#supplementary-information
#' 
#' @return Returns an updated Seurat object that includes the assigned cell cycle phases.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#' 
#' @export

cellCycle = function(sc_object = NULL,
                     counts_matrix = NULL,
                     data_matrix = NULL,
                     cellcycleCutoff = NULL,
                     cellTypeOrders = NULL,
                     output.dir = NULL,
                     Org=NULL){

  counts_matrix <- counts_matrix[, rownames(sc_object@meta.data)]
  rownames(counts_matrix) <- toupper(rownames(counts_matrix))
  #data("genecode_geneSymbolandEnsembleID")
  if(Org == 'mmu'){
     load("../data/mouseGeneSymbolandEnsembleID.rdata")
     rownames(counts_matrix) <- mapvalues(rownames(counts_matrix),
                                          from = mouseGeneSymbolandEnsembleID$geneName,
                                          to = mouseGeneSymbolandEnsembleID$ensemblIDNoDot,
                                          warn_missing = FALSE)
  }else if(Org == 'hsa'){
     load("../data/humanGeneSymbolandEnsembleID.rdata")
     rownames(counts_matrix) <- mapvalues(rownames(counts_matrix),
                                          from = humanGeneSymbolandEnsembleID$geneName,
                                          to = humanGeneSymbolandEnsembleID$ensemblIDNoDot,
                                          warn_missing = FALSE) 
  }
    
  scran_object <- SingleCellExperiment(list(counts=as.matrix(counts_matrix)))
    if(Org == 'mmu'){
     Org.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
        }else if(Org == 'hsa'){
     Org.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    }else{
    stop("Org should be 'mmu' or 'hsa'.")
  }
  assigned <- cyclone(scran_object, pairs=Org.pairs)
  saveRDS(assigned, file=paste0(gsub("/Step7.Assign_cell_cycles/$", "/RDSfiles", output.dir), '/assigned.rds'))  

  data_matrix <- data_matrix[, rownames(sc_object@meta.data)]

  #data("CellCycleGeneSets")  
  load("../data/CellCycleGeneSets.rdata")
  if(Org == 'hsa'){
      cycle_genes$Gene <- toupper(cycle_genes$Gene)
  }  
  cycle_genes_G1 <- subset(cycle_genes, cycle_genes$Phase=='G1-S')
  cycle_genes_G2M <- subset(cycle_genes, cycle_genes$Phase=='G2-M')

  G1_cycle_genes <- subset(data_matrix, rownames(data_matrix)%in%cycle_genes_G1$Gene)
  G2M_cycle_genes <- subset(data_matrix, rownames(data_matrix)%in%cycle_genes_G2M$Gene)
  print('G1_cycle_genes')
  print(dim(G1_cycle_genes))
  print('G2M_cycle_genes')
  print(dim(G2M_cycle_genes))   
  G1_cycle_genes_colMeans <- colMeans(G1_cycle_genes)
  G1_cycle_genes_colMeans <- as.data.frame(G1_cycle_genes_colMeans)
  G2M_cycle_genes_colMeans <- colMeans(G2M_cycle_genes)
  G2M_cycle_genes_colMeans <- as.data.frame(G2M_cycle_genes_colMeans)
  G1AndG2M_cycle_genes_colMeans <- G1_cycle_genes_colMeans
  G1AndG2M_cycle_genes_colMeans$G2M_cycle_genes_colMeans <- G2M_cycle_genes_colMeans$G2M_cycle_genes_colMeans
  G1AndG2M_cycle_genes_colMeans$G1G2Score <- rowMeans(G1AndG2M_cycle_genes_colMeans)
  G1AndG2M_cycle_genes_colMeans$cellTypes <- mapvalues(rownames(G1AndG2M_cycle_genes_colMeans),
                                                       from=rownames(sc_object@meta.data),
                                                       to=as.character(sc_object@meta.data$selectLabels),
                                                       warn_missing=FALSE)

  if(is.null(cellcycleCutoff)){
        pdf(paste0(output.dir, '/G1G2Score.pdf'), width = 5, height = 5)
         print(ggplot(G1AndG2M_cycle_genes_colMeans, aes(G1G2Score)) +
               geom_histogram(aes(y=after_stat(density)),
                       binwidth=0.05,
                       colour="black",
                       fill="white") +
               geom_density(alpha=0.2, fill="#666666", colour=c('#08519c'), size=1) + mytheme)
         dev.off()

        png(paste0(output.dir, '/G1G2Score.png'), width = 500, height = 500)
         print(ggplot(G1AndG2M_cycle_genes_colMeans, aes(G1G2Score)) +
               geom_histogram(aes(y=after_stat(density)),
                       binwidth=0.05,
                       colour="black",
                       fill="white") +
               geom_density(alpha=0.2, fill="#666666", colour=c('#08519c'), size=1) + mytheme)
         dev.off()
  }else{
        pdf(paste0(output.dir, '/G1G2Score.pdf'), width = 5, height = 5)
         print(ggplot(G1AndG2M_cycle_genes_colMeans, aes(G1G2Score)) +
               geom_histogram(aes(y=after_stat(density)),
                       binwidth=0.05,
                       colour="black",
                       fill="white") +
               geom_density(alpha=0.2, fill="#666666", colour=c('#08519c'), size=1) +
               geom_vline(xintercept=c(cellcycleCutoff, cellcycleCutoff),
                          linetype="longdash", 
                          linewidth=1) + mytheme)
         dev.off()

        png(paste0(output.dir, '/G1G2Score.png'), width = 500, height = 500)
         print(ggplot(G1AndG2M_cycle_genes_colMeans, aes(G1G2Score)) +
               geom_histogram(aes(y=after_stat(density)),
                       binwidth=0.05,
                       colour="black",
                       fill="white") +
               geom_density(alpha=0.2, fill="#666666", colour=c('#08519c'), size=1) +
               geom_vline(xintercept=c(cellcycleCutoff, cellcycleCutoff),
                          linetype="longdash", 
                          linewidth=1) + mytheme)
         dev.off()
  }

  # assign cell cycle
  sc_object@meta.data$cellCycle <- assigned$phases

  G1AndG2M_cycle_genes_colMeans$phase <- 'cycling'
  if(!is.null(cellcycleCutoff)){
        G1AndG2M_cycle_genes_colMeans$phase[which(G1AndG2M_cycle_genes_colMeans$G1G2Score<=cellcycleCutoff)] <- 'quiesent'
  }
  G1AndG2M_cycle_genes_colMeans$scran <- sc_object@meta.data$cellCycle
  G1AndG2M_cycle_genes_colMeans$scranAndScore <- G1AndG2M_cycle_genes_colMeans$scran
  G1AndG2M_cycle_genes_colMeans$scranAndScore[which((G1AndG2M_cycle_genes_colMeans$scran=='G1')&(G1AndG2M_cycle_genes_colMeans$phase=='quiesent'))] <- 'G0'
  if(is.null(cellTypeOrders)){cellTypeOrders <- unique(sc_object@meta.data$selectLabels)}
  G1AndG2M_cycle_genes_colMeans$cellTypes <- factor(sc_object@meta.data$selectLabels, levels=cellTypeOrders)
  G1AndG2M_cycle_genes_colMeans$count <- 1
  if(length(which(is.na(G1AndG2M_cycle_genes_colMeans$scran)==TRUE))>0){
    G1AndG2M_cycle_genes_colMeansNoNa <- G1AndG2M_cycle_genes_colMeans[-which(is.na(G1AndG2M_cycle_genes_colMeans$scran)==TRUE),]
    }else{
    G1AndG2M_cycle_genes_colMeansNoNa <- G1AndG2M_cycle_genes_colMeans
    }

  pdf(paste0(output.dir, '/CellCycle.pdf'), width = (0.5*length(unique(G1AndG2M_cycle_genes_colMeans$cellTypes))), height = 5)
     print(ggplot(G1AndG2M_cycle_genes_colMeansNoNa, aes(cellTypes,count,fill=scranAndScore))+
       geom_bar(stat="identity",position="fill")+
       scale_fill_manual(values = c(
       G0="#abd9e9",G1="#2c7bb6",
       G2M="#d7191c",S="#fdae61"
       ))+guides(fill=guide_legend(title=NULL))+mytheme+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust =0.5),
                                                              axis.ticks.length = unit(0.5,'cm')))
  dev.off()

  png(paste0(output.dir, '/CellCycle.png'), width = (50*length(unique(G1AndG2M_cycle_genes_colMeans$cellTypes))), height = 500)
     print(ggplot(G1AndG2M_cycle_genes_colMeansNoNa, aes(cellTypes,count,fill=scranAndScore))+
       geom_bar(stat="identity",position="fill")+
       scale_fill_manual(values = c(
       G0="#abd9e9",G1="#2c7bb6",
       G2M="#d7191c",S="#fdae61"
       ))+guides(fill=guide_legend(title=NULL))+mytheme+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust =0.5),
                                                              axis.ticks.length = unit(0.5,'cm')))
  dev.off()

  return(sc_object)
}

#' Run Monocle Analysis
#'
#' Function \code{run_monocle} runs the Monocle analysis on single-cell RNA-seq data to infer cell trajectories and explore gene expression changes along these trajectories.
#' 
#' @param cellData A matrix of gene expression values, where columns represent cells and rows represent genes.
#' @param phenoData A data frame containing cell metadata, such as cell labels or other relevant information.
#' @param featureData A data frame containing information about features (genes) in the dataset.
#' @param lowerDetectionLimit The lower detection limit for gene expression. Genes with expression values below this limit will be treated as non-detected.
#' @param expressionFamily The family of the expression distribution used in Monocle analysis.
#' @param cellTypes A character vector specifying cell types or labels used for coloring in trajectory plots.
#' @param output.dir The directory where the Monocle analysis results and plots will be saved.
#' 
#' @details
#' This function performs Monocle analysis on single-cell RNA-seq data. It involves several steps:
#' 
#' 1. Creating a Monocle cellDataSet using the provided cellData, phenoData, and featureData.
#' 2. Estimating size factors, dispersions, and detecting highly variable genes.
#' 3. Performing differential gene expression analysis to identify genes associated with cell state changes.
#' 4. Ordering cells along the inferred trajectories and reducing dimensionality.
#' 5. Generating and saving trajectory plots, including cell trajectory by "State" and by "Cell Types."
#' 
#' @return None. This function generates and saves Monocle analysis results and trajectory plots.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#' 
#' @export

run_monocle = function(cellData = NULL,
                       phenoData = NULL,
                       featureData = NULL,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = VGAM::negbinomial.size(),
                       cellTypes=NULL,
                       monocle.orders=NULL,
                       monocle.colors = NULL,
                       output.dir=getwd()){
  if(is.null(monocle.orders)){monocle.orders <- unique(phenoData[,cellTypes])}  
  phenoData[,cellTypes] <- factor(phenoData[,cellTypes], levels=monocle.orders)
  monocleobject <- newCellDataSet(cellData = cellData,
                                  phenoData =  new("AnnotatedDataFrame", data = phenoData),
                                  featureData =  new("AnnotatedDataFrame", data = featureData),
                                  lowerDetectionLimit = lowerDetectionLimit,
                                  expressionFamily = expressionFamily)

  monocleobject <- estimateSizeFactors(monocleobject)
  monocleobject <- estimateDispersions(monocleobject)
  monocleobject <- detectGenes(monocleobject, min_expr = 0.1)
  pdata <- pData(monocleobject)
  fdata <- fData(monocleobject)
  expressed_genes <- row.names(subset(fData(monocleobject),num_cells_expressed >= 10))

  diff_test_res <- differentialGeneTest(monocleobject[expressed_genes,],
                                        fullModelFormulaStr = paste0('~', cellTypes),
                                        cores = 1)
  diff_test_res_order <- diff_test_res[order(diff_test_res[,4], decreasing = F),]
  diff_test_res_order <- subset(diff_test_res_order, qval < 1e-10)
  ordering_genes <- rownames(diff_test_res_order)[1:1000]
  monocleobject <- setOrderingFilter(monocleobject, ordering_genes)
  #plot_ordering_genes(monocleobject)
  monocleobject <- reduceDimension(monocleobject, max_components = 10, method = 'DDRTree')
  monocleobject <- orderCells(monocleobject)
  saveRDS(monocleobject, file = paste0(gsub("/Step12.Construct_trajectories/monocle2/$", "/RDSfiles", output.dir), '/monocleobject.rds'))
    
  pdf(paste0(output.dir,'/monocle_cell_trajectory_state.pdf'), width = 8, height = 8)
    print(plot_cell_trajectory(monocleobject, color_by = "State", cell_size=2))
  dev.off()
   
  png(paste0(output.dir,'/monocle_cell_trajectory_state.png'), width = 800, height = 800)
    print(plot_cell_trajectory(monocleobject, color_by = "State", cell_size=2))
  dev.off()  

  if(is.null(monocle.colors)){
      pdf(paste0(output.dir,'/monocle_cell_trajectory_cellTypes.pdf'), width = 8, height = 8)
        print(plot_cell_trajectory(monocleobject, color_by = cellTypes, cell_size=2))
      dev.off()
        
      png(paste0(output.dir,'/monocle_cell_trajectory_cellTypes.png'), width = 800, height = 800)
        print(plot_cell_trajectory(monocleobject, color_by = cellTypes, cell_size=2)) 
      dev.off()
  }else{
      pdf(paste0(output.dir,'/monocle_cell_trajectory_cellTypes.pdf'), width = 8, height = 8)
        print(plot_cell_trajectory(monocleobject, color_by = cellTypes, cell_size=2) + scale_color_manual(values = monocle.colors))
      dev.off()
        
      png(paste0(output.dir,'/monocle_cell_trajectory_cellTypes.png'), width = 800, height = 800)
        print(plot_cell_trajectory(monocleobject, color_by = cellTypes, cell_size=2) + scale_color_manual(values = monocle.colors))
      dev.off()
  }
}

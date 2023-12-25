#' Run Slingshot Analysis
#' 
#' Function \code{run_slingshot} performs Slingshot analysis on single-cell RNA-seq data to infer cell trajectories and visualize lineage relationships.
#' 
#' @param slingshot.PCAembeddings A matrix containing the PCA embeddings of the single-cell data, typically obtained from dimensionality reduction techniques like PCA.
#' @param slingshot.cellTypes A character vector specifying cell types or labels for each cell.
#' @param slingshot.start.clus A character vector specifying the initial cluster(s) from which cell trajectories should start.
#' @param slingshot.end.clus A character vector specifying the target cluster(s) where cell trajectories should end.
#' @param slingshot.colors A vector of colors corresponding to cell types for plotting. If not provided, default colors will be used.
#' @param output.dir The directory where Slingshot analysis results and plots will be saved.
#' 
#' @details
#' Slingshot is used to infer cell trajectories and lineage relationships from single-cell RNA-seq data. This function performs the following steps:
#' 
#' 1. Constructs a Slingshot object using PCA embeddings, cell types, start clusters, and end clusters.
#' 2. Computes and plots the trajectory curves.
#' 3. Computes and plots pseudotime values along the trajectory.
#' 
#' For more details on Slingshot, refer to the package documentation and vignettes. https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html and https://bustools.github.io/BUS_notebooks_R/slingshot.html
#' 
#' @return None. This function generates and saves Slingshot analysis results and trajectory plots.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#' 
#' @export

# slingshot--------------------------------------------------------------------------------------------------------------------------------
run_slingshot = function(slingshot.PCAembeddings = NULL,
                         slingshot.cellTypes = NULL,
                         slingshot.start.clus = NULL,
                         slingshot.end.clus = NULL,
                         slingshot.colors = NULL,
                         output.dir = getwd()){
  # avoid system is computationally singular
  label_counts <- table(slingshot.cellTypes)
  unique_labels <- names(label_counts[label_counts <= 20])
  positions <- which(slingshot.cellTypes %in% unique_labels)
  if(length(positions)>0){
      slingshot.PCAembeddings <- slingshot.PCAembeddings[-positions,]
      slingshot.cellTypes <- slingshot.cellTypes[-positions]  
  }
    
  slingshot_object <- getLineages(data = slingshot.PCAembeddings,
                                  clusterLabels = slingshot.cellTypes,
                                  start.clus = slingshot.start.clus,
                                  end.clus = slingshot.end.clus)

  crv_slingshot_object <- getCurves(slingshot_object, approx_points = 500)

  if(is.null(slingshot.colors)){
    slingshot.colors <- viridis(length(unique(slingshot.cellTypes)),
                                option = "viridis")[as.numeric(factor(slingshot.cellTypes))]
  }
  pdf(paste0(output.dir,'/slingshot_curve.pdf'), width = 8, height = 8)
      plot(slingshot.PCAembeddings[, 1:2],
           col = slingshot.colors,
           asp = 1,
           pch = 16)
      lines(SlingshotDataSet(crv_slingshot_object), lwd = 2, col = 'black')
      legend("topright", legend = unique(slingshot.cellTypes), col = unique(slingshot.colors), pch = 16)
  dev.off()

  png(paste0(output.dir,'/slingshot_curve.png'), width = 800, height = 800)
      plot(slingshot.PCAembeddings[, 1:2],
           col = slingshot.colors,
           asp = 1,
           pch = 16)
      lines(SlingshotDataSet(crv_slingshot_object), lwd = 2, col = 'black')
      legend("topright", legend = unique(slingshot.cellTypes), col = unique(slingshot.colors), pch = 16)
  dev.off()  

  colors <- rev(viridis(100, option = "plasma"))
  pseudotime <- slingPseudotime(crv_slingshot_object)
  pseudotime[is.na(pseudotime)] <- 0
  pseudotime_rowmeans <- apply(pseudotime, 1, function(row) {
      non_zero_vals <- row[row != 0]
      if (length(non_zero_vals) > 0) {
        mean(non_zero_vals)
      } else {
        0 
      }
  })
  plotcol <- colors[cut(pseudotime_rowmeans, breaks=100)]

  pdf(paste0(output.dir,'/slingshot_pseudotime.pdf'), width = 8, height = 8)
      plot(slingshot.PCAembeddings[, 1:2], col = plotcol, pch=16, asp = 1)
  dev.off()

  png(paste0(output.dir,'/slingshot_pseudotime.png'), width = 800, height = 800)
      plot(slingshot.PCAembeddings[, 1:2], col = plotcol, pch=16, asp = 1)
  dev.off()  
}

#' This function executes the workflow of doubletFinder.
#'
#' Function \code{doublerFinderwraper} takes a Seurat object as input, automatically identifies doublets 
#' within it, and then adds a column 'doublet' to the meta.data to distinguish between singlets and doublets.
#' 
#' @param seuratObject The input Seurat object.
#' @param percentage Assuming 'percentage' doublet formation rate - tailor for your dataset. The default value is 0.05. 
#' @param doublerFinderwraper.PCs Which dimensions to use as input features for doubletFinder.
#' @param doublerFinderwraper.pN The percentage of real-artifical data for doubletFinder.
#' @param doublerFinderwraper.pK The pK parameter controls the doublet cell detection by determining the number of nearest neighbors and influencing the calculation of pANN scores and the final cell classification results. Adjusting the pK value allows optimization of the doublet cell detection process based on specific data and analysis requirements.
#'
#' @details
#' 
#' @return A seurat object with 'singlet' and 'doublet' labels in the column named 'doublet' in meta.data.
#' 
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
#' @export

## use doubletFinder to remove doublet---------------------------------------------------------------------------------------------------------
doublerFinderwraper = function(seuratObject = NULL,
                               percentage = NULL,
                               doublerFinderwraper.PCs = NULL,
                               doublerFinderwraper.pN = NULL,
                               doublerFinderwraper.pK = NULL
                               ){

  # pK Identification (no ground-truth)
  sweep.res.seuratObject <- paramSweep_v3(seuratObject,
                                          PCs = doublerFinderwraper.PCs,
                                          sct = FALSE)
  sweep.stats_seuratObject <- summarizeSweep(sweep.res.seuratObject,
                                             GT = FALSE)
  bcmvn_seuratObject <- find.pK(sweep.stats_seuratObject)

  # homotypic Doublet Proportion Estimate
  seuratObject.annotations <- seuratObject@meta.data$seurat_clusters
  seuratObject.homotypic.prop <- modelHomotypic(seuratObject.annotations)
  seuratObject.nExp_poi <- round(percentage*nrow(seuratObject@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
  seuratObject.nExp_poi.adj <- round(seuratObject.nExp_poi*(1-seuratObject.homotypic.prop))

  # run DoubletFinder with varying classification stringencies
  seuratObject <- doubletFinder_v3(seuratObject,
                                   PCs = doublerFinderwraper.PCs,
                                   pN = doublerFinderwraper.pN,
                                   pK = doublerFinderwraper.pK,
                                   nExp = seuratObject.nExp_poi,
                                   reuse.pANN = FALSE,
                                   sct = FALSE)
  # print(seuratObject.nExp_poi.adj)

  matching_col <- grepl("DF.classifications", colnames(seuratObject@meta.data))
  seuratObject@meta.data$doublet <- seuratObject@meta.data[,matching_col]

  return(seuratObject)
}

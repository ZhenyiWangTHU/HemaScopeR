#' The quality control step for single scRNA-Seq dataset.
#'
#' Function \code{QC_single_scRNASeq} perform the quality control step for single scRNA-Seq dataset,
#' which includes normalizing data, finding variable features, scaling data, doing principle component analysis,
#' calculating K nearest neighbors, clustering, and removing doublets.
#' 
#' @param sc_object The input Seurat object.
#' @param datasetID The ID of the dataset (e.g. 'ctrl' or 'case', ...). Please do not include ':' in the datasetID.
#' @param output.dir The path of directory for saving putputs.
#' @param Step2_Quality_Control.RemoveDoublets Logical, whether to remove doublets during quality control.
#' @param nFeature_RNA.limit The cutoff of the minimum number of detected genes in each cell.
#' @param percent.mt.limit The cutoff of the maximum percentage of mitochondria genes in each cell. 
#' @param scale.factor The scale factor for the 'data' slot in the seurat object.
#' @param nfeatures The number of selected highly variable features for down stream analysis.
#' @param ndims The number of principle components in PCA.
#' @param vars.to.regress Variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito. (ScaleData in Seurat)
#' @param PCs Which dimensions to use as input features.(RunTSNE and RunUMAP in Seurat)
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. (FindClusters in Seurat)
#' @param n.neighbors Defines k for the k-nearest neighbor algorithm. (FindNeighbors in Seurat)
#' @param percentage  Assuming 'percentage' doublet formation rate - tailor for your dataset. The default value is 0.05. 
#' @param doublerFinderwraper.PCs Which dimensions to use as input features for doubletFinder.
#' @param doublerFinderwraper.pN The percentage of real-artifical data for doubletFinder.
#' @param doublerFinderwraper.pK The pK parameter controls the doublet cell detection by determining the number of nearest neighbors and influencing the calculation of pANN scores and the final cell classification results. Adjusting the pK value allows optimization of the doublet cell detection process based on specific data and analysis requirements.
#' 
#' @details
#' The function \code{QC_single_scRNASeq} takes a Seurat object representing a single scRNA-Seq dataset as input and performs the following quality control steps:
#' 1. Normalize data using the LogNormalize method.
#' 2. Find variable features using the vst method.
#' 3. Scale data using the identified variable features and specified variables to regress out.
#' 4. Perform principal component analysis (PCA) on the scaled data.
#' 5. Find K nearest neighbors based on PCA dimensions.
#' 6. Perform clustering analysis based on the found neighbors.
#' 7. Optionally, remove doublets using doubletFinder.
#'
#' The resulting Seurat object contains data after quality control.
#'
#' @return Return a Seurat object containing data after quality control.
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#' @export

QC_single_scRNASeq = function(sc_object=NULL,
                              datasetID = NULL,
                              output.dir=NULL,
                              Step2_Quality_Control.RemoveDoublets = FALSE,
                              nFeature_RNA.limit=NULL,
                              percent.mt.limit=NULL,
                              scale.factor=NULL,
                              nfeatures=NULL,
                              ndims=NULL,
                              vars.to.regress=NULL,
                              PCs=NULL,
                              resolution = resolution,
                              n.neighbors=NULL,
                              percentage = NULL,
                              doublerFinderwraper.PCs = NULL,
                              doublerFinderwraper.pN = NULL,
                              doublerFinderwraper.pK = NULL){
    
  sc_object@meta.data$datasetID <- datasetID
  sc_object <- RenameCells(object = sc_object, add.cell.id = datasetID)  
    
  pdf(paste0(output.dir, '/sc_object_','nFeature_nCount_percentMt.pdf'), width = 12, height = 4)
    print(VlnPlot(sc_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0, ncol = 3))
  dev.off()

  png(paste0(output.dir, '/sc_object_','nFeature_nCount_percentMt.png'), width = 1200, height = 400)
    print(VlnPlot(sc_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0, ncol = 3))
  dev.off()  

  sc_object <- subset(sc_object, subset = ((nFeature_RNA > nFeature_RNA.limit) & (percent.mt < percent.mt.limit)))

  sc_object.plot1 <- FeatureScatter(sc_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  sc_object.plot2 <- FeatureScatter(sc_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

  pdf(paste0(output.dir, '/sc_object_','FeatureScatters.pdf'), width = 8, height = 4)
    print(CombinePlots(plots = list(sc_object.plot1, sc_object.plot2)))
  dev.off()

  png(paste0(output.dir, '/sc_object_','FeatureScatters.png'), width = 800, height = 400)
    print(CombinePlots(plots = list(sc_object.plot1, sc_object.plot2)))
  dev.off()  

  sc_object <- NormalizeData(sc_object, normalization.method = "LogNormalize", scale.factor = scale.factor)
  sc_object <- FindVariableFeatures(sc_object, selection.method = "vst", nfeatures = nfeatures)

  # identify the 10 most highly variable genes
  sc_object.top10 <- head(VariableFeatures(sc_object), 10)

  # plot variable features with and without labels
  sc_object.plot1 <- VariableFeaturePlot(sc_object)
  sc_object.plot2 <- LabelPoints(plot = sc_object.plot1, points = sc_object.top10, repel = TRUE)

  pdf(paste0(output.dir, '/sc_object_','VariableFeaturePlot.pdf'), width = 4, height = 4)
    print(sc_object.plot2)
  dev.off()

  png(paste0(output.dir, '/sc_object_','VariableFeaturePlot.png'), width = 400, height = 400)
    print(sc_object.plot2)
  dev.off()  

  sc_object <- ScaleData(sc_object, features = rownames(sc_object), vars.to.regress =  vars.to.regress)

  # PCA
  sc_object <- RunPCA(sc_object, npcs = ndims, features = VariableFeatures(object = sc_object), verbose=FALSE)

  pdf(paste0(output.dir, '/sc_object_','ElbowPlot.pdf'), width = 4, height = 4)
    print(ElbowPlot(sc_object, ndims = ndims))
  dev.off()

  png(paste0(output.dir, '/sc_object_','ElbowPlot.png'), width = 400, height = 400)
    print(ElbowPlot(sc_object, ndims = ndims))
  dev.off()  

  # FindNeighbors
  sc_object <- FindNeighbors(sc_object, dims = PCs, k.param=n.neighbors, force.recalc = TRUE)

  # Clustering
  sc_object <- FindClusters(sc_object, resolution = resolution)

  # run TSNE
  sc_object <- RunTSNE(object = sc_object, dims=PCs)

  # run UMAP
  sc_object <- RunUMAP(object = sc_object, dims=PCs, n.neighbors=n.neighbors, min.dist=1, seed.use=2023L)

  # remove doublet
  sc_object <- doublerFinderwraper(seuratObject = sc_object,
                                   percentage = percentage,
                                   doublerFinderwraper.PCs = doublerFinderwraper.PCs,
                                   doublerFinderwraper.pN = doublerFinderwraper.pN,
                                   doublerFinderwraper.pK = doublerFinderwraper.pK)
  color_vector <- c('#bdbdbd', '#a50f15')
  names(color_vector) <- c('Singlet', 'Doublet')
 
  pdf(paste0(output.dir, '/sc_object_','findDoublets.pdf'),width=6,height=5)
  print(DimPlot(sc_object, reduction = "umap", label = FALSE, pt.size = 0.5, cols = color_vector, group.by='doublet', raster = FALSE)+
    ggtitle(paste(as.character(length(which(sc_object@meta.data$doublet == 'Doublet'))), 'doublets and',
                  as.character(length(which(sc_object@meta.data$doublet == 'Singlet'))), 'singlets')))
  dev.off()

  png(paste0(output.dir, '/sc_object_','findDoublets.png'),width=600,height=500)
    print(DimPlot(sc_object, reduction = "umap", label = FALSE, pt.size = 0.5, cols = color_vector, group.by='doublet', raster = FALSE)+
      ggtitle(paste(as.character(length(which(sc_object@meta.data$doublet == 'Doublet'))), 'doublets and',
                  as.character(length(which(sc_object@meta.data$doublet == 'Singlet'))), 'singlets')))
  dev.off()  

  if(Step2_Quality_Control.RemoveDoublets){
      print(paste('remove', as.character(length(which(sc_object@meta.data$doublet == 'Doublet'))), 'doublets.'))
      sc_object <- subset(sc_object, subset = doublet=='Singlet')
  }  

  saveRDS(sc_object,file = paste0(gsub("/Step2.Quality_control/$", "/RDSfiles", output.dir), '/datasets.before.batch.removal.rds'))
  saveRDS(sc_object, file=paste0(gsub("/Step2.Quality_control/$", "/RDSfiles", output.dir), '/sc_objectAfterQC.RDS'))

  return(sc_object)
}
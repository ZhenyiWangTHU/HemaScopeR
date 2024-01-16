#' The quality control step for multiple scRNA-Seq datasets.
#'
#' Function \code{QC_multiple_scRNASeq} perform the quality control step for multiple scRNA-Seq datasets,
#' which includes normalizing data, finding variable features, scaling data, doing principle component analysis,
#' calculating K nearest neighbors, clustering, removing doublets, and integrating multiple datasets.
#' 
#' @param seuratObjects The input Seurat objects.
#' @param datasetID The ID of each dataset (e.g. c('ctrl', 'case')). Please do not include ':' in the datasetID.
#' @param output.dir The path of directory for saving putputs.
#' @param Step2_Quality_Control.RemoveBatches A logical value indicating whether to perform batch removal. Default is FALSE.
#' @param Step2_Quality_Control.RemoveDoublets A logical value indicating whether to remove doublets. Default is FALSE.
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
#' 
#' @return Return a Seurat object containing data after quality control.
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#' @export

QC_multiple_scRNASeq = function(seuratObjects = NULL,
                                datasetID = NULL,
                                output.dir = NULL,
                                Step2_Quality_Control.RemoveBatches = FALSE,
                                Step2_Quality_Control.RemoveDoublets = FALSE,
                                nFeature_RNA.limit = NULL,
                                percent.mt.limit = NULL,
                                scale.factor = NULL,
                                nfeatures = NULL,
                                ndims = NULL,
                                vars.to.regress = NULL,
                                PCs = NULL,
                                resolution = NULL,
                                n.neighbors = NULL,
                                percentage = NULL,
                                doublerFinderwraper.PCs = NULL,
                                doublerFinderwraper.pN = NULL,
                                doublerFinderwraper.pK = NULL
                                ){
  for(i in 1:length(seuratObjects)){
    seuratObjects[[i]]@meta.data$datasetID <- datasetID[i]
  }

  for(i in 1:length(seuratObjects)){
    seuratObject.temp <- seuratObjects[[i]]
    ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
    seuratObject.temp <- NormalizeData(seuratObject.temp, normalization.method = "LogNormalize", scale.factor = scale.factor)
    seuratObject.temp <- FindVariableFeatures(seuratObject.temp, selection.method = "vst", nfeatures = nfeatures)
    seuratObject.temp <- ScaleData(seuratObject.temp, features = rownames(seuratObject.temp) ,vars.to.regress =  vars.to.regress)
    seuratObject.temp <- RunPCA(seuratObject.temp, features = VariableFeatures(object = seuratObject.temp), verbose = FALSE)
    seuratObject.temp <- FindNeighbors(seuratObject.temp, dims = PCs, k.param = n.neighbors)
    seuratObject.temp <- FindClusters(seuratObject.temp, resolution = resolution)
    # remove doublet
    seuratObject.temp <- doublerFinderwraper(seuratObject = seuratObject.temp,
                                             percentage = percentage,
                                             doublerFinderwraper.PCs = doublerFinderwraper.PCs,
                                             doublerFinderwraper.pN = doublerFinderwraper.pN,
                                             doublerFinderwraper.pK = doublerFinderwraper.pK)
    seuratObjects[[i]] <- seuratObject.temp
  }

  seuratObjectList.temp <- c()
  for(i in 2:length(seuratObjects)){
      seuratObjectList.temp <- c(seuratObjectList.temp, seuratObjects[[i]])
  }  
  # merge multiple Seurat objects
  merged.data <- merge(x = seuratObjects[[1]],
                       y = seuratObjectList.temp,
                       add.cell.ids = datasetID,
                       project = 'merged.data')
    
  # get the longest datasetID
  datasetID_lengths <- sapply(unique(merged.data@meta.data$datasetID), nchar)
  max_length <- max(datasetID_lengths)  
  pdf(paste0(output.dir, '/datasets.before.batch.removal_', 'nFeature_nCount_percentMt.pdf'), width = (12+max_length*0.1), height = 5)
    print(VlnPlot(merged.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0, ncol = 3, group.by = "datasetID"))
  dev.off()

  png(paste0(output.dir, '/datasets.before.batch.removal_', 'nFeature_nCount_percentMt.png'), width = (1200+max_length*10), height = 500)
    print(VlnPlot(merged.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0, ncol = 3, group.by = "datasetID"))
  dev.off()  

  merged.data <- subset(merged.data, subset = nFeature_RNA > nFeature_RNA.limit & percent.mt < percent.mt.limit)

  # before batch removal
  merged.data <- NormalizeData(merged.data, normalization.method = "LogNormalize", scale.factor = scale.factor)
  merged.data <- FindVariableFeatures(merged.data, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)

  # Identify the 10 most highly variable genes
  merged.datatop10 <- head(VariableFeatures(merged.data), 10)

  # plot variable features with and without labels
  merged.dataplot1 <- VariableFeaturePlot(merged.data)
  merged.dataplot2 <- LabelPoints(plot = merged.dataplot1, points = merged.datatop10, repel = TRUE)
  pdf(paste0(output.dir, '/datasets.before.batch.removal_variablefeatures.pdf'), width = 8, height = 5)
    print(merged.dataplot2)
  dev.off()

  png(paste0(output.dir, '/datasets.before.batch.removal_variablefeatures.png'), width = 800, height = 500)
    print(merged.dataplot2)
  dev.off()
    
  merged.data <- ScaleData(merged.data, features = rownames(merged.data) ,vars.to.regress = vars.to.regress)
  merged.data <- RunPCA(merged.data, npcs = ndims, features = VariableFeatures(object = merged.data), verbose=FALSE)

  pdf(paste0(output.dir, '/datasets.before.batch.removal_ElbowPlot.pdf'), width = 5, height = 5)
    print(ElbowPlot(merged.data, ndims = ndims))
  dev.off()

  png(paste0(output.dir, '/datasets.before.batch.removal_ElbowPlot.png'), width = 500, height = 500)
    print(ElbowPlot(merged.data, ndims = ndims))
  dev.off()  

  # visualization
  merged.data <- RunTSNE(object = merged.data, dims=PCs, check_duplicates = FALSE)
  merged.data <- RunUMAP(object = merged.data, dims=PCs,
                         n.neighbors = n.neighbors, min.dist=1, seed.use=2023L)
    
  pdf(paste0(output.dir, '/datasets.before.batch.removal_tsne_datasetID.pdf'), width = (5+max_length*0.1), height = 5)
    print(DimPlot(merged.data, reduction = "tsne", label = FALSE, pt.size = 0.1, group.by='datasetID', raster = FALSE))
  dev.off()

  png(paste0(output.dir, '/datasets.before.batch.removal_tsne_datasetID.png'), width = (500+max_length*10), height = 500)
    print(DimPlot(merged.data, reduction = "tsne", label = FALSE, pt.size = 0.1, group.by='datasetID', raster = FALSE))
  dev.off()  

  pdf(paste0(output.dir, '/datasets.before.batch.removal_umap_datasetID.pdf'), width = (5+max_length*0.1), height = 5)
    print(DimPlot(merged.data, reduction = "umap", label = FALSE, pt.size = 0.1, group.by='datasetID', raster = FALSE))
  dev.off()

  png(paste0(output.dir, '/datasets.before.batch.removal_umap_datasetID.png'), width = (500+max_length*10), height = 500)
    print(DimPlot(merged.data, reduction = "umap", label = FALSE, pt.size = 0.1, group.by='datasetID', raster = FALSE))
  dev.off()  

  saveRDS(merged.data,
          file = paste0(gsub("/Step2.Quality_control/$", "/RDSfiles", output.dir), '/datasets.before.batch.removal.rds'))
  # DimPlot(merged.data, reduction = "tsne", group.by = "datasetID", label = FALSE, repel = TRUE, split.by = "datasetID") + theme(plot.title = element_text(hjust = 0.5))
  # DimPlot(merged.data, reduction = "umap", group.by = "datasetID", label = FALSE, repel = TRUE, split.by = "datasetID") + theme(plot.title = element_text(hjust = 0.5))

  # perform batch removal
  if(Step2_Quality_Control.RemoveBatches){
      print('Performing batch removal...')
      ifnb.list <- SplitObject(merged.data, split.by = "datasetID")
      # print(ifnb.list)
      ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
        x <- NormalizeData(x, verbose = FALSE)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
      })
    
      # perform integration
      merged.data.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = PCs, anchor.features = nfeatures)
      merged.data.combined <- IntegrateData(anchorset = merged.data.anchors, dims = PCs)
    
      # perform an integrated analysis
      DefaultAssay(merged.data.combined) <- "integrated"
    
      # run the standard workflow for visualization and clustering
      merged.data.combined <- ScaleData(merged.data.combined, verbose = FALSE, vars.to.regress = vars.to.regress)
      merged.data.combined <- RunPCA(merged.data.combined, npcs = ndims, features = VariableFeatures(object = merged.data.combined), verbose = FALSE)
    
      pdf(paste0(output.dir, '/datasets.after.batch.removal.ElbowPlot.pdf'), width = 5, height = 5)
       print(ElbowPlot(merged.data.combined, ndims = ndims))
      dev.off()

      png(paste0(output.dir, '/datasets.after.batch.removal.ElbowPlot.png'), width = 500, height = 500)
       print(ElbowPlot(merged.data.combined, ndims = ndims))
      dev.off()

      # FindNeighbors
      merged.data.combined <- FindNeighbors(merged.data.combined,
                                            dims = PCs,
                                            k.param=n.neighbors,
                                            force.recalc = TRUE)
    
      # Visualization
      merged.data.combined <- RunTSNE(object = merged.data.combined, dims=PCs, check_duplicates = FALSE)
      merged.data.combined <- RunUMAP(merged.data.combined, reduction = "pca", dims=PCs,
                                      n.neighbors=n.neighbors, min.dist=1, seed.use=2023L)
    
      pdf(paste0(output.dir, '/datasets.after.batch.removal_tsne_datasetID.pdf'), width = (5+max_length*0.1), height = 5)
        print(DimPlot(merged.data.combined, reduction = "tsne", group.by = "datasetID", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      png(paste0(output.dir, '/datasets.after.batch.removal_tsne_datasetID.png'), width = (500+max_length*10), height = 500)
        print(DimPlot(merged.data.combined, reduction = "tsne", group.by = "datasetID", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()
    
      pdf(paste0(output.dir, '/datasets.after.batch.removal_umap_datasetID.pdf'), width = (5+max_length*0.1), height = 5)
        print(DimPlot(merged.data.combined, reduction = "umap", group.by = "datasetID", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      png(paste0(output.dir, '/datasets.after.batch.removal_umap_datasetID.png'), width = (500+max_length*10), height = 500)
        print(DimPlot(merged.data.combined, reduction = "umap", group.by = "datasetID", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()
  }else{
      merged.data.combined <- merged.data
  }  

  # Doublet
  color_vector <- c('#bdbdbd', '#a50f15')
  names(color_vector) <- c('Singlet', 'Doublet')
  saveRDS(merged.data.combined, file=paste0(gsub("/Step2.Quality_control/$", "/RDSfiles", output.dir), '/multipleDataSetsAfterQC.rds'))
  pdf(paste0(output.dir, '/findDoublets.pdf'),width=6,height=5)
    print(DimPlot(merged.data.combined, reduction = "umap", label = FALSE, pt.size = 0.5, cols = color_vector, group.by='doublet', raster = FALSE)+
      ggtitle(paste(as.character(length(which(merged.data.combined@meta.data$doublet == 'Doublet'))), 'doublets and',
                    as.character(length(which(merged.data.combined@meta.data$doublet == 'Singlet'))), 'singlets')))
  dev.off()

  png(paste0(output.dir, '/findDoublets.png'),width=600,height=500)
    print(DimPlot(merged.data.combined, reduction = "umap", label = FALSE, pt.size = 0.5, cols = color_vector, group.by='doublet', raster = FALSE)+
      ggtitle(paste(as.character(length(which(merged.data.combined@meta.data$doublet == 'Doublet'))), 'doublets and',
                    as.character(length(which(merged.data.combined@meta.data$doublet == 'Singlet'))), 'singlets')))
  dev.off()  

  if(Step2_Quality_Control.RemoveDoublets){
    print(paste('remove', as.character(length(which(merged.data.combined@meta.data$doublet == 'Doublet'))), 'doublets.'))
    merged.data.combined <- subset(merged.data.combined, subset = doublet=='Singlet') 
  }  

  saveRDS(merged.data.combined, file=paste0(gsub("/Step2.Quality_control/$", "/RDSfiles", output.dir), '/multipleDataSetsAfterQC.rds'))
  return(merged.data.combined)
}

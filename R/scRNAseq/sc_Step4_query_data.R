#' Map Data from Query Object to Reference Object
#'
#' Function \code{mapDataToRef} performs data mapping from a query Seurat object to a reference Seurat object.
#'
#' @param ref_object A reference Seurat object containing the reference data.
#' @param ref_labels A vector or factor specifying the labels for the reference data.
#' @param query_object A query Seurat object containing the data to be mapped to the reference.
#' @param PCs A vector specifying the principal components to be used for mapping.
#' @param output.dir The directory where the mapping results and plots will be saved.
#'
#' @details
#' This function maps data from a query Seurat object to a reference Seurat object using the specified principal components (PCs). It identifies transfer anchors between the reference and query objects, transfers data, and adds predictions to the query object. The mapping results are visualized and saved as plots in the output directory. Note that ref_object and query_object should have sample type feature names. Refer to https://satijalab.org/seurat/archive/v3.0/integration.html
#' 
#' @return mapDataToRef returns the query Seurat object with added prediction scores and IDs.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#'
#' @export

mapDataToRef = function(ref_object=NULL,
                        ref_labels=NULL,
                        query_object=NULL,
                        PCs = NULL,
                        output.dir = NULL){
  ref.anchors <- FindTransferAnchors(reference = ref_object,
                                     query = query_object,
                                     dims = PCs)
  predictions <- TransferData(anchorset = ref.anchors,
                              refdata = ref_labels,
                              dims = PCs)
  query_object <- AddMetaData(query_object,
                              metadata = predictions)

  pdf(paste0(output.dir,'/predicted.id.pdf'), width = 8, height = 6)
    print(DimPlot(query_object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='predicted.id', raster = FALSE))
  dev.off()
    
  png(paste0(output.dir,'/predicted.id.png'), width = 800, height = 600)
    print(DimPlot(query_object, reduction = "umap", pt.size=0.5, label = FALSE, label.size=7, group.by='predicted.id', raster = FALSE))
  dev.off()  

  dir.create(paste0(output.dir,'/SeuratScores/'))
  for(i in paste0('prediction.score.', unique(make.names(ref_labels)))){
    if(!all(query_object@meta.data[,which(colnames(query_object@meta.data)==i)]==query_object@meta.data[1,which(colnames(query_object@meta.data)==i)])){
        pdf(paste0(output.dir,'/SeuratScores/', make.names(as.character(i)),'.pdf'),width=8,height=8)
          print(FeaturePlot(query_object, features = make.names(as.character(i)),cols=c("lightgrey",'red'),raster=FALSE))
        dev.off()
          
        png(paste0(output.dir,'/SeuratScores/', make.names(as.character(i)),'.png'),width=800,height=800)
          print(FeaturePlot(query_object, features = make.names(as.character(i)),cols=c("lightgrey",'red'),raster=FALSE))
        dev.off()  
    }  
  }

  saveRDS(query_object, file = paste0(gsub("/Step4.Identify_Cell_Types/$", "/RDSfiles", output.dir),'/sc_object.With.PredictionScore.And.PredictionID.RDS'))
  return(query_object)
}

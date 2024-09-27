#' Dimension reduction by phateR.
#'
#' Function \code{run_phateR} performs dimension reduction using phateR algorithm.
#'
#' @param sc_object A Seurat object containing single-cell RNA-seq data.
#' @param output.dir The directory where the dimension reduction plot will be saved.
#' @param pythonPath The path to the Python executable if required for running phateR.
#' @param phate.knn The number of nearest neighbors to consider in the phateR algorithm. Default 50.
#' @param phate.npca The number of principal components to use in the phateR algorithm. Default 20.
#' @param phate.t The t-value for the phateR algorithm, which controls the level of exploration. Default 10.
#' @param phate.ndim The number of dimensions for the output embedding in the phateR algorithm. Default 2.
#'
#' @details
#' This function performs dimension reduction using the phateR algorithm on a Seurat object.
#' It generates a dimension reduction plot and saves it as a PDF and PNG file in the specified output directory.
#' 
#' @return run_phateR returns the dimension reduction plot in both PDF and PNG formats.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#'
#' @export

run_phateR = function(sc_object = NULL,
                      output.dir = NULL,
                      pythonPath=NULL,
                      phate.knn = 50,
                      phate.npca = 20,
                      phate.t = 10,
                      phate.ndim = 2){
  if(is.null(pythonPath)==FALSE){reticulate::use_python(pythonPath)}

  phate.object <- phate(subset(GetAssayData(object = sc_object, slot = "scale.data"),
                               rownames(GetAssayData(object = sc_object, slot = "scale.data"))%in% VariableFeatures(object = sc_object))%>%t(),
                        knn = phate.knn,
                        npca = phate.npca,
                        t = phate.t,
                        ndim = phate.ndim)
  saveRDS(phate.object, file = paste0(gsub("/Step5.Visualization/$", "/RDSfiles", output.dir),'/phate.object.RDS'))
    
  pdf(paste0(output.dir, '/phateR.pdf'), width = 8, height = 8)
   print(ggplot(phate.object) + geom_point(aes(x=PHATE1, y=PHATE2, color=Idents(sc_object)))+mytheme)
  dev.off()
    
  png(paste0(output.dir, '/phateR.png'), width = 800, height = 800)
   print(ggplot(phate.object) + geom_point(aes(x=PHATE1, y=PHATE2, color=Idents(sc_object)))+mytheme)
  dev.off()  
}
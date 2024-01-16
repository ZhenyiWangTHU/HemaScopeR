#' Prepare Data for scVelo Analysis
#'
#' Function \code{prepareDataForScvelo} prepares the data for scVelo analysis by formatting and saving the necessary input files.
#' 
#' @param sc_object A Seurat object containing the single-cell RNA-seq data.
#' @param loom.files.path A character vector specifying the path(s) to the loom files for scVelo analysis.
#' @param scvelo.reduction A character specifying the reduction method used for scVelo analysis (default is 'pca').
#' @param scvelo.column A character specifying the column in the Seurat object metadata containing cell types.
#' @param output.dir A character specifying the directory where the scVelo input files will be saved.
#' 
#' @details
#' This function takes a Seurat object, cell embeddings, and cell type information as input and performs the following steps:
#' 
#' 1. Formats the Seurat object metadata, including cell types and sample names.
#' 2. Extracts the spliced, unspliced, and ambiguous count matrices from the Seurat object.
#' 3. Combines the metadata and cell embeddings.
#' 4. Writes the necessary input files for scVelo analysis, including cell embeddings, count matrices, and metadata.
#' 
#' @return None. This function saves the required input files for scVelo analysis in the output.dir.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#' 
#' @export

prepareDataForScvelo = function(sc_object = NULL,
                                loom.files.path = NULL,
                                scvelo.reduction = 'pca',
                                scvelo.column = NULL,
                                output.dir = NULL){
  sc_object <- RenameCells(sc_object, new.names = gsub('-\\d+$', '', colnames(sc_object)))  
  cell.embeddings <- Embeddings(object = sc_object, reduction = scvelo.reduction)[, 1:2] %>% as.data.frame()
  if(length(loom.files.path)>1){
        splicedmat.list <- list()
        umat.list <- list()
        amat.list <- list()
        scvelo_combined.obs.list <- list()
        loom.files.path.length <- length(loom.files.path)
      for(i in 1:loom.files.path.length){
        ldat.temp <- ReadVelocity(file = loom.files.path[i])
        ldat.seurat_object.temp <- as.Seurat(x = ldat.temp)
        ldat.seurat_object.temp <- RenameCells(ldat.seurat_object.temp, new.names = gsub(':', '_', colnames(ldat.seurat_object.temp)))  
        ldat.seurat_object.temp <- RenameCells(ldat.seurat_object.temp, new.names = gsub('x$', '', colnames(ldat.seurat_object.temp)))    
        ldat.seurat_object.temp@meta.data$samplenames <- colnames(ldat.seurat_object.temp)
        ldat.seurat_object.temp <- subset(ldat.seurat_object.temp, subset = samplenames%in%colnames(sc_object))
        ldat.seurat_object.temp@meta.data$cellTypes <- plyr::mapvalues(rownames(ldat.seurat_object.temp@meta.data),
                                                                       from=rownames(sc_object@meta.data),
                                                                       to=sc_object@meta.data[,scvelo.column],
                                                                       warn_missing=FALSE)
        splicedmat.temp <- ldat.seurat_object.temp@assays$spliced@counts
        umat.temp <- ldat.seurat_object.temp@assays$unspliced@counts
        amat.temp <- ldat.seurat_object.temp@assays$ambiguous@counts  
        scvelo_combined.obs.temp <- data.frame(cellTypes = ldat.seurat_object.temp@meta.data$cellTypes, 
                                               samplenames = ldat.seurat_object.temp@meta.data$samplenames)
        rownames(scvelo_combined.obs.temp) <- scvelo_combined.obs.temp$samplenames  
        #colnames(scvelo_combined.obs.temp) <- c("cellTypes", "samplenames")  
        # splicedmat.list <- c(splicedmat.list, splicedmat.temp)
        # umat.list <- c(umat.list, umat.temp)
        # amat.list <- c(amat.list, amat.temp)
        # scvelo_combined.obs.list <- c(scvelo_combined.obs.list, scvelo_combined.obs.temp) 
          
        splicedmat.list[[i]] <- splicedmat.temp
        umat.list[[i]] <- umat.temp
        amat.list[[i]] <- amat.temp
        scvelo_combined.obs.list[[i]] <- scvelo_combined.obs.temp
      }
        #splicedmat <- cbind_all(splicedmat.list)
        #umat <- cbind_all(umat.list)
        #amat <- cbind_all(amat.list)
        #scvelo_combined.obs <- rbind_all(scvelo_combined.obs.list)
        splicedmat <- do.call(cbind, splicedmat.list)
        umat <- do.call(cbind, umat.list)
        amat <- do.call(cbind, amat.list)
        scvelo_combined.obs <- do.call(rbind, scvelo_combined.obs.list)
        cell.embeddings$cellNames <- factor(rownames(cell.embeddings), levels = colnames(splicedmat))
        cell.embeddings <- cell.embeddings[order(cell.embeddings$cellNames, decreasing = FALSE),]
        cell.embeddings <- cell.embeddings[,-3]
        write.csv(cell.embeddings, file=paste0(output.dir, '/cell.embeddings.csv'), quote=FALSE, row.names = FALSE)
        writeMM(t(splicedmat), file = paste0(output.dir, "/splicedmat.mtx"))
        writeMM(t(umat), file = paste0(output.dir, "/umat.mtx"))
        writeMM(t(amat), file = paste0(output.dir, "/amat.mtx"))
        write.csv(scvelo_combined.obs, file = paste0(output.dir, "/scvelo_combined.obs.csv"), row.names = FALSE)
        write.csv(rownames(splicedmat), file= paste0(output.dir, '/geneInfo.csv'), quote=FALSE,row.names = FALSE)
  }else if(length(loom.files.path)==1){
        ldat <- ReadVelocity(file = loom.files.path)
        ldat.seurat_object <- as.Seurat(x = ldat)
        ldat.seurat_object <- RenameCells(ldat.seurat_object, new.names = gsub(':', '_', colnames(ldat.seurat_object)))  
        ldat.seurat_object <- RenameCells(ldat.seurat_object, new.names = gsub('x$', '', colnames(ldat.seurat_object)))    
        ldat.seurat_object@meta.data$samplenames <- colnames(ldat.seurat_object)
        ldat.seurat_object <- subset(ldat.seurat_object, subset = samplenames%in%colnames(sc_object))
        ldat.seurat_object@meta.data$cellTypes <- plyr::mapvalues(rownames(ldat.seurat_object@meta.data),
                                                                  from=rownames(sc_object@meta.data),
                                                                  to=sc_object@meta.data[,scvelo.column],
                                                                  warn_missing=FALSE)
        splicedmat <- ldat.seurat_object@assays$spliced@counts
        umat <- ldat.seurat_object@assays$unspliced@counts
        amat <- ldat.seurat_object@assays$ambiguous@counts
        scvelo_combined.obs <- as.data.frame(cbind(ldat.seurat_object@meta.data$cellTypes, 
                                                   ldat.seurat_object@meta.data$samplenames))
        colnames(scvelo_combined.obs) <- c("cellTypes", "samplenames")
        cell.embeddings$cellNames <- factor(rownames(cell.embeddings), levels = colnames(splicedmat))
        cell.embeddings <- cell.embeddings[order(cell.embeddings$cellNames, decreasing = FALSE),]
        cell.embeddings <- cell.embeddings[,-3]
        write.csv(cell.embeddings, file=paste0(output.dir, '/cell.embeddings.csv'), quote=FALSE, row.names = FALSE)
        writeMM(t(splicedmat), file = paste0(output.dir, "/splicedmat.mtx"))
        writeMM(t(umat), file = paste0(output.dir, "/umat.mtx"))
        writeMM(t(amat), file = paste0(output.dir, "/amat.mtx"))
        write.csv(scvelo_combined.obs, file = paste0(output.dir, "/scvelo_combined.obs.csv"), row.names = FALSE)
        write.csv(rownames(splicedmat), file= paste0(output.dir, '/geneInfo.csv'), quote=FALSE,row.names = FALSE)
  }else{
        print('Please input the directory of loom files.')
  }  
 }
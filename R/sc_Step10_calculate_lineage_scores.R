#' Calculate Lineage Scores
#' 
#' This function calculates lineage scores based on provided expression data and cell type annotations. It also generates visualizations of lineage scores and gene expression patterns.
#' 
#' @param expression_matrix A data frame or matrix representing the expression data, where rows are cells and columns are genes.
#' @param cellTypes A character vector specifying cell type annotations for each cell. e.g. c("HSC","HSC","HSC","MPP1","MPP2","MPP2","MPP2","MPP3","MPP4","MEP","GMP", "MKP","CMP","CLP","LSKlow", ...)
#' @param cellTypes_orders A character vector specifying the order of cell types for plotting. e.g. c("HSC","MPP1","MPP2","MPP3","MPP4","MEP","GMP", "MKP","CMP","CLP","LSKlow")
#' @param cellTypes_colors A character vector specifying the colors to use for cell type groups. e.g. c("HSC" = '#006d2c',"MPP1" = '#4292c6',"MPP2"= '#810f7c',"MPP3" = '#fec44f',"MPP4" = '#dd3497',"MEP" = '#ef6548',"GMP" = '#993404', "MKP" = '#8073ac',"CMP" = '#8c510a',"CLP" = '#b2df8a',"LSKlow" = '#d53e4f').
#' @param groups A character vector specifying groups or clusters within each cell type.
#' @param groups_orders A character vector specifying the order of groups or clusters for plotting.
#' @param groups_colors A character vector specifying the colors to use for group or cluster annotations. e.g. c('group1'='#d73027','group2'='#2171b5')
#' @param lineage.genelist A list of gene sets representing lineage markers.
#' @param lineage.names A character vector specifying the names of the lineages.
#' @param Org A character specifying the organism ('mmu' for mouse, 'hsa' for human).
#' @param output.dir The path to the directory where the resulting visualizations and scores will be saved.
#' 
#' @details
#' This function calculates lineage scores for specified gene sets based on the provided expression data. It then generates a heatmap of lineage scores and a heatmap of gene expression patterns. You can customize the order and colors of cell types, groups, and clusters for visualization. The resulting visualizations and lineage scores are saved to the specified output directory.
#' 
#' @return None. This function generates and saves visualizations and lineage scores as PDF and CSV files.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#' 
#' @export

lineageScores = function(expression_matrix = NULL,
                         cellTypes = NULL,
                         cellTypes_orders = NULL,
                         cellTypes_colors = NULL,
                         groups = NULL,
                         groups_orders = NULL,
                         groups_colors = NULL,
                         lineage.genelist = NULL,
                         lineage.names= NULL,
                         Org = NULL,
                         output.dir = NULL){
  names(lineage.genelist) <- lineage.names
  cellTypes_groups_orders <- c()
  for (i in cellTypes_orders) {
    for (j in groups_orders) {
      combined_str <- paste(i, j, sep = '_')
      cellTypes_groups_orders <- c(cellTypes_groups_orders, combined_str)
    }
  }
  cellTypes_groups <- paste(cellTypes, groups, sep='_')
  cellorders <- colnames(expression_matrix)

  lineages_signatures <- c()
  for(i in 1:length(lineage.genelist)){
    lineages_signatures <- c(lineages_signatures, lineage.genelist[[i]])
  }

  # check the feature names
  #data("genecode_geneSymbolandEnsembleID")
  if(Org == 'mmu'){
     load("../data/mouseGeneSymbolandEnsembleID.rdata")
     genecode_geneSymbolandEnsembleID <- mouseGeneSymbolandEnsembleID
  }else if(Org == 'hsa'){
     load("../data/humanGeneSymbolandEnsembleID.rdata")
     genecode_geneSymbolandEnsembleID <- humanGeneSymbolandEnsembleID
  }
        
  if(all(grepl('^ENSMUSG|^ENSG', rownames(expression_matrix)))){
    # if the features are ensemble id, we need convert the rownames from ensembleID to gene name,
    # e.g. 'ENSG00000102145' -> 'Gata1'.
    lineages_signatures <- unique(lineages_signatures)
    lineages_signatures_ensembleId <- mapvalues(x = lineages_signatures,
                                                from = genecode_geneSymbolandEnsembleID$geneName,
                                                to= genecode_geneSymbolandEnsembleID$ensemblIDNoDot,
                                                warn_missing = FALSE)
    expression_matrix <- subset(expression_matrix, rownames(expression_matrix)%in%lineages_signatures_ensembleId)
    rownames(expression_matrix) <- mapvalues(x = rownames(expression_matrix),
                                             from = genecode_geneSymbolandEnsembleID$ensemblIDNoDot,
                                             to= genecode_geneSymbolandEnsembleID$geneName,
                                             warn_missing = FALSE)
    }else{
      lineages_signatures <- unique(lineages_signatures)
      expression_matrix <- subset(expression_matrix, rownames(expression_matrix)%in%lineages_signatures)
    }

  # lineage score matrix
  lineage_signatures_scores <- matrix(data = 0,
                                      nrow=ncol(expression_matrix),
                                      ncol=length(lineage.genelist))
  rownames(lineage_signatures_scores) <- colnames(expression_matrix)
  colnames(lineage_signatures_scores) <- names(lineage.genelist)

  # calculate scores
  for(i in 1:length(lineage.genelist)){
    sub_matrix <- subset(expression_matrix, rownames(expression_matrix) %in% unique(lineage.genelist[[i]]))
    for(j in 1:nrow(lineage_signatures_scores)){
      lineage_signatures_scores[j,i] <- sum(sub_matrix[,j])/nrow(sub_matrix)
    }
  }
  write.csv(lineage_signatures_scores, paste0(output.dir,'/lineage_signatures_scores.csv'), quote =FALSE)

  # visualization
  # plot the heatmap for the lineage score matrix-----------------------------------------------------------
  matrixForheatmap <- lineage_signatures_scores
  # order the matrix for heatmap
  i_zero <- c()
  for(i in 1:length(cellTypes_groups_orders)){
      if(length(which(cellTypes_groups==cellTypes_groups_orders[i]))==0){
          i_zero <- c(i_zero, i)
      }
  }

  if(length(i_zero)>0){
      cellTypes_groups_orders <- cellTypes_groups_orders[-i_zero]
  }  
 
  if(length(which(cellTypes_groups==cellTypes_groups_orders[1]))>1){
    matrixForheatmap_ordered <- matrixForheatmap[which(cellTypes_groups==cellTypes_groups_orders[1]), ]  
    matrixForheatmap_ordered <- matrixForheatmap_ordered[hclust(dist(matrixForheatmap_ordered))$order, ]
    }else if(length(which(cellTypes_groups==cellTypes_groups_orders[1]))==1){
      matrixForheatmap_ordered <- matrixForheatmap[which(cellTypes_groups==cellTypes_groups_orders[1]), ] 
      matrixForheatmap_ordered <- as.data.frame(matrix(matrixForheatmap_ordered, nrow = 1, ncol = length(matrixForheatmap_ordered)))
      rownames(matrixForheatmap_ordered) <- rownames(matrixForheatmap)[which(cellTypes_groups==cellTypes_groups_orders[1])]
      colnames(matrixForheatmap_ordered) <- colnames(matrixForheatmap)
  }  

  for(i in cellTypes_groups_orders[2:length(cellTypes_groups_orders)]){
    matrixForheatmap_temp <- matrixForheatmap[which(cellTypes_groups==i), ]
    if(length(which(cellTypes_groups==i))>1){
      matrixForheatmap_temp <- matrixForheatmap_temp[hclust(dist(matrixForheatmap_temp))$order,]
      matrixForheatmap_ordered <- rbind(matrixForheatmap_ordered, matrixForheatmap_temp)
    }else if(length(which(cellTypes_groups==i))==1){
      matrixForheatmap_temp <- as.data.frame(matrix(matrixForheatmap_temp, nrow = 1, ncol = length(matrixForheatmap_temp)))
      rownames(matrixForheatmap_temp) <- rownames(matrixForheatmap)[which(cellTypes_groups==i)]
      colnames(matrixForheatmap_temp) <- colnames(matrixForheatmap_ordered)
      matrixForheatmap_ordered <- rbind(matrixForheatmap_ordered, matrixForheatmap_temp)
    }
  }
 
  matrixForheatmap_ordered <- t(matrixForheatmap_ordered)

  annotation_col_C_andcluster = data.frame(cellTypes = factor(plyr::mapvalues(colnames(matrixForheatmap_ordered),
                                                                              from = cellorders,
                                                                              to = cellTypes),
                                                              levels = cellTypes_orders),
                                           group = factor(plyr::mapvalues(colnames(matrixForheatmap_ordered),
                                                                          from = cellorders,
                                                                          to = groups),
                                                          levels = groups_orders))

  rownames(annotation_col_C_andcluster) = colnames(matrixForheatmap_ordered)
    
  if((is.null(cellTypes_colors) == FALSE)&(is.null(cellTypes_orders) == FALSE)&(is.null(groups_colors) == FALSE)&(is.null(groups_orders) == FALSE)){
      
    if((length(cellTypes_colors)==length(cellTypes_orders))&(length(groups_colors)==length(groups_orders))){
      names(cellTypes_colors) <- levels(annotation_col_C_andcluster$cellTypes)
      names(groups_colors) <- levels(annotation_col_C_andcluster$group)  
      ann_colors_C = list(cellTypes = cellTypes_colors,group = groups_colors)
    }else{ann_colors_C = NULL}
      
  }else{ann_colors_C = NULL}

  matrixForheatmap_ordered_scaled <- scale(matrixForheatmap_ordered,
                                           center = T,
                                           scale = T)

  if(length(which(colSums(is.na(matrixForheatmap_ordered_scaled))>0))>0){
    annotation_col_C_andcluster <- annotation_col_C_andcluster[-which(colSums(is.na(matrixForheatmap_ordered_scaled))>0), ]
    matrixForheatmap_ordered_scaled <- matrixForheatmap_ordered_scaled[,-which(colSums(is.na(matrixForheatmap_ordered_scaled))>0)]
  }

    pheatmap::pheatmap(as.matrix(matrixForheatmap_ordered_scaled),
                               cluster_rows = F,
                               cluster_cols =F,
                               scale = "none",
                               legend_breaks = ceiling(seq(min(matrixForheatmap_ordered_scaled),
                                                          max(matrixForheatmap_ordered_scaled),0.01)),
                               color = colorRampPalette(colors = c("turquoise1","black","gold"))(length(seq(min(matrixForheatmap_ordered_scaled),
                                                                                                            max(matrixForheatmap_ordered_scaled),0.01))),
                               show_colnames = F,
                               show_rownames = T,
                               annotation_col  = annotation_col_C_andcluster,
                               annotation_colors = ann_colors_C,
                               #fontsize =0.5,
                               filename=paste0(output.dir,'/lineageScoresHeatmap.pdf'),
                               height=0.4*nrow(matrixForheatmap_ordered_scaled),
                               annotation_legend = FALSE,
                               weight=10
                               #cellheight=0.3,
                               #cellwidth=0.01
    )
    
    pheatmap::pheatmap(as.matrix(matrixForheatmap_ordered_scaled),
                               cluster_rows = F,
                               cluster_cols =F,
                               scale = "none",
                               legend_breaks = ceiling(seq(min(matrixForheatmap_ordered_scaled),
                                                          max(matrixForheatmap_ordered_scaled),0.01)),
                               color = colorRampPalette(colors = c("turquoise1","black","gold"))(length(seq(min(matrixForheatmap_ordered_scaled),
                                                                                                            max(matrixForheatmap_ordered_scaled),0.01))),
                               show_colnames = F,
                               show_rownames = T,
                               annotation_col  = annotation_col_C_andcluster,
                               annotation_colors = ann_colors_C,
                               #fontsize =0.5,
                               filename=paste0(output.dir,'/lineageScoresHeatmap.png'),
                               height=0.4*nrow(matrixForheatmap_ordered_scaled),
                               annotation_legend = FALSE,
                               weight=10
                               #cellheight=0.3,
                               #cellwidth=0.01
    )

  # plot the heatmap for the gene expression matrix------------------------------------------------------------------------------------------
  expression_matrix_ordered <- expression_matrix[order(factor(rownames(expression_matrix),
                                                          levels=factor(lineages_signatures, levels=lineages_signatures))),
                                                 order(factor(colnames(expression_matrix),
                                                          levels=factor(colnames(matrixForheatmap_ordered),
                                                                        levels=colnames(matrixForheatmap_ordered))))]
  expression_matrix_ordered <- scale(expression_matrix_ordered, center = T, scale = T)
  expression_matrix_ordered[expression_matrix_ordered > 2] <- 2
  expression_matrix_ordered[expression_matrix_ordered < (-2)] <- (-2)

  # cluster and order cells in each block and then merge them
  expression_matrix_ordered.temp <- subset(expression_matrix_ordered, rownames(expression_matrix_ordered)%in%lineage.genelist[[1]])
  expression_matrix_ordered.temp <- expression_matrix_ordered.temp[hclust(dist(expression_matrix_ordered.temp))$order,]
  for(i in 2:length(lineage.genelist)){
    sub_matrix_ordered <- subset(expression_matrix_ordered, rownames(expression_matrix_ordered)%in%lineage.genelist[[i]])
    sub_matrix_ordered <- sub_matrix_ordered[hclust(dist(sub_matrix_ordered))$order,]
    expression_matrix_ordered.temp <- rbind(expression_matrix_ordered.temp, sub_matrix_ordered)
  }
  expression_matrix_ordered <- expression_matrix_ordered.temp

  if(length(which(colSums(is.na(expression_matrix_ordered))>0))>0){
    expression_matrix_ordered <- expression_matrix_ordered[,-which(colSums(is.na(expression_matrix_ordered))>0)]
  }

    pheatmap::pheatmap(as.matrix(expression_matrix_ordered),
           cluster_rows = F,
           cluster_cols =F,
           scale = "none" ,
           legend_breaks= ceiling(seq(min(expression_matrix_ordered),
                                      max(expression_matrix_ordered),0.01)),
           color = colorRampPalette(colors = c("#377eb8","#deebf7","#e41a1c"))(length(seq(min(expression_matrix_ordered),
                                                                                          max(expression_matrix_ordered),0.01))),
           show_colnames = F,
           show_rownames = T,
           annotation_col  = annotation_col_C_andcluster,
           annotation_colors = ann_colors_C,
           #fontsize =15,
           filename=paste0(output.dir,'/lineageGenesHeatmap.pdf'),
           height=0.2*nrow(expression_matrix_ordered),
           weight=10
           #cellheight=30,
           #cellwidth=0.5
    )
    
    pheatmap::pheatmap(as.matrix(expression_matrix_ordered),
           cluster_rows = F,
           cluster_cols =F,
           scale = "none" ,
           legend_breaks= ceiling(seq(min(expression_matrix_ordered),
                                      max(expression_matrix_ordered),0.01)),
           color = colorRampPalette(colors = c("#377eb8","#deebf7","#e41a1c"))(length(seq(min(expression_matrix_ordered),
                                                                                          max(expression_matrix_ordered),0.01))),
           show_colnames = F,
           show_rownames = T,
           annotation_col  = annotation_col_C_andcluster,
           annotation_colors = ann_colors_C,
           #fontsize =15,
           filename=paste0(output.dir,'/lineageGenesHeatmap.png'),
           height=0.2*nrow(expression_matrix_ordered),
           weight=10
           #cellheight=30,
           #cellwidth=0.5
    )

}

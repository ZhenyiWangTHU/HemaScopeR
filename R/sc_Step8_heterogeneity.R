#' Measure Cell Heterogeneity Using Spearman Correlation
#' 
#' This function calculates the cell heterogeneity based on Spearman correlation using the input expression matrix and cell type group annotations.
#' 
#' @param expression_matrix A numeric matrix representing the expression data, where rows are genes and columns are cells.
#' @param cell_types_groups A data frame specifying cell type annotations for each cell, including cell type labels and group information.
#' @param output.dir The path to the directory where the resulting plots will be saved.
#' 
#' @details
#' This function quantifies cell heterogeneity by computing Spearman correlation coefficients between cells within the same cell type groups. It generates plots visualizing the heterogeneity measurements for different cell types and groups. The `expression_matrix` should be appropriately preprocessed and filtered before using this function.
#' 
#' @return None. This function generates and saves plots to the specified output directory.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#' 
#' @export

heterogeneity = function(expression_matrix = NULL,
                         cell_types_groups = NULL,
                         cellTypeOrders = NULL,
                         output.dir = NULL){
  normalizedData_rowMeans <- rowMeans(expression_matrix)
  normalizedData_rowMeans <- normalizedData_rowMeans[which(normalizedData_rowMeans > 0)] %>% names()
  expression_matrix <- subset(expression_matrix, rownames(expression_matrix) %in% normalizedData_rowMeans)
  expression_matrix <- t(expression_matrix)
  cell_types_groups <- as.data.frame(cell_types_groups)
  if(length(unique(cell_types_groups[,2])) == 2){
      cell_types_groups$clusterAndgroup <- paste(cell_types_groups[,1], cell_types_groups[,2], sep='_')
      distance_values <- c()
      cell_types_groups.vector <- c()
      cellTypes <- c()
      group <- c()
      cell_types_groups.unique <- unique(cell_types_groups)
      for(i in 1:nrow(cell_types_groups.unique)){
        sub_expression_matrix <- expression_matrix[which(cell_types_groups$clusterAndgroup == cell_types_groups.unique$clusterAndgroup[i]), ]
        # spearman correlation
        sub_expression_matrix <- t(sub_expression_matrix)
        distance_values_temp <- cor(x=sub_expression_matrix, y=NULL, method="spearman")
        distance_values_temp <- as.numeric(distance_values_temp[upper.tri(distance_values_temp)])

        cell_types_groups_temp <- as.character(rep(paste(cell_types_groups.unique[i,1],
                                                         cell_types_groups.unique[i,2],
                                                         sep='_'),
                                                times = length(distance_values_temp)))
        cellTypes_temp <- as.character(rep(cell_types_groups.unique[i,1],
                                           times = length(distance_values_temp)))
        group_temp <- as.character(rep(cell_types_groups.unique[i,2],
                                       times = length(distance_values_temp)))

        distance_values <- c(distance_values, distance_values_temp)
        cell_types_groups.vector <- c(cell_types_groups.vector, cell_types_groups_temp)
        cellTypes <- c(cellTypes, cellTypes_temp)
        group <- c(group, group_temp)
      }

      distance_stat <- data.frame(x=cell_types_groups.vector,
                                  y=distance_values,
                                  stringsAsFactors=FALSE)
      colnames(distance_stat) <- c('orig.ident', 'distance')
      distance_stat$group <- group
      distance_stat$cellTypes <- cellTypes
      distance_stat$cellTypes <- factor(distance_stat$cellTypes, levels=cellTypeOrders)

      # get the longest datasetID
      datasetID_lengths <- sapply(unique(group), nchar)
      max_length <- max(datasetID_lengths)
      cellTypes_length <- length(unique(cellTypes))
      
      pdf(file = paste0(output.dir,'/Heterogeneity_spearmanCorrelation.pdf'),
          width = (8+max_length*0.1+cellTypes_length*0.5),
          height = 8)
          print(ggplot(data=distance_stat, aes(x=cellTypes,
                                               y=distance,
                                               fill=group),
                       color='#ffffff')+
                stat_boxplot(geom ='errorbar', width = 0.6,position =position_dodge(0.8))+
                geom_boxplot(outlier.shape = NA, position=position_dodge(0.8))+
                scale_fill_manual(values = c(alpha('#e41a1c',1), alpha('#377eb8',1)))+
                xlab(label='cell types')+ylab(label='Spearman correlation')+
                ggpubr::stat_compare_means(aes(group=group), method = 'wilcox.test', label = "p.signif")+
                ylim((min(distance_stat$distance)-0.05),(max(distance_stat$distance)+0.05))+mytheme+
                  theme(axis.text.x = element_text(size = 20,
                                                   family = "sans",
                                                   color = "black",
                                                   face = "bold",
                                                   vjust = 0,
                                                   hjust = 0.5))
                  )
      dev.off()

      png(file = paste0(output.dir,'/Heterogeneity_spearmanCorrelation.png'),
          width = (800 + max_length*10+cellTypes_length*50),
          height = 800)
          print(ggplot(data=distance_stat, aes(x=cellTypes,
                                               y=distance,
                                               fill=group),
                       color='#ffffff')+
                stat_boxplot(geom ='errorbar', width = 0.6,position =position_dodge(0.8))+
                geom_boxplot(outlier.shape = NA, position=position_dodge(0.8))+
                scale_fill_manual(values = c(alpha('#e41a1c',1), alpha('#377eb8',1)))+
                xlab(label='cell types')+ylab(label='Spearman correlation')+
                ggpubr::stat_compare_means(aes(group=group), method = 'wilcox.test', label = "p.signif")+
                ylim((min(distance_stat$distance)-0.05),(max(distance_stat$distance)+0.05))+mytheme+
                  theme(axis.text.x = element_text(size = 20,
                                                   family = "sans",
                                                   color = "black",
                                                   face = "bold",
                                                   vjust = 0,
                                                   hjust = 0.5))
              )
      dev.off()  
  }else{
      cell_types_groups$clusterAndgroup <- paste(cell_types_groups[,1], cell_types_groups[,2], sep='_')
      distance_values <- c()
      cell_types_groups.vector <- c()
      cellTypes <- c()
      group <- c()
      cell_types_groups.unique <- unique(cell_types_groups)
      for(i in 1:nrow(cell_types_groups.unique)){
        sub_expression_matrix <- expression_matrix[which(cell_types_groups$clusterAndgroup == cell_types_groups.unique$clusterAndgroup[i]), ]
        # spearman correlation
        sub_expression_matrix <- t(sub_expression_matrix)
        distance_values_temp <- cor(x=sub_expression_matrix, y=NULL, method="spearman")
        distance_values_temp <- as.numeric(distance_values_temp[upper.tri(distance_values_temp)])

        cell_types_groups_temp <- as.character(rep(paste(cell_types_groups.unique[i,1],
                                                         cell_types_groups.unique[i,2],
                                                         sep='_'),
                                                times = length(distance_values_temp)))
        cellTypes_temp <- as.character(rep(cell_types_groups.unique[i,1],
                                           times = length(distance_values_temp)))
        group_temp <- as.character(rep(cell_types_groups.unique[i,2],
                                       times = length(distance_values_temp)))

        distance_values <- c(distance_values, distance_values_temp)
        cell_types_groups.vector <- c(cell_types_groups.vector, cell_types_groups_temp)
        cellTypes <- c(cellTypes, cellTypes_temp)
        group <- c(group, group_temp)
      }

      distance_stat <- data.frame(x=cell_types_groups.vector,
                                  y=distance_values,
                                  stringsAsFactors=FALSE)
      colnames(distance_stat) <- c('orig.ident', 'distance')
      distance_stat$group <- group
      distance_stat$cellTypes <- cellTypes
      distance_stat$cellTypes <- factor(distance_stat$cellTypes, levels=cellTypeOrders)
      # get the longest datasetID
      datasetID_lengths <- sapply(unique(group), nchar)
      max_length <- max(datasetID_lengths)
      cellTypes_length <- length(unique(cellTypes))
      
      pdf(file = paste0(output.dir,'/Heterogeneity_spearmanCorrelation.pdf'),
          width = (8 + max_length*0.1+cellTypes_length*0.5),
          height = 8)
          print(ggplot(data=distance_stat, aes(x=cellTypes,
                                               y=distance,
                                               fill=group),
                       color='#ffffff')+
                stat_boxplot(geom ='errorbar', width = 0.6,position =position_dodge(0.8))+
                geom_boxplot(outlier.shape = NA, position=position_dodge(0.8))+
                #scale_fill_manual(values = c(alpha('#e41a1c',1), alpha('#377eb8',1)))+
                xlab(label='cell types')+ylab(label='Spearman correlation')+
                #ggpubr::stat_compare_means(aes(group=group), method = 'wilcox.test', label = "p.signif")+
                ylim((min(distance_stat$distance)-0.05),(max(distance_stat$distance)+0.05))+mytheme+
                  theme(axis.text.x = element_text(size = 20,
                                                   family = "sans",
                                                   color = "black",
                                                   face = "bold",
                                                   vjust = 0,
                                                   hjust = 0.5))
                  )
      dev.off()

      png(file = paste0(output.dir,'/Heterogeneity_spearmanCorrelation.png'),
          width = (800 + max_length*10+cellTypes_length*50),
          height = 800)
          print(ggplot(data=distance_stat, aes(x=cellTypes,
                                               y=distance,
                                               fill=group),
                       color='#ffffff')+
                stat_boxplot(geom ='errorbar', width = 0.6,position =position_dodge(0.8))+
                geom_boxplot(outlier.shape = NA, position=position_dodge(0.8))+
                # scale_fill_manual(values = c(alpha('#e41a1c',1), alpha('#377eb8',1)))+
                xlab(label='cell types')+ylab(label='Spearman correlation')+
                # ggpubr::stat_compare_means(aes(group=group), method = 'wilcox.test', label = "p.signif")+
                ylim((min(distance_stat$distance)-0.05),(max(distance_stat$distance)+0.05))+mytheme+
                  theme(axis.text.x = element_text(size = 20,
                                                   family = "sans",
                                                   color = "black",
                                                   face = "bold",
                                                   vjust = 0,
                                                   hjust = 0.5))
              )
      dev.off()  
  }
}

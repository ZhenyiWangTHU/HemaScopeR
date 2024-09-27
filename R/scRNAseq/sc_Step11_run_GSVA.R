#' Run GSVA Analysis
#'
#' Function \code{run_GSVA} performs GSVA (Gene Set Variation Analysis) analysis on a given single-cell RNA-seq dataset.
#' 
#' @param sc_object A Seurat object containing the single-cell RNA-seq data.
#' @param GSVA.genelist A list of gene sets for GSVA analysis.
#' @param GSVA.cellTypes A character vector specifying the cell types or labels for each cell.
#' @param GSVA.cellTypes.orders A character vector specifying the order of cell types for visualization.
#' @param GSVA.cellGroups A character vector specifying the cell groups or conditions for each cell.
#' @param GSVA.identify.cellType.features Logical. If TRUE, identify cell type-specific features.
#' @param GSVA.identify.diff.features Logical. If TRUE, identify differentially expressed features between cell groups.
#' @param GSVA.comparison.design A list specifying the experimental design for differential GSVA analysis.
#' @param OrgDB An organism-specific annotation database (OrgDb) for gene symbol conversion. e.g. org.Mm.eg.db or org.Hs.eg.db.
#' @param output.dir The directory where the GSVA results and output files will be saved.
#' 
#' @details
#' This function runs GSVA analysis, which calculates enrichment scores for gene sets in each cell using the provided gene list.
#' It also performs differential GSVA analysis between specified cell groups and generates heatmaps of the results. GSVA was referred to
#' https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf;http://revigo.irb.hr/;https://www.slideshare.net/dgrapov/proteomics-workshop-2014-lab-dmitry-grapov
#' 
#' @return None. This function saves GSVA results and output files in the specified output directory.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#' 
#' @export

run_GSVA = function(sc_object=NULL,
                    GSVA.genelist=NULL,
                    GSVA.cellTypes=NULL,
                    GSVA.cellTypes.orders=NULL,
                    GSVA.cellGroups=NULL,
                    GSVA.identify.cellType.features=TRUE,
                    GSVA.identify.diff.features=FALSE,
                    GSVA.comparison.design=NULL,
                    OrgDB=NULL,
                    output.dir=NULL){
    
  if(GSVA.identify.cellType.features){
      dir.create(paste0(output.dir, '/GSVA.identify.cellType.features/'))
      filtered.data <- sc_object
      filtered.data@meta.data$cellTypes <- GSVA.cellTypes
      Idents(filtered.data) <- filtered.data@meta.data$selectLabels
      filtered.data@meta.data$pseudoGroup <- filtered.data@meta.data$selectLabels

      # table count
      cellTypeCounts <- table(filtered.data@meta.data$selectLabels)
      # >3
      cellTypesOverThree <- names(cellTypeCounts[cellTypeCounts > 3])
      filtered.data <- subset(filtered.data, subset = selectLabels %in% cellTypesOverThree)
      
      for(i in unique(Idents(filtered.data))){
        temp_index <- which(filtered.data@meta.data$pseudoGroup==i)
        C1_index <- base::sample(temp_index, size = floor(length(temp_index)/3), replace = FALSE)
        filtered.data@meta.data$pseudoGroup[C1_index] <- 'Cluster1'
    
        C2_C3 <- setdiff(temp_index, C1_index)
        C2_index <- base::sample(C2_C3, size = floor(length(C2_C3)/2), replace = FALSE)
        filtered.data@meta.data$pseudoGroup[C2_index] <- 'Cluster2'
    
        C3_index <- setdiff(C2_C3, C2_index)
        filtered.data@meta.data$pseudoGroup[C3_index] <- 'Cluster3'
      }
    
      averaged.filtered.data <- AverageExpression(object=filtered.data,
                                                  return.seurat = FALSE,
                                                  slot = "data",
                                                  add.ident = "pseudoGroup")
    
      averaged.filtered.data <- averaged.filtered.data[[1]]
    
      entrezID_symbol <- AnnotationDbi::select(OrgDB,
                                               keys = rownames(averaged.filtered.data),
                                               columns = c("ENTREZID","SYMBOL"),
                                               keytype = "SYMBOL")
    
      entrezID_symbol <- entrezID_symbol[-which(is.na(entrezID_symbol$ENTREZID)==TRUE),]
      averaged.filtered.data <- subset(averaged.filtered.data, rownames(averaged.filtered.data)%in%entrezID_symbol$SYMBOL)
      rownames(averaged.filtered.data) <- mapvalues(rownames(averaged.filtered.data),
                                                    from = entrezID_symbol$SYMBOL,
                                                    to = entrezID_symbol$ENTREZID)

      # GSVA KEGG
      gsva.result <- gsva(expr=as.matrix(averaged.filtered.data),
                          gset.idx.list=GSVA.genelist,
                          kcdf="Gaussian",
                          parallel.sz=1,
                          min.sz=2)
      saveRDS(gsva.result, file=paste0(gsub("/Step11.GSVA/$", "/RDSfiles", output.dir), '/gsva.result.rds'))
      cellTypes_vec <- unique(filtered.data@meta.data$selectLabels)
      gene_set_scores <- gsva.result
      write.csv(gene_set_scores, file = paste0(output.dir,'/GSVA.identify.cellType.features/All_GSVA_result.csv'))
      gsva.result_vec <- c()
      
      for(i in cellTypes_vec){
        grouP <- c(rep("A", 3), rep("B", (length(unique(filtered.data@meta.data$selectLabels))-1)*3)) %>% as.factor()
        desigN <- model.matrix(~ grouP + 0)
          
        rownames(desigN) <- c(paste(rep(i,3), c('Cluster1','Cluster2','Cluster3'), sep='_'),
                              setdiff(colnames(gene_set_scores), paste(rep(i,3), c('Cluster1','Cluster2','Cluster3'), sep='_')))
          
        comparE <- makeContrasts(grouPB - grouPA, levels=desigN)
        # order the score matrix
        gsva_score <- gene_set_scores[, which(colnames(gene_set_scores) %in% rownames(desigN))]
        colnames_gsva_score <- factor(colnames(gsva_score), levels = rownames(desigN))
        gsva_score <- gsva_score[, order(colnames_gsva_score, decreasing = FALSE)]
          
        if(ncol(gsva_score) < 6){next}
        # limma
        fiT <- lmFit(gsva_score, desigN)
        fiT2 <- contrasts.fit(fiT, comparE)
        fiT3 <- eBayes(fiT2)
        keggDiff <- topTable(fiT3, coef=1, number=nrow(gsva_score))
        keggDiff <- keggDiff[order(keggDiff$logFC, decreasing =TRUE),]
        keggDiff <- keggDiff[which(keggDiff$adj.P.Val < 0.05),]
        gsva.result_vec <- unique(c(gsva.result_vec, rownames(keggDiff)[1:5]))

        if(nrow(keggDiff)>0){
            write.csv(keggDiff,file = paste0(output.dir,
                                             '/GSVA.identify.cellType.features/GSVA_', 
                                             as.character(make.names(i)),
                                             '_result.csv'))
        }
      }
      
        select_pathways <- toupper(gsva.result_vec)
        select_pathways <- select_pathways[!is.na(select_pathways)]

        if(length(unique(select_pathways) >= 2)){
              gsva.result <- as.data.frame(subset(gsva.result, rownames(gsva.result)%in%unique(select_pathways)))
              gsva.result$pathways <- factor(rownames(gsva.result), levels = unique(select_pathways))
              gsva.result <- gsva.result[order(gsva.result$pathways, decreasing = FALSE),]
              gsva.result <- gsva.result[, -which(colnames(gsva.result)=='pathways')]
              if(is.null(GSVA.cellTypes.orders)){
                  GSVA.cellTypes.orders <- unique(Idents(filtered.data))
                  GSVA.cellTypes.orders <- lapply(GSVA.cellTypes.orders, function(element) {
                                                  paste(rep(element, each = 3), c('Cluster1', 'Cluster2', 'Cluster3'), sep = '_')
                                                }) %>% unlist()
              }else{
                  GSVA.cellTypes.orders <- lapply(GSVA.cellTypes.orders, function(element) {
                                                  paste(rep(element, each = 3), c('Cluster1', 'Cluster2', 'Cluster3'), sep = '_')
                                                }) %>% unlist()              
              }
              colnames_gsva.result <- factor(colnames(gsva.result), levels = GSVA.cellTypes.orders)
              gsva.result <- gsva.result[,order(colnames_gsva.result, decreasing = FALSE)]
              rownames(gsva.result) <- tolower(rownames(gsva.result))
              ordered_genes_expression_matrix <- gsva.result   
              pheatmap::pheatmap(as.matrix(ordered_genes_expression_matrix),
                                 cluster_rows = TRUE,
                                 cluster_cols =F,
                                 scale = "row" ,
                                 legend_breaks= ceiling(seq(min(ordered_genes_expression_matrix),
                                                            max(ordered_genes_expression_matrix),0.01)),
                                 color = colorRampPalette(colors = c("#377eb8","#deebf7","#e41a1c"))(length(seq(min(ordered_genes_expression_matrix),max(ordered_genes_expression_matrix),0.01))),
                                 border_color = NA,
                                 breaks= seq(min(ordered_genes_expression_matrix),
                                             max(ordered_genes_expression_matrix),
                                             by=0.01),
                                 show_colnames = T, show_rownames = T,
                                 #annotation_col  = annotation_col,
                                 #annotation_colors = ann_colors_C,
                                 fontsize =3,
                                 filename=paste0(output.dir,'/GSVA.identify.cellType.features/GSVA_Heatmap.pdf'),
                                 #height=15,
                                 #weight=15
                                 cellwidth=3,
                                 cellheight=3
                                )
        
              pheatmap::pheatmap(as.matrix(ordered_genes_expression_matrix),
                                 cluster_rows = TRUE,
                                 cluster_cols =F,
                                 scale = "row" ,
                                 legend_breaks= ceiling(seq(min(ordered_genes_expression_matrix),
                                                            max(ordered_genes_expression_matrix),0.01)),
                                 color = colorRampPalette(colors = c("#377eb8","#deebf7","#e41a1c"))(length(seq(min(ordered_genes_expression_matrix),max(ordered_genes_expression_matrix),0.01))),
                                 border_color = NA,
                                 breaks= seq(min(ordered_genes_expression_matrix),
                                             max(ordered_genes_expression_matrix),
                                             by=0.01),
                                 show_colnames = T, show_rownames = T,
                                 #annotation_col  = annotation_col,
                                 #annotation_colors = ann_colors_C,
                                 fontsize =3,
                                 filename=paste0(output.dir,'/GSVA.identify.cellType.features/GSVA_Heatmap.png'),
                                 #height=15,
                                 #weight=15
                                 cellwidth=3,
                                 cellheight=3
                                )  
          }
  }

  if(GSVA.identify.diff.features&(!is.null(GSVA.comparison.design))){
      dir.create(paste0(output.dir, '/GSVA.identify.diff.features/'))
      filtered.data <- sc_object
      filtered.data@meta.data$selectLabels <- GSVA.cellTypes
      filtered.data@meta.data$group <- GSVA.cellGroups
      filtered.data@meta.data$cellTypes_group <- paste(filtered.data@meta.data$selectLabels,
                                                       filtered.data@meta.data$group,
                                                       sep='_')
      Idents(filtered.data) <- filtered.data@meta.data$cellTypes_group
      filtered.data@meta.data$pseudoGroup <- filtered.data@meta.data$cellTypes_group
      # table count
      cellTypeCounts <- table(filtered.data@meta.data$selectLabels)
      # >3
      cellTypesOverThree <- names(cellTypeCounts[cellTypeCounts > 3])
      filtered.data <- subset(filtered.data, subset = selectLabels %in% cellTypesOverThree)
    
      for(i in unique(Idents(filtered.data))){
        temp_index <- which(filtered.data@meta.data$pseudoGroup==i)
        C1_index <- base::sample(temp_index, size = floor(length(temp_index)/3), replace = FALSE)
        filtered.data@meta.data$pseudoGroup[C1_index] <- 'Cluster1'
    
        C2_C3 <- setdiff(temp_index, C1_index)
        C2_index <- base::sample(C2_C3, size = floor(length(C2_C3)/2), replace = FALSE)
        filtered.data@meta.data$pseudoGroup[C2_index] <- 'Cluster2'
    
        C3_index <- setdiff(C2_C3, C2_index)
        filtered.data@meta.data$pseudoGroup[C3_index] <- 'Cluster3'
      }
    
      averaged.filtered.data <- AverageExpression(object=filtered.data,
                                                  return.seurat = FALSE,
                                                  slot = "data",
                                                  add.ident = "pseudoGroup")
    
      averaged.filtered.data <- averaged.filtered.data[[1]]
    
      entrezID_symbol <- AnnotationDbi::select(OrgDB,
                                               keys = rownames(averaged.filtered.data),
                                               columns = c("ENTREZID","SYMBOL"),
                                               keytype = "SYMBOL")
    
      entrezID_symbol <- entrezID_symbol[-which(is.na(entrezID_symbol$ENTREZID)==TRUE),]
      averaged.filtered.data <- subset(averaged.filtered.data, rownames(averaged.filtered.data)%in%entrezID_symbol$SYMBOL)
      rownames(averaged.filtered.data) <- mapvalues(rownames(averaged.filtered.data),
                                                    from = entrezID_symbol$SYMBOL,
                                                    to = entrezID_symbol$ENTREZID)

      # GSVA KEGG
      gsva.result <- gsva(expr=as.matrix(averaged.filtered.data),
                          gset.idx.list=GSVA.genelist,
                          kcdf="Gaussian",
                          parallel.sz=1,
                          min.sz=2)
      saveRDS(gsva.result, file=paste0(gsub("/Step11.GSVA/$", "/RDSfiles", output.dir), '/gsva.result.rds'))
      cellTypes_vec <- unique(filtered.data@meta.data$selectLabels)
      gene_set_scores <- gsva.result
      write.csv(gene_set_scores, file = paste0(output.dir,'/GSVA.identify.diff.features/All_GSVA_result.csv'))

      gsva.result_vec <- c()

      for(j in 1:length(GSVA.comparison.design)){
          GSVA.comparison.design.temp <- GSVA.comparison.design[[j]]
          GSVA.comparison.design.temp.1 <- unlist(GSVA.comparison.design.temp[1])
          GSVA.comparison.design.temp.2 <- unlist(GSVA.comparison.design.temp[2])
          for(i in cellTypes_vec){
          grouP <- c(rep("A", 3*length(GSVA.comparison.design.temp.1)), 
                     rep("B", 3*length(GSVA.comparison.design.temp.2))) %>% as.factor()
          desigN <- model.matrix(~ grouP + 0)
          rownames(desigN) <- paste(rep(i, (3*length(GSVA.comparison.design.temp.1) + 3*length(GSVA.comparison.design.temp.2))),
                                    c(rep(GSVA.comparison.design.temp.1, each=3), rep(GSVA.comparison.design.temp.2, each=3)),
                                    c(rep(c('Cluster1','Cluster2','Cluster3'), length(GSVA.comparison.design.temp.1)),   
                                    rep(c('Cluster1','Cluster2','Cluster3'), length(GSVA.comparison.design.temp.2))),
                                    sep='_')
              
          desigN <- desigN[which(rownames(desigN)%in% intersect(rownames(desigN), colnames(gene_set_scores))), ]
          if(length(unique(desigN[,1]))==1){next}
          # order the score matrix
          gsva_score <- gene_set_scores[, which(colnames(gene_set_scores)%in%intersect(rownames(desigN), colnames(gene_set_scores)))]
          if(is.null(dim(gsva_score))){next}
              
          comparE <- makeContrasts(grouPB - grouPA, levels = desigN)
        
          colnames(gsva_score) <- factor(colnames(gsva_score), levels = rownames(desigN))
          gsva_score <- gsva_score[,order(colnames(gsva_score), decreasing = TRUE)]
          if(ncol(gsva_score) < 6){next}
          # limma
 
          fiT <- lmFit(gsva_score, desigN)
          fiT2 <- contrasts.fit(fiT, comparE)
          fiT3 <- eBayes(fiT2)
          keggDiff <- topTable(fiT3, coef=1, number=nrow(gsva_score))
          keggDiff <- keggDiff[order(keggDiff$logFC, decreasing =TRUE),]
          keggDiff <- keggDiff[which(keggDiff$adj.P.Val < 0.05),]
          gsva.result_vec <- unique(c(gsva.result_vec, rownames(keggDiff)[1:10]))
          if(nrow(keggDiff)>0){    
              write.csv(keggDiff,file = paste0(output.dir,
                                             '/GSVA.identify.diff.features/GSVA_', as.character(make.names(i)), '_',
                                             as.character(paste(GSVA.comparison.design.temp.2, collapse = ' And ')),
                                             ' vs ',
                                             as.character(paste(GSVA.comparison.design.temp.1, collapse = ' And ')),
                                             '_result.csv'))
            }
          }
          select_pathways <- toupper(gsva.result_vec)
          select_pathways <- select_pathways[!is.na(select_pathways)]
          if(length(unique(select_pathways) >= 2)){
              gsva.result <- as.data.frame(subset(gsva.result, rownames(gsva.result)%in%unique(select_pathways)))
              gsva.result$pathways <- factor(rownames(gsva.result), levels = unique(select_pathways))
              gsva.result <- gsva.result[order(gsva.result$pathways, decreasing = FALSE),]
              gsva.result <- gsva.result[, -which(colnames(gsva.result)=='pathways')]
              if(is.null(GSVA.cellTypes.orders)){
                  GSVA.cellTypes.orders <- unique(Idents(filtered.data))
                  GSVA.cellTypes.orders <- lapply(GSVA.cellTypes.orders, function(element) {
                                                  paste(rep(element, each = 3), c('Cluster1', 'Cluster2', 'Cluster3'), sep = '_')
                                                }) %>% unlist()
              }else{
                  GSVA.cellTypes.orders <- lapply(GSVA.cellTypes.orders, function(element) {
                                                  paste(rep(element, each = 3), c('Cluster1', 'Cluster2', 'Cluster3'), sep = '_')
                                                }) %>% unlist()              
              }
              colnames_gsva.result <- factor(colnames(gsva.result), levels = GSVA.cellTypes.orders)
              gsva.result <- gsva.result[,order(colnames_gsva.result, decreasing = FALSE)]
              rownames(gsva.result) <- tolower(rownames(gsva.result))
              ordered_genes_expression_matrix <- gsva.result
        
              pheatmap::pheatmap(as.matrix(ordered_genes_expression_matrix),
                                 cluster_rows = TRUE,
                                 cluster_cols =F,
                                 scale = "row" ,
                                 legend_breaks= ceiling(seq(min(ordered_genes_expression_matrix),
                                                            max(ordered_genes_expression_matrix),0.01)),
                                 color = colorRampPalette(colors = c("#377eb8","#deebf7","#e41a1c"))(length(seq(min(ordered_genes_expression_matrix),max(ordered_genes_expression_matrix),0.01))),
                                 breaks= seq(min(ordered_genes_expression_matrix),
                                             max(ordered_genes_expression_matrix),
                                             by=0.01),
                                 show_colnames = T, show_rownames = T,
                                 #annotation_col  = annotation_col,
                                 #annotation_colors = ann_colors_C,
                                 fontsize =10,
                                 filename=paste0(output.dir,'/GSVA.identify.diff.features/',as.character(paste(GSVA.comparison.design.temp.2, collapse = ' And ')),
                                                 ' vs ',as.character(paste(GSVA.comparison.design.temp.1, collapse = ' And ')),' GSVA_Heatmap.pdf'),
                                 #height=15,
                                 #weight=15,
                                 cellwidth=20,
                                 cellheight=20)
        
              pheatmap::pheatmap(as.matrix(ordered_genes_expression_matrix),
                                 cluster_rows = TRUE,
                                 cluster_cols =F,
                                 scale = "row" ,
                                 legend_breaks= ceiling(seq(min(ordered_genes_expression_matrix),
                                                            max(ordered_genes_expression_matrix),0.01)),
                                 color = colorRampPalette(colors = c("#377eb8","#deebf7","#e41a1c"))(length(seq(min(ordered_genes_expression_matrix),max(ordered_genes_expression_matrix),0.01))),
                                 breaks= seq(min(ordered_genes_expression_matrix),
                                             max(ordered_genes_expression_matrix),
                                             by=0.01),
                                 show_colnames = T, show_rownames = T,
                                 #annotation_col  = annotation_col,
                                 #annotation_colors = ann_colors_C,
                                 fontsize =10,
                                 filename=paste0(output.dir,'/GSVA.identify.diff.features/',as.character(paste(GSVA.comparison.design.temp.2, collapse = ' And ')),
                                                 ' vs ',as.character(paste(GSVA.comparison.design.temp.1, collapse = ' And ')),' GSVA_Heatmap.png'),
                                 #height=15,
                                 #weight=15,
                                 cellwidth=20,
                                 cellheight=20)  
          }
        }  
     }
}

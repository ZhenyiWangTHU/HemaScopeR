#' Run SCENIC Analysis
#'
#' Function \code{run_SCENIC} performs SCENIC (Single-Cell Regulatory Network Inference and Clustering) analysis on a given single-cell RNA-seq dataset.
#' 
#' @param countMatrix A matrix containing the raw counts of the single-cell RNA-seq data.
#' @param cellTypes A character vector specifying the cell types or labels for each cell.
#' @param datasetID A character vector specifying the dataset IDs for each cell.
#' @param cellTypes_colors A named vector of colors for cell type visualization.
#' @param cellTypes_orders A character vector specifying the desired order of cell types.
#' @param groups_colors A named vector of colors for grouping visualization.
#' @param groups_orders A character vector specifying the desired order of groups.
#' @param Org A character vector specifying the organism ('mmu' for mouse or 'hsa' for human).
#' @param output.dir The directory where the SCENIC results and output files will be saved.
#' @param pythonPath The path to the Python environment to use for the analysis.
#' 
#' @details
#' This function runs SCENIC analysis, including the construction of a co-expression network, gene filtering, correlation, and the GENIE3 algorithm to infer regulatory networks.
#' It also saves various output files for downstream analysis and visualization.
#' 
#' @return None. This function saves SCENIC results and output files in the specified output directory.
#' 
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#' 
#' @export

run_SCENIC = function(countMatrix = NULL,
                      cellTypes = NULL,
                      datasetID = NULL,
                      cellTypes_colors = NULL,
                      cellTypes_orders = NULL,
                      groups_colors = NULL,
                      groups_orders = NULL,
                      Org = NULL,
                      output.dir = NULL,
                      pythonPath = NULL,
                      databasePath = NULL){
  if(Org == 'mmu'){
    org="mgi"
    dbDir=databasePath
    dbs=c("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
          "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
    data("motifAnnotations_mgi_v9", package="RcisTarget")
    #rename the motif annnotion by attributing it to the variable that is in the error
    motifAnnotations_mgi <<- motifAnnotations_mgi_v9
  }else if(Org == 'hsa'){
    org="hgnc"
    dbDir=databasePath
    dbs=c("hg19-tss-centered-10kb-7species.mc9nr.feather",
          "hg19-500bp-upstream-7species.mc9nr.feather")
    data("motifAnnotations_hgnc_v9", package="RcisTarget")
    #rename the motif annnotion by attributing it to the variable that is in the error
    motifAnnotations_hgnc <<- motifAnnotations_hgnc_v9
  }else{
    stop("Org should be 'mmu' or 'hsa'.")
  }
    
  setwd(output.dir)
  # input count matrix
  exprMat <- as.matrix(countMatrix)
  # input cell info
  cellInfo <- data.frame(cellTypes)
  colnames(cellInfo) <- 'CellTypes'
  dir.create(paste0(output.dir, "/int"))
  print(paste0(output.dir, "/int"))
  
  saveRDS(cellInfo, file=paste0(output.dir, "/int/cellInfo.Rds"))
  # Initialize SCENIC settings
  myDatasetTitle = "scRNASeq"
  scenicOptions <- initializeScenic(org=org,
                                    dbDir=dbDir,
                                    dbs=dbs,
                                    datasetTitle=myDatasetTitle)

  scenicOptions@inputDatasetInfo$cellInfo <- paste0(output.dir, "/int/cellInfo.Rds")
  # scenicOptions@inputDatasetInfo$colVars <- paste0(output.dir, "/int/colVars.Rds")

  scenicOptions@settings$verbose <- TRUE
  # scenicOptions@settings$nCores <- parallel::detectCores(logical=F)
  scenicOptions@settings$nCores <- 1
  scenicOptions@settings$seed <- 2023

  # Save to use at a later time...
  saveRDS(scenicOptions, file=paste0(output.dir, "/int/scenicOptions.Rds"))

  # Co-expression network
  # Gene filter/selection
  # (Adjust minimum values according to your dataset)
  genesKept <- geneFiltering(exprMat,
                             scenicOptions=scenicOptions,
                             minCountsPerGene=0,
                             minSamples=ncol(exprMat)*.005)

  exprMat_filtered <- exprMat[genesKept, ]

  # Correlation
  runCorrelation(exprMat_filtered, scenicOptions)
    
  # GRNBoost
  # export mat for GRNBoost
  exportsForArboreto(exprMat=exprMat_filtered,
                     scenicOptions=scenicOptions,
                     dir = paste0(output.dir, "/int/"))
  
  # run pyscenic_GRNboost
  if(is.null(pythonPath)==FALSE){ reticulate::use_python(pythonPath) }else{print('Please set the path of Python.')}  
 
  reticulate::py_run_string(paste0("import os\nos.chdir('", output.dir, "')"))   
  reticulate::py_run_file(file.path(system.file(package = "HemaScopeR"), "python/run_pyscenic_grnboost.py"), convert = FALSE)
   
  # import GRNBoost results
  GRNBoost_output <- read.delim("./03.grnboot2_network_links.tsv", header=TRUE)
  colnames(GRNBoost_output) <- c("TF", "Target", "weight")
  GRNBoost_output$weight <- as.numeric(GRNBoost_output$weight)
  saveRDS(GRNBoost_output, file=paste0(output.dir, "/int/1.4_GENIE3_linkList.Rds"))

  # # GENIE3
  # # Optional: add log (if it is not logged/normalized already)
  exprMat_filtered <- log2(exprMat_filtered + 1)
  # # Run GENIE3
  # runGenie3(exprMat_filtered, scenicOptions)

  # Build and score the GRN (runSCENIC_…)
  runSCENIC_1_coexNetwork2modules(scenicOptions)
  runSCENIC_2_createRegulons(scenicOptions)
  #scenicOptions@settings$nCores <- 1
  #saveRDS(scenicOptions, file=paste0(output.dir, "/int/scenicOptions.Rds"))  
  runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)  
  scenicOptions <- readRDS(file=paste0(output.dir, "/int/scenicOptions.Rds"))
  scenicOptions@settings$nCores <- 10  
  runSCENIC_4_aucell_binarize(scenicOptions, skipHeatmaps = TRUE)

  # Exploring/interpreting the results
  loadInt(scenicOptions)
  # auc in each single cell
  geneset_auc <- readRDS(paste0(output.dir, "/int/3.4_regulonAUC.Rds"))
  #geneset_auc <- as.data.frame(slot(geneset_auc@assays, name='data')[[1]])
  geneset_auc <- as.data.frame(geneset_auc@assays@data)
  geneset_auc <- t(geneset_auc)
  geneset_auc <- geneset_auc[c(-1,-2), ]
  write.table(geneset_auc,
              paste0(output.dir, "/int/05.geneset_auc_per_cell.txt"),
              quote=F,
              sep="\t")

  # regulon activity
  regulon_activity <- readRDS(paste0(output.dir, "/int/4.1_binaryRegulonActivity.Rds"))
  regulon_activity <- as.data.frame(regulon_activity)
  regulon_activity <- t(regulon_activity)
  write.table(regulon_activity,
              paste0(output.dir, "/int/06.regulon_activity_per_cell.txt"),
              quote=F,
              sep="\t")

  # auc threshold
  auc_threshold <- read.table(file = paste0(output.dir, "/int/3.5_AUCellThresholds_Info.tsv"),
                              sep = '\t',
                              header = TRUE)

  auc_threshold <- subset(auc_threshold, select=c(regulon, threshold, nCellsAssigned, nGenes))
  write.table(auc_threshold,
              paste0(output.dir, "/int/07.auc_threshold_per_regulon.txt"),
              quote=F,
              sep="\t")

  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

  # heatmap of TFs in each cell
  #regulonAUC <- as.data.frame(slot(regulonAUC@assays, name='data')[[1]])
  regulonAUC <- as.data.frame(regulonAUC@assays@data)
  regulonAUC <- regulonAUC[, -c(1,2)]
  ordergenes_expression_matrix <- regulonAUC
  if(is.null(cellTypes_orders)){cellTypes_orders <- unique(cellTypes)}
  if(is.null(groups_orders)){groups_orders <- unique(datasetID)}  
  annotation_col_C_andcluster = data.frame(cellTypes = factor(cellTypes, levels=cellTypes_orders),
                                           group = factor(datasetID, levels=groups_orders))
  rownames(annotation_col_C_andcluster) = colnames(ordergenes_expression_matrix)

  if((is.null(cellTypes_colors) == FALSE)&(is.null(cellTypes_orders) == FALSE)&(is.null(groups_colors) == FALSE)&(is.null(groups_orders) == FALSE)){
    if((length(cellTypes_colors)==length(cellTypes_orders))&(length(groups_colors)==length(groups_orders))){
      names(cellTypes_colors) <- levels(annotation_col_C_andcluster$cellTypes)
      names(groups_colors) <- levels(annotation_col_C_andcluster$group)  
      ann_colors_C = list(
        cellTypes = cellTypes_colors,
        group = groups_colors)
    }else{ann_colors_C = NULL}
  }else{ann_colors_C = NULL}

  if(is.null(cellTypes_orders)){cellTypes_orders <- unique(cellTypes)}
  ordergenes_expression_matrix <- t(ordergenes_expression_matrix)
  ordered_genes_expression_matrix <- subset(ordergenes_expression_matrix, cellTypes%in%cellTypes_orders[1])
  for (i in 2:length(cellTypes_orders)) {
    subcluster <- subset(ordergenes_expression_matrix, cellTypes%in%cellTypes_orders[i])
    ordered_genes_expression_matrix <- rbind(ordered_genes_expression_matrix, subcluster)
  }
  ordered_genes_expression_matrix <- t(ordered_genes_expression_matrix)
  annotation_col_C_andcluster$cellTypes <- as.character(annotation_col_C_andcluster$cellTypes)
  annotation_col_C_andcluster$cellTypes <- factor(annotation_col_C_andcluster$cellTypes)
  saveRDS(ordered_genes_expression_matrix,
          file = paste0(output.dir, '/TF.matrix.Rds'))
  pheatmap(as.matrix(ordered_genes_expression_matrix), cluster_rows = T, cluster_cols =F,
           scale = "row" ,
           legend_breaks= ceiling(seq(min(ordered_genes_expression_matrix),
                                      max(ordered_genes_expression_matrix),0.01)),
           color = colorRampPalette(colors = c("#4575b4","#e0f3f8","#d73027"))(length(seq(min(ordered_genes_expression_matrix),
                                                                                          max(ordered_genes_expression_matrix),0.01))),
           breaks= seq(min(ordered_genes_expression_matrix),
                       max(ordered_genes_expression_matrix),
                       by=0.01),
           show_colnames = F,
           show_rownames = T,
           annotation_col  = annotation_col_C_andcluster,
           annotation_colors = ann_colors_C,
           fontsize =1,
           treeheight_row=10,
           filename=paste0(output.dir, "/int/TFs_Heatmap.pdf"),
           height=7,
           weight=15
  )

  pheatmap(as.matrix(ordered_genes_expression_matrix), cluster_rows = T, cluster_cols =F,
           scale = "row" ,
           legend_breaks= ceiling(seq(min(ordered_genes_expression_matrix),
                                      max(ordered_genes_expression_matrix),0.01)),
           color = colorRampPalette(colors = c("#4575b4","#e0f3f8","#d73027"))(length(seq(min(ordered_genes_expression_matrix),
                                                                                          max(ordered_genes_expression_matrix),0.01))),
           breaks= seq(min(ordered_genes_expression_matrix),
                       max(ordered_genes_expression_matrix),
                       by=0.01),
           show_colnames = F,
           show_rownames = T,
           annotation_col  = annotation_col_C_andcluster,
           annotation_colors = ann_colors_C,
           fontsize =1,
           treeheight_row=10,
           filename=paste0(output.dir, "/int/TFs_Heatmap.png"),
           height=7,
           weight=15
  )
  # TF info
  #  regulons <- loadInt(scenicOptions, "regulons")
  #  head(regulons)
  #  aucell_regulons <- loadInt(scenicOptions, "aucell_regulons")
  #  head(aucell_regulons)
  regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
  write.csv(regulonTargetsInfo, file=paste0(output.dir, "/int/regulonTargetsInfo.csv"))
}
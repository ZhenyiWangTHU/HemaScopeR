#' run_CellChat - Function for Cell-Cell Communication Analysis with CellChat
#' 
#' This function performs cell-cell communication analysis using the CellChat package. It takes expression data and cluster labels as input, identifies cell-cell communication networks, and visualizes the results, including interaction networks and signaling pathways.
#'
#' @param data.input A matrix of expression data, where rows represent genes and columns represent cells. Row names should be in the format of gene symbols.
#' @param labels A vector of cluster labels for each cell, corresponding to the columns of data.input.
#' @param cell.orders A character vector specifying the order of cell types or clusters in the analysis.
#' @param cell.colors A character vector specifying colors for cell types or clusters in the analysis.
#' @param sample.names A vector of sample or cell names, corresponding to the columns of data.input.
#' @param Org A string indicating the organism used in the analysis. It should be either "mmu" (mouse) or "hsa" (human).
#' @param sorting A logical value indicating whether to consider cell population size in communication analysis.
#' @param output.dir The directory path where the output files, including plots and data, will be saved.
#'
#' @details The function takes expression data, cluster labels, and other parameters to perform cell-cell communication analysis using the CellChat package. It includes the following steps:
#'
#' 1. Data input and preprocessing.
#' 2. Initialization of a CellChat object.
#' 3. Set the ligand-receptor interaction database based on the specified organism.
#' 4. Preprocess the expression data for cell-cell communication analysis.
#' 5. Identify overexpressed genes and interactions.
#' 6. Project data based on protein-protein interaction networks.
#' 7. Inference of cell-cell communication network.
#' 8. Visualization of the communication network.
#' 9. Systems analysis of cell-cell communication network.
#'
#' The results are saved as PDF files containing various network visualizations and analysis plots.
#'
#' @return The function does not return any values directly, but it saves PDF files with various network visualizations, plots, and data files that can be used for further analysis.
#'
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#'
#' @export
run_CellChat = function(data.input=NULL,
                        labels = NULL,
                        cell.orders = NULL,
                        cell.colors = NULL,
                        sample.names = NULL,
                        Org = NULL,
                        sorting = FALSE,
                        output.dir = NULL){
  setwd(output.dir)  
  # Part I: Data input & processing and initialization of CellChat object

  # the row names of data.input should be title format gene symbol
  data.input <- normalizeData(data.raw=data.input)
  if(is.null(cell.orders)){cell.orders <- unique(labels)} 
  identity <- data.frame(labels = factor(labels, levels = cell.orders),
                         row.names = sample.names)

  # Create a CellChat object
  cellchat <- createCellChat(object = data.input,
                             meta = identity,
                             group.by = 'labels')

  # number of cells in each cell group
  groupSize <- as.numeric(table(cellchat@idents))

  # Set the ligand-receptor interaction database
  if(Org =='mmu'){
    CellChatDB <- CellChatDB.mouse
  }else if(Org =='hsa'){
    CellChatDB <- CellChatDB.human
  }else{
    stop('The Org should be "mmu" or "hsa".')
  }

  # use Secreted Signaling for cell-cell communication analysis
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  CellChatDB.use <- CellChatDB
  # set the used database in the object
  cellchat@DB <- CellChatDB.use

  # Preprocessing the expression data for cell-cell communication analysis
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat)

  # do parallel
  # future::plan("multicore", workers = parallel::detectCores(logical=F))
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  if(Org =='mmu'){
    cellchat <- projectData(cellchat, PPI.mouse)
  }else if(Org =='hsa'){
    cellchat <- projectData(cellchat, PPI.human)
  }else{
    stop('The Org should be "mmu" or "hsa".')
  }

  # Part II: Inference of cell-cell communication network
  # Compute the communication probability and infer cellular communication network
  if(sorting==TRUE){population.size <- FALSE}else{population.size <- TRUE}
  cellchat <- computeCommunProb(cellchat, population.size = population.size, raw.use = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  # Extract the inferred cellular communication network as a data frame
  # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors.
  # Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
  df.net <- subsetCommunication(cellchat)

  # Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)

  # Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)

  groupSize <- as.numeric(table(cellchat@idents))
  pdf(paste0(output.dir,'/circle plot of the number of interactions.pdf'), width = 8, height = 8)
      print(netVisual_circle(cellchat@net$count,
                             color.use = cell.colors,
                             vertex.weight = groupSize,
                             weight.scale = T,
                             label.edge= F,
                             title.name = "Number of interactions"))
  dev.off()

  png(paste0(output.dir,'/circle plot of the number of interactions.png'), width = 800, height = 800)
      print(netVisual_circle(cellchat@net$count,
                             color.use = cell.colors,
                             vertex.weight = groupSize,
                             weight.scale = T,
                             label.edge= F,
                             title.name = "Number of interactions"))
  dev.off()  

  pdf(paste0(output.dir,'/circle plot of the interaction weights or strength.pdf'), width = 8, height = 8)
      print(netVisual_circle(cellchat@net$weight,
                             color.use = cell.colors,
                             vertex.weight = groupSize,
                             weight.scale = T,
                             label.edge= F,
                             title.name = "Interaction weights/strength"))
  dev.off()

  png(paste0(output.dir,'/circle plot of the interaction weights or strength.png'), width = 800, height = 800)
      print(netVisual_circle(cellchat@net$weight,
                             color.use = cell.colors,
                             vertex.weight = groupSize,
                             weight.scale = T,
                             label.edge= F,
                             title.name = "Interaction weights/strength"))
  dev.off()  

  pdf(paste0(output.dir,'/count_heatmap.pdf'), width = 8, height = 8)
      print(netVisual_heatmap(cellchat,
                              color.use = cell.colors,
                              measure = "count",
                              color.heatmap = c("#08306b","#821A4D")))
  dev.off()

  png(paste0(output.dir,'/count_heatmap.png'), width = 800, height = 800)
      print(netVisual_heatmap(cellchat,
                              color.use = cell.colors,
                              measure = "count",
                              color.heatmap = c("#08306b","#821A4D")))
  dev.off()  

  pdf(paste0(output.dir,'/weight_heatmap.pdf'), width = 8, height = 8)
      print(netVisual_heatmap(cellchat,
                              color.use = cell.colors,
                              measure = "weight",
                              color.heatmap = c("#08306b","#821A4D")))
  dev.off()

  png(paste0(output.dir,'/weight_heatmap.png'), width = 800, height = 800)
      print(netVisual_heatmap(cellchat,
                              color.use = cell.colors,
                              measure = "weight",
                              color.heatmap = c("#08306b","#821A4D")))
  dev.off()  

  # Automatically save the plots of the all inferred network for quick exploration
  # Access all the signaling pathways showing significant communications
  pathways.show.all <- cellchat@netP$pathways
  # check the order of cell identity to set suitable vertex.receiver
  # levels(cellchat@idents)

  for (i in 1:length(pathways.show.all)) {
    # Visualize communication network associated with both signaling pathway and individual L-R pairs
    #pdf(paste0(output.dir,'/',pathways.show.all[i],'_circle.pdf'), width = 8, height = 8)
        netVisual(cellchat, color.use = cell.colors, signaling = pathways.show.all[i], layout = "circle", out.format = c("png"))
    #dev.off()

    # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
    pdf(paste0(output.dir,'/',pathways.show.all[i],'_conctribution.pdf'), width = 8, height = 8)
        print(netAnalysis_contribution(cellchat, signaling = pathways.show.all[i]))
    dev.off()

    png(paste0(output.dir,'/',pathways.show.all[i],'_conctribution.png'), width = 800, height = 800)
        print(netAnalysis_contribution(cellchat, signaling = pathways.show.all[i]))
    dev.off()  

    # pdf(paste0(output.dir,'/',as.character(pathways.show.all[i]),'_netVisual_aggregate_circle.pdf'), width = 8, height = 8)
    #     print(netVisual_aggregate(cellchat,
    #                               signaling = pathways.show.all[i],
    #                               layout = "circle",
    #                               edge.width.max = 10,
    #                               signaling.name = pathways.show.all[i]))
    # dev.off()

    pdf(paste0(output.dir,'/',as.character(pathways.show.all[i]),'_LRexpression.pdf'), width = 8, height = 8)
        print(plotGeneExpression(cellchat, color.use = cell.colors, signaling = pathways.show.all[i]))
    dev.off()

    png(paste0(output.dir,'/',as.character(pathways.show.all[i]),'_LRexpression.png'), width = 800, height = 800)
        print(plotGeneExpression(cellchat, color.use = cell.colors, signaling = pathways.show.all[i]))
    dev.off()  
  }

  # Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
  # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  pdf(paste0(output.dir,'/cell-cell communication mediated by multiple ligand-receptors.pdf'), width = 10, height = 10)
    print(netVisual_bubble(cellchat,      
                   sources.use = unique(factor(labels)),
                   targets.use = unique(factor(labels)),
                   remove.isolate = FALSE))
  dev.off()

  png(paste0(output.dir,'/cell-cell communication mediated by multiple ligand-receptors.png'), width = 1000, height = 1000)
    print(netVisual_bubble(cellchat,
                   sources.use = unique(factor(labels)),
                   targets.use = unique(factor(labels)),
                   remove.isolate = FALSE))
  dev.off()  

  # Part IV: Systems analysis of cell-cell communication network
  # Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
  # Compute and visualize the network centrality scores
  # Compute the network centrality scores
  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  # Visualize the dominant senders (sources) and receivers (targets) in a 2D space
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  pdf(paste0(output.dir,'/Signaling role analysis.pdf'), width = 8, height = 8)
    print(netAnalysis_signalingRole_scatter(cellchat, color.use = cell.colors))
  dev.off()

  png(paste0(output.dir,'/Signaling role analysis.png'), width = 800, height = 800)
    print(netAnalysis_signalingRole_scatter(cellchat, color.use = cell.colors))
  dev.off()  

  # Identify signals contributing most to outgoing or incoming signaling of certain cell groups
  pdf(paste0(output.dir,'/outgoing signaling.pdf'), width = 6, height = 10)
    print(netAnalysis_signalingRole_heatmap(cellchat, color.use = cell.colors, pattern = "outgoing",height = 20))
  dev.off()

  png(paste0(output.dir,'/outgoing signaling.png'), width = 600, height = 1000)
    print(netAnalysis_signalingRole_heatmap(cellchat, color.use = cell.colors, pattern = "outgoing",height = 20))
  dev.off()  

  pdf(paste0(output.dir,'/incoming signaling.pdf'), width = 6, height = 10)
    print(netAnalysis_signalingRole_heatmap(cellchat, color.use = cell.colors, pattern = "incoming",height = 20))
  dev.off()

  png(paste0(output.dir,'/incoming signaling.png'), width = 600, height = 1000)
    print(netAnalysis_signalingRole_heatmap(cellchat, color.use = cell.colors, pattern = "incoming",height = 20))
  dev.off()  

  # Part V: Save the CellChat object
  saveRDS(cellchat, file = paste0(gsub("/Step14.Cell_cell_interection/$", "/RDSfiles", output.dir),"/afterCellchat.rds"))
}

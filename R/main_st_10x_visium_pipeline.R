#' #' The pipeline for analyzing 10x Visium data.
#' #'
#' #' Function \code{st_10x_visium_pipeline} encompasses the complete downstream
#' #' analysis pipeline for 10x Visium data within HemascopeR.
#' #' @param input.data.dirs A character of data path, where there are filtered_feature_bc_matrix.h5
#' #' and spatial folder
#' #' @param output.dir A character of path to store the results and figures
#' #' @param sampleName A character of the name of the sample
#' #' @param min.gene An integer representing the minimum number of genes detected in a spot
#' #' @param max.gene An integer representing the maximum number of genes detected in a spot
#' #' @param min.nUMI An integer representing the minimum number of nUMI detected in a spot
#' #' @param max.nUMI An integer representing the maximum number of nUMI detected in a spot
#' #' @param min.spot An integer representing the minimum number of spots expressing each gene
#' #' @param verbose verbose as `Seurat`
#' #' @param normalization.method 'SCTransform' or other choice of `normalization.method` in `NormalizaData`
#' #' @param n.dim.used An integer representing the number of dimension used for finding neighbors
#' #' @param resolution An float parameter of `FindClusters` of `Seurat`
#' #' @param species A character representing the species of sample, `human` or `mouse`
#' #' @param bool.remove.mito A bool value indicating whether removing mitochondrial gene
#' #' @param bool.regress.cycling A bool value indicating whether regressing cycling
#' #' @param bool.regress.cycling.standard A bool value indicating the method using in
#' #' regressing cycling. See more at `https://satijalab.org/seurat/articles/cell_cycle_vignette.html`
#' #' @param s.features A gene vector used to calculate the S phase score. If NULL,
#' #' the default gene set will be used
#' #' @param g2m.features A gene vector used to calculate the G2/M phase score. If NULL,
#' #' the default gene set will be used
#' #' @param genReport A bool value indicating whether generating the report
#' #' @details
#' #' This workflow encompasses data quality control (spot QC and gene QC),
#' #' normalization, PCA dimensionality reduction, clustering and visualization
#' #' @return Return the user with a well-formatted HTML analysis report.
#' #' Additionally, meticulously organize and provide the analyzed results and
#' #' publication-quality vector images for each step of the analysis workflow.
#' #' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#' #' @import Seurat
#' #' @import stringr
#' #' @import ggplot2
#' #' @import patchwork
#' #' @import RColorBrewer
#' #' @import SeuratDisk
#' #' @import dplyr
#' #' @import knitr
#' #' @import kableExtra
#' #'
#' #'
#'
#' st_10x_visium_pipeline_old <- function(
#'         input.data.dir,
#'         output.dir,
#'         sampleName = 'Hema_ST',
#'
#'         # For Step1 Loading
#'         rds.file = FALSE,
#'         filename = "filtered_feature_bc_matrix.h5",
#'         assay = "Spatial",
#'         slice = "slice1",
#'         filter.matrix = TRUE,
#'         to.upper = FALSE,
#'
#'         # For Step2 QC
#'         Step2_QC = TRUE,
#'         min.gene = 200,
#'         min.nUMI = 500,
#'         max.gene = Inf,
#'         max.nUMI = Inf,
#'         min.spot = 0, # for genes
#'         bool.remove.mito = FALSE,
#'         species = 'mouse', # human or mosue
#'
#'         # For Step3 Clustering
#'         Step3_Clustering = TRUE,
#'         normalization.method = 'SCTransform',
#'         npcs = 50,
#'         pcs.used = 1:10,
#'         resolution = 0.8,
#'
#'         # For Step4 Find DEGs
#'         Step4_Find_DEGs = TRUE,
#'         only.pos = TRUE,
#'         min.pct = 0.25,
#'         logfc.threshold = 0.25,
#'         test.use = 'wilcox',
#'
#'         # For Step5 SVF
#'         Step5_SVFs = TRUE,
#'         selection.method = 'moransi',
#'         n.top.show = 10,
#'         n.col.show = 5,
#'
#'         # For Step6 Interaction
#'         Step6_Interaction = TRUE,
#'         commot.signaling_type = 'Secreted Signaling',
#'         commot.database = 'CellChat',
#'         commot.min_cell_pct = 0.05,
#'         commot.dis_thr = 500,
#'         commot.n_permutations = 100,
#'
#'         # For Step7 CNV analysis
#'         Step7_CNV = TRUE,
#'         copykat.genome = NULL,
#'         copykat.LOW.DR = 0.05,
#'         copykat.UP.DR = 0.1,
#'         copykat.win.size = 25,
#'         copykat.distance = "euclidean",
#'         copykat.n.cores = 1,
#'
#'         # For Step8 Deconvolution
#'         Step8_Deconvolution = TRUE,
#'         cell2loc.sc.h5ad.dir = NULL,
#'         cell2loc.sc.max.epoch = 1000,
#'         cell2loc.st.max.epoch = 10000,
#'         cell2loc.use.gpu = TRUE,
#'         cell2loc.use.dataset = 'LymphNode',
#'
#'         # For Step9 Cellcycle
#'         Step9_Cellcycle = TRUE,
#'         s.features = NULL,
#'         g2m.features = NULL,
#'
#'         # For Step10 Nich
#'         Step10_Niche = TRUE,
#'         coexistence.method = 'correlation',
#'         Niche.cluster.n = 4,
#'
#'         # settings
#'         condaenv = 'r-reticulate',
#'         verbose = FALSE,
#'         genReport = TRUE
#' ){
#'     ### Param ###
#'     SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
#'     FeatureColors.bi <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu')))
#'     FeatureColors.one <- colorRampPalette(colors = brewer.pal(n = 9, name = 'YlOrRd'))
#'
#'     if(!dir.exists(file.path(output.dir, sampleName))){
#'         dir.create(file.path(output.dir, sampleName))
#'     }else{
#'         warning(paste0('The new results will overwrite the file under ',
#'                        file.path(output.dir, sampleName)))
#'     }
#'     output.dir <- file.path(output.dir, sampleName)
#'
#'     print(paste0('The results will be saved in ', output.dir))
#'
#'     #### Step1: Loading data ####
#'     print('Loading data...')
#'     st_obj <- st_Loading_Data(
#'         input.data.dir = input.data.dir,
#'         output.dir = file.path(output.dir, 'Step1_Loading_Data'),
#'         sampleName = sampleName,
#'
#'         rds.file = rds.file,
#'         filename = filename,
#'         assay = assay,
#'         slice = slice,
#'         filter.matrix = filter.matrix,
#'         to.upper = to.upper
#'     )
#'
#'     #### Step2: QC ####
#'     if(Step2_QC){
#'         print('Performing QC of genes and spots...')
#'         st_obj <- QC_Spatial(
#'             st_obj = st_obj,
#'             output.dir = file.path(output.dir, 'Step2_QC'),
#'             min.gene = min.gene,
#'             min.nUMI = min.nUMI,
#'             max.gene = max.gene,
#'             max.nUMI = max.nUMI,
#'             min.spot = min.spot,
#'             species = species,
#'             bool.remove.mito = bool.remove.mito,
#'             SpatialColors = SpatialColors
#'         )
#'     }
#'
#'
#'     #### Step3: Normalization, PCA and Clustering ####
#'     if(Step3_Clustering){
#'         print('Performing normalization, PCA and clustering...')
#'         st_obj <- st_Clustering(
#'             st_obj = st_obj,
#'             output.dir = file.path(output.dir, 'Step3_Clustering'),
#'             normalization.method = normalization.method,
#'             npcs = npcs,
#'             pcs.used = pcs.used,
#'             resolution = resolution,
#'             verbose = verbose
#'         )
#'     }
#'
#'     #### Step4: Differential expressed genes ####
#'     if(Step4_Find_DEGs){
#'         print('Finding differential expressed genes in each cluster...')
#'         st.markers <- st_Find_DEGs(
#'             st_obj = st_obj,
#'             output.dir = file.path(output.dir, 'Step4_Find_DEGs'),
#'             ident.label = 'seurat_clusters',
#'             only.pos = only.pos,
#'             min.pct = min.pct,
#'             logfc.threshold = logfc.threshold,
#'             test.use = test.use,
#'             verbose = verbose
#'         )
#'     }
#'
#'     #### Step5: Spatially variable features ####
#'     if(Step5_SVFs){
#'         print('Identifying spatially variable features...')
#'         st_obj <- st_SpatiallyVariableFeatures(
#'             st_obj = st_obj,
#'             output.dir = file.path(output.dir, 'Step5_SpatiallyVariableFeatures'),
#'             assay = st_obj@active.assay,
#'             selection.method = selection.method,
#'             n.top.show = n.top.show,
#'             n.col = n.col.show,
#'             verbose = verbose
#'         )
#'     }
#'
#'     #### Step6: Spatial interaction ####
#'     if(Step6_Interaction){
#'         print('Performing spatial interaction analysis using Commot...')
#'         interaction_path = file.path(output.dir, 'Step6_Interaction')
#'         if(!dir.exists(interaction_path)){
#'             dir.create(interaction_path)
#'         }
#'         write.csv(st_obj@meta.data,
#'                   file.path(interaction_path, 'metadata.csv'),
#'                   row.names = TRUE)
#'         # st_interaction(
#'         #     st_data_path = input.data.dir,
#'         #     metadata_path = file.path(interaction_path, 'metadata.csv'),
#'         #     label_key = 'seurat_clusters',
#'         #     save_path = interaction_path,
#'         #     species = species,
#'         #     signaling_type = commot.signaling_type,
#'         #     database = commot.database,
#'         #     min_cell_pct = commot.min_cell_pct,
#'         #     dis_thr = commot.dis_thr,
#'         #     n_permutations = commot.n_permutations,
#'         #     condaenv = condaenv
#'         # )
#'         st_Interaction(
#'             st_data_path = file.path(output.dir, 'Step1_Loading_Data'),
#'             metadata_path = file.path(interaction_path, 'metadata.csv'),
#'             library_id = sampleName,
#'             label_key = 'seurat_clusters',
#'             save_path = interaction_path,
#'             species = species,
#'             signaling_type = commot.signaling_type,
#'             database = commot.database,
#'             min_cell_pct = commot.min_cell_pct,
#'             dis_thr = commot.dis_thr,
#'             n_permutations = commot.n_permutations,
#'             condaenv = condaenv
#'         )
#'     }
#'
#'     #### Step7: CNV analysis ####
#'     if(Step7_CNV){
#'         print('Performing CNV analysis using copykat...')
#'         st_obj <- st_CNV(
#'             st_obj = st_obj,
#'             save_path = file.path(output.dir, 'Step7_CNV_analysis'),
#'             assay = assay,
#'             LOW.DR = copykat.LOW.DR,
#'             UP.DR = copykat.UP.DR,
#'             win.size = copykat.win.size,
#'             distance = copykat.distance,
#'             genome = copykat.genome,
#'             n.cores = copykat.n.cores,
#'             species = species
#'         )
#'     }
#'
#'     #### Step8: Deconvolution ####
#'     if(Step8_Deconvolution){
#'         print('Performing deconvolution using cell2location...')
#'         st_obj <- st_Deconvolution(
#'             st.data.dir = file.path(output.dir, 'Step1_Loading_Data'),
#'             library_id = sampleName,
#'             sc.h5ad.dir = cell2loc.sc.h5ad.dir,
#'             st_obj = st_obj,
#'             save_path = file.path(output.dir, 'Step8_Deconvolution'),
#'             sc.labels.key = 'seurat_clusters',
#'             species = species,
#'             sc.max.epoch = cell2loc.sc.max.epoch,
#'             st.max.epoch = cell2loc.st.max.epoch,
#'             use.gpu = cell2loc.use.gpu,
#'             condaenv = condaenv,
#'             use.Dataset = cell2loc.use.dataset
#'         )
#'
#'         if(is.null(cell2loc.sc.h5ad.dir)){
#'             cell2loc.sc.h5ad.dir <- 'Default data'
#'         }
#'     }
#'
#'     #### Step9: Cell cycle analysis ####
#'     if(Step9_Cellcycle){
#'         print('Performing cell cycle analysis...')
#'         st_obj <- st_Cell_cycle(
#'             st_obj = st_obj,
#'             save_path = file.path(output.dir, 'Step9_Cellcycle'),
#'             s.features = s.features,
#'             g2m.features = g2m.features,
#'             species = species,
#'             FeatureColors.bi = FeatureColors.bi
#'         )
#'
#'         g2m.features = unlist(
#'             ifelse(is.null(g2m.features),
#'                    ifelse(species == 'human',
#'                           yes = list(cc.genes.updated.2019$g2m.genes),
#'                           no = list(convertHumanGene(cc.genes.updated.2019$g2m.genes))),
#'                    list(g2m.features))
#'         )
#'         s.features = unlist(
#'             ifelse(is.null(s.features),
#'                    ifelse(species == 'human',
#'                           yes = list(cc.genes.updated.2019$s.genes),
#'                           no = list(convertHumanGene(cc.genes.updated.2019$s.genes))),
#'                    list(s.features))
#'         )
#'     }
#'
#'     #### Step10: Niche analysis ####
#'     if(Step10_Niche){
#'         print('Performing niche analysis...')
#'         if(!file.exists(file.path(output.dir, 'Step8_Deconvolution', 'cell2loc_res.csv'))){
#'             stop('Please run step 8 deconvolution first.')
#'         }
#'         tmp <- read.csv(file.path(output.dir, 'Step8_Deconvolution', 'cell2loc_res.csv'),
#'                         row.names = 1)
#'         features <- colnames(tmp)
#'
#'         if(!all(features %in% names(st_obj@meta.data))){
#'             common.barcodes <- intersect(colnames(st_obj), rownames(tmp))
#'             tmp <- tmp[common.barcodes, ]
#'             st_obj <- st_obj[, common.barcodes]
#'             st_obj <- AddMetaData(st_obj,
#'                                   metadata = tmp)
#'         }
#'
#'         rm(tmp)
#'         gc()
#'
#'         st_obj <- st_NicheAnalysis(
#'             st_obj,
#'             features = features,
#'             save_path = file.path(output.dir, 'Step10_NicheAnalysis'),
#'             coexistence.method = coexistence.method,
#'             kmeans.n = Niche.cluster.n,
#'             st_data_path = file.path(output.dir, 'Step1_Loading_Data'),
#'             slice = slice,
#'             species = species,
#'             condaenv = condaenv
#'         )
#'     }
#'
#'     #### Save data ####
#'     output.dir.final <- file.path(output.dir, 'Data_and_report')
#'     if(!dir.exists(output.dir.final)){
#'         dir.create(output.dir.final)
#'     }
#'     saveRDS(st_obj, file.path(output.dir.final, 'st_object.rds'))
#'     suppressMessages(suppressWarnings(
#'         SaveH5Seurat(st_obj,
#'                      filename = file.path(output.dir.final, 'st_object.h5Seurat'),
#'                      overwrite = TRUE)
#'     ))
#'     suppressMessages(suppressWarnings(
#'         Convert(file.path(output.dir.final, 'st_object.h5Seurat'), dest = "h5ad",
#'                 overwrite = TRUE)
#'     ))
#'
#'     #### Generate the report ####
#'     print('Generating the report...')
#'     if(genReport){
#'         library(kableExtra)
#'         knitr::knit(file.path(system.file(package = "HemaScopeR"), "rmd/st_base.Rmd"),
#'                     file.path(output.dir.final, 'st_pipeline.md'))
#'         markdown::markdownToHTML(file.path(output.dir.final, 'st_pipeline.md'),
#'                                  file.path(output.dir.final, 'st_pipeline.html'))
#'     }
#' }
#'
#'
#'
#' # ### Regress cell-cycle ###
#' # if(bool.regress.cycling){
#' #
#' #     print('Regressing out cell cycle scores...')
#' #
#' #     st_obj@active.assay <- 'Spatial'
#' #     st_obj <- st_obj %>%
#' #         NormalizeData(assay = 'Spatial', verbose = verbose) %>%
#' #         FindVariableFeatures(verbose = verbose)
#' #     if(bool.regress.cycling.standard){
#' #         st_obj <- ScaleData(st_obj,
#' #                             vars.to.regress = c('S.Score', 'G2M.Score'),
#' #                             features = rownames(st_obj),
#' #                             verbose = verbose)
#' #     }else{
#' #         st_obj$CC.Difference <- st_obj$S.Score - st_obj$G2M.Score
#' #         st_obj <- ScaleData(st_obj,
#' #                             vars.to.regress = 'CC.Difference',
#' #                             features = rownames(st_obj),
#' #                             verbose = verbose)
#' #     }
#' #     st_obj <- st_obj %>%
#' #         RunPCA(assay = 'Spatial', verbose = verbose) %>%
#' #         FindNeighbors(reduction = "pca", dims = 1:n.dim.used, verbose = verbose) %>%
#' #         FindClusters(resolution = resolution, verbose = verbose) %>%
#' #         RunUMAP(reduction = "pca", dims = 1:n.dim.used, verbose = verbose)
#' #
#' #     n.cluster.regress <- length(unique(st_obj$seurat_clusters))
#' #
#' #     ### Visualization ###
#' #     suppressMessages(suppressWarnings(
#' #         p.cluster.regress.spatial <- SpatialDimPlot(st_obj,
#' #                                                     group.by = 'seurat_clusters',
#' #                                                     stroke = NA) +
#' #             scale_fill_manual(name = 'Cluster',
#' #                               values = getDefaultClusterColor(n.cluster.regress)) +
#' #             theme(legend.position = 'right',
#' #                   legend.key = element_blank()) +
#' #             guides(fill=guide_legend(override.aes = list(size=4)))
#' #     ))
#' #     saveImage(output.dir.figure,
#' #               p.cluster.regress.spatial,
#' #               'cluster_regress_spatial',
#' #               height = 4,
#' #               width = 4)
#' #
#' #     p.cluster.regress.umap <- DimPlot(st_obj, reduction = 'umap') +
#' #         scale_color_manual(name = 'Cluster',
#' #                            values = getDefaultClusterColor(n.cluster.regress)) +
#' #         theme(legend.position = 'right',
#' #               legend.key = element_blank()) +
#' #         guides(color=guide_legend(override.aes = list(size=4)))
#' #     saveImage(output.dir.figure,
#' #               p.cluster.regress.umap,
#' #               'cluster_regress_UMAP',
#' #               height = 4,
#' #               width = 4)
#' #
#' #     ### Differential expressed genes ###
#' #     st_obj@active.ident <- st_obj$seurat_clusters
#' #     st_obj.markers <- FindAllMarkers(st_obj, only.pos = TRUE,
#' #                                      min.pct = 0.25, logfc.threshold = 0.25,
#' #                                      verbose = verbose)
#' #     st_obj.markers.top5 <- st_obj.markers %>% group_by(cluster) %>%
#' #         top_n(n = 5, wt = .data[[grep('avg_log', colnames(st_obj.markers), value = T)]])
#' #     st_obj.markers.top5 <- st_obj.markers.top5[!duplicated(st_obj.markers.top5$gene), ]
#' #     p.degs.regress.dot <- DotPlot(st_obj, features = st_obj.markers.top5$gene,
#' #                                   cols=c("lightgrey",'red'),
#' #                                   group.by = 'seurat_clusters') +
#' #         theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5)) +
#' #         scale_y_discrete(limits = rev(levels(st_obj$seurat_clusters)))
#' #     saveImage(output.dir.figure,
#' #               p.degs.regress.dot,
#' #               'DEGs_regress_dot',
#' #               height = 3+0.2*n.cluster.regress,
#' #               width = 2+0.25*nrow(st_obj.markers.top5))
#' # }

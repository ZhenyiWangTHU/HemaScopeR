#' Comprehensive Pipeline for Analyzing 10x Single-Cell RNA-Seq Data
#'
#' The function \code{scRNASeq_10x_pipeline} provides a comprehensive analysis
#' pipeline for 10x single-cell RNA-Seq data within the HemascopeR package.
#'
#' @param input.data.dirs A character vector containing directory paths where 10x single-cell RNA-Seq data is stored. Each element corresponds to a path for a dataset.
#' @param project.names A character vector specifying project names for each dataset in input.data.dirs. Its length should match the number of datasets in input.data.dirs, providing project names for each dataset. Please note that the 'project.names' parameter should be the same as the 'id' parameter in cellranger count.
#' @param output.dir The directory where analysis results and reports will be stored.
#' @param pythonPath The path to the Python environment to use for the analysis.
#' @param databasePath The path to the database required for the analysis.
#' @param gene.column The column number in the 10x data where gene information is stored. (Default: 2)
#' @param min.cells The minimum number of cells for a gene to be considered. (Default: 10)
#' @param min.feature The minimum number of features (genes) for a cell to be considered. (Default: 200)
#' @param mt.pattern The regular expression pattern to identify mitochondrial genes. (Default: '^mt-')
#' @param nFeature_RNA.limit The minimum number of features (genes) in the cell. (Default: 200)
#' @param percent.mt.limit The maximum percentage of mitochondrial genes allowed. (Default: 20)
#' @param scale.factor The factor for data scaling. (Default: 10000)
#' @param nfeatures The number of highly variable features (genes) to consider. (Default: 3000)
#' @param ndims The number of dimensions for downstream analysis. (Default: 50)
#' @param vars.to.regress Variables to regress during data normalization. (Default: NULL)
#' @param PCs The principal components (PCs) to use for downstream analysis. (Default: 1:35)
#' @param resolution The resolution parameter for Louvain clustering. (Default: 0.4)
#' @param n.neighbors The number of neighbors for downstream analysis. (Default: 50)
#' @param doublet.percentage The percentage used for doublet removal. (Default: 0.04)
#' @param doublerFinderwraper.PCs The PCs used in the doublet removal process. (Default: 1:20)
#' @param doublerFinderwraper.pN The pN parameter for doublet removal. (Default: 0.25)
#' @param doublerFinderwraper.pK The pK parameter for doublet removal. (Default: 0.1)
#' @param phate.knn The number of nearest neighbors for PhateR analysis. (Default: 50)
#' @param phate.npca The number of principal components for PhateR analysis. (Default: 20)
#' @param phate.t The t parameter for PhateR analysis. (Default: 10)
#' @param phate.ndim The number of dimensions for PhateR analysis. (Default: 2)
#' @param min.pct The minimum percentage for a gene to be considered in differential gene detection. (Default: 0.25)
#' @param marker.genes A character vector specifying marker genes for violin plot visualization. (Default: NULL)
#' @param logfc.threshold The log-fold change threshold for differential gene detection. (Default: 0.25)
#' @param ViolinPlot.cellTypeOrders The order of cell types for violin plot visualization. (Default: NULL)
#' @param ViolinPlot.cellTypeColors The colors for cell types in violin plot visualization. (Default: NULL)
#' @param Org The organism for analysis (e.g., 'mmu' for mouse, 'hsa' for human). (Default: NULL)
#' @param lineage.genelist A list of gene names for defining cell lineages. (Default: NULL)
#' @param lineage.names A character vector specifying names for each lineage in lineage.genelist. (Default: NULL)
#' @param groups_colors A named list specifying colors for groups. (Default: NULL)
#' @param your_openai_API_key The openai_key parameter for gptcelltype. The OpenAI key obtained from https://platform.openai.com/account/api-keys The default is NA, which will resulting outputing the prompt itself. If an actual key is provided, then the output will be the celltype annotations from the GPT model specified by the user.
#' @param tissuename The tissuename parameter for gptcelltype. Optional input of tissue name. Default = ''.
#' @param gptmodel The model parameter for gptcelltype. A valid GPT-4 or GPT-3.5 model name list on https://platform.openai.com/docs/models. Default is 'gpt-3.5'.
#' @param slingshot.start.clus The starting cluster for slingshot analysis. (Default: NULL)
#' @param slingshot.end.clus The ending cluster for slingshot analysis. (Default: NULL)
#' @param slingshot.colors The colors for data points in slingshot visualization. (Default: NULL)
#' @param loom.files.path The path of loom files for scVelo analysis. (Default: NULL)
#' @param cellcycleCutoff The cell cycle cutoff for cell cycle assignment (Default: NULL).
#' @param sorting Boolean indicating whether the cells were sorted. (Default: FALSE)
#' @param ncores The number of CPU cores to use for parallel processing. (Default: 1)
#' @param Whether_load_previous_results Boolean indicating whether to load previous analysis results. (Default: FALSE)
#' @param Step1_Input_Data Boolean indicating whether to perform Step 1: Input Data. (Default: TRUE)
#' @param Step1_Input_Data.type The type of the input data, e.g. '10x' for cellranger-count software output, 'Seurat' for Seurat object in Rds format, and 'Matrix' for gene expression matrix in Rds format. (Default: NULL)
#' @param Step2_Quality_Control Boolean indicating whether to perform Step 2: Quality Control. (Default: TRUE)
#' @param Step2_Quality_Control.RemoveBatches Boolean indicating whether to remove batches during quality control. (Default: TRUE)
#' @param Step2_Quality_Control.RemoveDoublets Boolean indicating whether to remove doublets during quality control. (Default: TRUE)
#' @param Step3_Clustering Boolean indicating whether to perform Step 3: Clustering. (Default: TRUE)
#' @param Step4_Identify_Cell_Types Boolean indicating whether to perform Step 4: Identify Cell Types. (Default: TRUE)
#' @param Step4_Use_Which_Labels The character specifies which labels to use for cell type identification in Step 4, including 'clustering', 'abcCellmap.1', 'abcCellmap.2', 'abcCellmap.3', 'abcCellmap.4', 'HematoMap', and 'changeLabels'. (Default: 'clustering')
#' @param Step4_Cluster_Labels Characters. If Step4_Use_Which_Labels parameter was set to 'changeLabels', please set Step4_Cluster_Labels parameter as the characters indicating the clustering labels. (Default: NULL)
#' @param Step4_Changed_Labels Characters. If Step4_Use_Which_Labels parameter was set to 'changeLabels', please set Step4_Changed_Labels parameter as the characters indicating the changed labels. Note that the length of the character vectors of Step4_Cluster_Labels should be equal to that of Step4_Changed_Labels. (Default: NULL)
#' @param Step4_run_sc_CNV Boolean indicating whether to run single-cell copy number variation analysis in Step 4. (Default: FALSE)
#' @param Step5_Visualization Boolean indicating whether to perform Step 5: Visualization. (Default: TRUE)
#' @param Step6_Find_DEGs Boolean indicating whether to perform Step 6: Find Differentially Expressed Genes (DEGs). (Default: TRUE)
#' @param Step7_Assign_Cell_Cycle Boolean indicating whether to perform Step 7: Assign Cell Cycle. (Default: TRUE)
#' @param Step8_Calculate_Heterogeneity Boolean indicating whether to perform Step 8: Calculate Heterogeneity. (Default: TRUE)
#' @param Step9_Violin_Plot_for_Marker_Genes Boolean indicating whether to perform Step 9: Violin Plot for Marker Genes. (Default: TRUE)
#' @param Step10_Calculate_Lineage_Scores Boolean indicating whether to perform Step 10: Calculate Lineage Scores. (Default: TRUE)
#' @param Step11_GSVA Boolean indicating whether to perform Step 11: Gene Set Variation Analysis (GSVA). (Default: TRUE)
#' @param Step11_GSVA.identify.cellType.features Boolean indicating whether to identify cell type-specific features in GSVA. (Default: TRUE)
#' @param Step11_GSVA.identify.diff.features Boolean indicating whether to identify differentially expressed features in GSVA. (Default: TRUE)
#' @param Step11_GSVA.comparison.design A list specifying the comparison design for GSVA. (Default: NULL)
#' @param Step12_Construct_Trajectories Boolean indicating whether to perform Step 12: Construct Trajectories. (Default: TRUE)
#' @param Step12_Construct_Trajectories.clusters A vector indicating the clusters used for trajectory construction in Step 12. (Default: NULL)
#' @param Step12_Construct_Trajectories.monocle Boolean indicating whether to use Monocle for trajectory construction in Step 12. (Default: FALSE)
#' @param Step12_Construct_Trajectories.slingshot Boolean indicating whether to use slingshot for trajectory construction in Step 12. (Default: FALSE)
#' @param Step12_Construct_Trajectories.scVelo Boolean indicating whether to use scVelo for trajectory construction in Step 12. (Default: FALSE)
#' @param Step13_TF_Analysis Boolean indicating whether to perform Step 13: Transcription Factor (TF) Analysis. (Default: TRUE)
#' @param Step13_TF_Analysis.cellTypes_colors A character vector indicating cell type-specific colors in TF analysis in Step 13. (Default: NULL)
#' @param Step13_TF_Analysis.groups_colors A character vector indicating group-specific colors in TF analysis in Step 13. (Default: NULL)
#' @param Step14_Cell_Cell_Interaction Boolean indicating whether to perform Step 14: Cell-Cell Interaction. (Default: TRUE)
#' @param Step15_Generate_the_Report Boolean indicating whether to perform Step 15: Generate the Report. (Default: TRUE)
#'
#' @details
#' This pipeline offers a comprehensive workflow for the analysis of 10x single-cell
#' RNA-Seq data. It includes the following steps:
#'
#' 1. Input Data: Loading and preprocessing of single-cell RNA-Seq data from specified
#' directories.
#'
#' 2. Quality Control: Performing data quality control, filtering low-quality cells
#' and features, and calculating metrics like the percentage of mitochondrial genes.
#'
#' 3. Clustering: Identifying cell clusters and visualizing them using t-SNE and UMAP
#' dimensionality reduction.
#'
#' 4. Automated Cell Type Identification: Automatically assigning cell types to clusters
#' based on marker genes and mapping single-cell data to an internal reference atlas of
#' hematopoietic cells.
#'
#' 5. Visualization: Run phateR for additional dimensionality reduction and visualization.
#'
#' 6. Differential Gene Detection: Identifying differentially expressed genes (DEGs) between
#' cell clusters.
#'
#' 7. Cell Cycle Assignment: Assigning cells to specific cell cycle phases.
#'
#' 8. Assessment of Heterogeneity: Assessing heterogeneity within cell populations.
#'
#' 9. Creating violin plots to explore marker genes and cluster characteristics.
#'
#' 10. Prediction of Lineage Scores: Predicting lineage scores for individual cells.
#'
#' 11. Gene Set Variation Analysis (GSVA): Calculating gene set variation scores for cell types.
#'
#' 12. Cell Differentiation Trajectory Prediction: Predicting cell differentiation trajectories
#' using Monocle2, Slingshot and scVelo.
#'
#' 13. Transcription Factor (TF) Analysis: Analyzing transcription factors and their activity.
#'
#' 14. Cell-Cell Interactions: Analyzing cell-cell interactions and communication.
#'
#' 15. Generate a well-formatted HTML analysis report.
#'
#' Note that, if you will use Step 5 or scVelo in Step 12, please ensure that you have installed required Python-packages in your Python environment and specify
#' the Python path using the 'pythonPath' parameter.
#'
#' @return Returns the user with a well-formatted HTML analysis report and organizes the analyzed.
#' results and publication-quality vector images for each step of the analysis workflow.
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
#' @export

scRNASeq_10x_pipeline = function(
                                 # input and output
                                 input.data.dirs = NULL,
                                 project.names = NULL,
                                 output.dir = NULL,
                                 pythonPath = python.path.sc(),
                                 databasePath = NULL,
                                 # quality control and preprocessing
                                 gene.column = 2,
                                 min.cells = 10,
                                 min.feature = 200,
                                 mt.pattern = '^mt-',
                                 nFeature_RNA.limit = 200,
                                 percent.mt.limit = 20,
                                 scale.factor = 10000,
                                 nfeatures = 3000,
                                 ndims = 50,
                                 vars.to.regress = NULL,
                                 PCs = 1:35,
                                 resolution = 0.4,
                                 n.neighbors = 50,
                                 # remove doublets
                                 doublet.percentage = 0.04,
                                 doublerFinderwraper.PCs = 1:20,
                                 doublerFinderwraper.pN = 0.25,
                                 doublerFinderwraper.pK = 0.1,
                                 # phateR
                                 phate.knn = 50,
                                 phate.npca = 20,
                                 phate.t = 10,
                                 phate.ndim = 2,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25,
                                 # visualization
                                 marker.genes = NULL,
                                 ViolinPlot.cellTypeOrders = NULL,
                                 ViolinPlot.cellTypeColors = NULL,
                                 Org = NULL,
                                 lineage.genelist = NULL,
                                 lineage.names = NULL,
                                 groups_colors =NULL,
                                 # GPT cell type annotation
                                 your_openai_API_key = '',
                                 tissuename = '',
                                 gptmodel = 'gpt-3.5',
                                 #slingshot
                                 slingshot.start.clus = NULL,
                                 slingshot.end.clus = NULL,
                                 slingshot.colors = NULL,
                                 loom.files.path = NULL,
                                 # cell cycle
                                 cellcycleCutoff = NULL,
                                 # cell chat
                                 sorting = FALSE,
                                 ncores = 1,
                                 # Verbose = FALSE,
                                 # activeEachStep
                                 Whether_load_previous_results = FALSE,
                                 Step1_Input_Data = TRUE,
                                 Step1_Input_Data.type = NULL,
                                 Step2_Quality_Control = TRUE,
                                 Step2_Quality_Control.RemoveBatches = TRUE,
                                 Step2_Quality_Control.RemoveDoublets = TRUE,
                                 Step3_Clustering = TRUE,
                                 Step4_Identify_Cell_Types = TRUE,
                                 Step4_Use_Which_Labels = NULL,
                                 Step4_Cluster_Labels = NULL,
                                 Step4_Changed_Labels = NULL,
                                 Step4_run_sc_CNV = TRUE,
                                 Step5_Visualization = TRUE,
                                 Step6_Find_DEGs = TRUE,
                                 Step7_Assign_Cell_Cycle = TRUE,
                                 Step8_Calculate_Heterogeneity = TRUE,
                                 Step9_Violin_Plot_for_Marker_Genes = TRUE,
                                 Step10_Calculate_Lineage_Scores = TRUE,
                                 Step11_GSVA = TRUE,
                                 Step11_GSVA.identify.cellType.features=TRUE,
                                 Step11_GSVA.identify.diff.features=FALSE,
                                 Step11_GSVA.comparison.design = NULL,
                                 Step12_Construct_Trajectories = TRUE,
                                 Step12_Construct_Trajectories.clusters = NULL,
                                 Step12_Construct_Trajectories.monocle = TRUE,
                                 Step12_Construct_Trajectories.slingshot = TRUE,
                                 Step12_Construct_Trajectories.scVelo = TRUE,
                                 Step13_TF_Analysis = TRUE,
                                 Step13_TF_Analysis.cellTypes_colors = NULL,
                                 Step13_TF_Analysis.groups_colors = NULL,
                                 Step14_Cell_Cell_Interaction = TRUE,
                                 Step15_Generate_the_Report = TRUE
                                 ){
  wdir <- getwd()

  if(is.null(pythonPath)==FALSE){ reticulate::use_python(pythonPath) }else{print('Please set the path of Python.')}

  if (!file.exists(paste0(output.dir, '/HemaScopeR_results/'))) {
    dir.create(paste0(output.dir, '/HemaScopeR_results/'))
  }

  output.dir <- paste0(output.dir,'/HemaScopeR_results/')

  if (!file.exists(paste0(output.dir, '/RDSfiles/'))) {
    dir.create(paste0(output.dir, '/RDSfiles/'))
  }

  previous_results_path <- paste0(output.dir, '/RDSfiles/')
  if (Whether_load_previous_results) {
       print('Loading the previous results...')
       Load_previous_results(previous_results_path = previous_results_path)
  }

  # Step1. Input data-----------------------------------------------------------------------------
  if(Step1_Input_Data == TRUE){
      print('Step1. Input data.')
      if (!file.exists(paste0(output.dir, '/Step1.Input_data/'))) {
          dir.create(paste0(output.dir, '/Step1.Input_data/'))
        }

      file.copy(from = input.data.dirs, to = paste0(output.dir,'/Step1.Input_data/'), recursive = TRUE)

      if(Step1_Input_Data.type == 'cellranger-count'){
          if(length(input.data.dirs) > 1){
            input.data.list <- c()
            for (i in 1:length(input.data.dirs)) {

                  sc_data.temp <- Read10X(data.dir = input.data.dirs[i],
                                          gene.column = gene.column)
                  sc_object.temp <- CreateSeuratObject(counts = sc_data.temp,
                                                       project = project.names[i],
                                                       min.cells = min.cells,
                                                       min.feature = min.feature)
                  sc_object.temp[["percent.mt"]] <- PercentageFeatureSet(sc_object.temp, pattern = mt.pattern)
                  input.data.list <- c(input.data.list, sc_object.temp)}

          }else{

                  sc_data <- Read10X(data.dir = input.data.dirs,
                                     gene.column = gene.column)
                  sc_object <- CreateSeuratObject(counts = sc_data,
                                                  project = project.names,
                                                  min.cells = min.cells,
                                                  min.feature = min.feature)
                  sc_object[["percent.mt"]] <- PercentageFeatureSet(sc_object, pattern = mt.pattern)

          }
      }else if(Step1_Input_Data.type == 'Seurat'){
          if(length(input.data.dirs) > 1){
            input.data.list <- c()
            for (i in 1:length(input.data.dirs)) {
                  sc_object.temp <- readRDS(input.data.dirs[i])
                  sc_object.temp[["percent.mt"]] <- PercentageFeatureSet(sc_object.temp, pattern = mt.pattern)
                  input.data.list <- c(input.data.list, sc_object.temp)
            }
          }else{
              sc_object <- readRDS(input.data.dirs)
              sc_object[["percent.mt"]] <- PercentageFeatureSet(sc_object, pattern = mt.pattern)
          }
      }else if(Step1_Input_Data.type == 'Matrix'){
            if(length(input.data.dirs) > 1){
            input.data.list <- c()
            for (i in 1:length(input.data.dirs)) {
                  sc_data.temp <- readRDS(input.data.dirs[i])
                  sc_object.temp <- CreateSeuratObject(counts = sc_data.temp,
                                                       project = project.names[i],
                                                       min.cells = min.cells,
                                                       min.feature = min.feature)
                  sc_object.temp[["percent.mt"]] <- PercentageFeatureSet(sc_object.temp, pattern = mt.pattern)
                  input.data.list <- c(input.data.list, sc_object.temp)}

          }else{

                  sc_data <- readRDS(input.data.dirs)
                  sc_object <- CreateSeuratObject(counts = sc_data,
                                                  project = project.names,
                                                  min.cells = min.cells,
                                                  min.feature = min.feature)
                  sc_object[["percent.mt"]] <- PercentageFeatureSet(sc_object, pattern = mt.pattern)

          }
      }else{
          stop('Please input data generated by the cellranger-count software, or a Seurat object, or a gene expression matrix. HemaScopeR does not support other formats of input data.')
      }
            # Get the names of all variables in the current environment
            variable_names <- ls()
            # Loop through the variable names and save them as RDS files
            for (var_name in variable_names) {
              var <- get(var_name)  # Get the variable by its name
              saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
            }
  }else{print('Skip Step1. Input data.')}

  # Step2. Quality Control----------------------------------------------------------------------
  if(Step2_Quality_Control == TRUE){
       print('Step2. Quality control.')
       if (!file.exists(paste0(output.dir, '/Step2.Quality_control/'))) {
          dir.create(paste0(output.dir, '/Step2.Quality_control/'))
       }

       if(length(input.data.dirs) > 1){
        # preprocess and quality control for multiple scRNA-Seq data sets
        sc_object <- QC_multiple_scRNASeq(seuratObjects = input.data.list,
                                          datasetID = project.names,
                                          output.dir = paste0(output.dir,'/Step2.Quality_control/'),
                                          Step2_Quality_Control.RemoveBatches = Step2_Quality_Control.RemoveBatches,
                                          Step2_Quality_Control.RemoveDoublets = Step2_Quality_Control.RemoveDoublets,
                                          nFeature_RNA.limit = nFeature_RNA.limit,
                                          percent.mt.limit = percent.mt.limit,
                                          scale.factor = scale.factor,
                                          nfeatures = nfeatures,
                                          ndims = ndims,
                                          vars.to.regress = vars.to.regress,
                                          PCs = PCs,
                                          resolution = resolution,
                                          n.neighbors = n.neighbors,
                                          percentage = doublet.percentage,
                                          doublerFinderwraper.PCs = doublerFinderwraper.PCs,
                                          doublerFinderwraper.pN = doublerFinderwraper.pN,
                                          doublerFinderwraper.pK = doublerFinderwraper.pK
                                          )

  }else{
        # preprocess and quality control for single scRNA-Seq data set
        sc_object <- QC_single_scRNASeq(sc_object = sc_object,
                                        datasetID = project.names,
                                        output.dir = paste0(output.dir,'/Step2.Quality_control/'),
                                        Step2_Quality_Control.RemoveDoublets = Step2_Quality_Control.RemoveDoublets,
                                        nFeature_RNA.limit = nFeature_RNA.limit,
                                        percent.mt.limit = percent.mt.limit,
                                        scale.factor = scale.factor,
                                        nfeatures = nfeatures,
                                        vars.to.regress = vars.to.regress,
                                        ndims = ndims,
                                        PCs = PCs,
                                        resolution = resolution,
                                        n.neighbors = n.neighbors,
                                        percentage = doublet.percentage,
                                        doublerFinderwraper.PCs = doublerFinderwraper.PCs,
                                        doublerFinderwraper.pN = doublerFinderwraper.pN,
                                        doublerFinderwraper.pK = doublerFinderwraper.pK)
   }

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step2. Quality control.')}

  # Step3. Clustering----------------------------------------------------------------------------------------
  if(Step3_Clustering == TRUE){
      print('Step3. Clustering.')
      if (!file.exists(paste0(output.dir, '/Step3.Clustering/'))) {
        dir.create(paste0(output.dir, '/Step3.Clustering/'))
      }

      if( (length(input.data.dirs) > 1) & Step2_Quality_Control.RemoveBatches ){graph.name <- 'integrated_snn'}else{graph.name <- 'RNA_snn'}
      sc_object <- FindNeighbors(sc_object, dims = PCs, k.param = n.neighbors, force.recalc = TRUE)
      sc_object <- FindClusters(sc_object, resolution = resolution, graph.name = graph.name)
      sc_object@meta.data$seurat_clusters <- as.character(as.numeric(sc_object@meta.data$seurat_clusters))

      # plot clustering
      pdf(paste0(paste0(output.dir,'/Step3.Clustering/'), '/sc_object ','tsne_cluster.pdf'), width = 6, height = 6)
       print(DimPlot(sc_object, reduction = "tsne", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      pdf(paste0(paste0(output.dir,'/Step3.Clustering/'), '/sc_object ','umap_cluster.pdf'), width = 6, height = 6)
       print(DimPlot(sc_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      png(paste0(paste0(output.dir,'/Step3.Clustering/'), '/sc_object ','tsne_cluster.png'), width = 600, height = 600)
       print(DimPlot(sc_object, reduction = "tsne", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      png(paste0(paste0(output.dir,'/Step3.Clustering/'), '/sc_object ','umap_cluster.png'), width = 600, height = 600)
       print(DimPlot(sc_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step3. Clustering.')}

  # Step4. Identify cell types automatically---------------------------------------------------------------------------
  if(Step4_Identify_Cell_Types == TRUE){
      print('Step4. Identify cell types automatically.')
      if (!file.exists(paste0(output.dir, '/Step4.Identify_Cell_Types/'))) {
        dir.create(paste0(output.dir, '/Step4.Identify_Cell_Types/'))
      }

      sc_object <- run_cell_annotation(object = sc_object,
                                       assay = 'RNA',
                                       species = Org,
                                       output.dir = paste0(output.dir,'/Step4.Identify_Cell_Types/'))

      if(Org == 'hsa'){
      load(paste0(databasePath,"/HematoMap.reference.rdata"))
        if(length(intersect(rownames(HematoMap.reference), rownames(sc_object))) < 1000){
              HematoMap.reference <- RenameGenesSeurat(obj = HematoMap.reference,
                                                       newnames = toupper(rownames(HematoMap.reference)),
                                                       gene.use = rownames(HematoMap.reference),
                                                       de.assay = "RNA",
                                                       lassays = "RNA")
          }

        if(sc_object@active.assay == 'integrated'){
            DefaultAssay(sc_object) <- 'RNA'
            sc_object <- mapDataToRef(ref_object = HematoMap.reference,
                                      ref_labels = HematoMap.reference@meta.data$CellType,
                                      query_object = sc_object,
                                      PCs = PCs,
                                      output.dir = paste0(output.dir, '/Step4.Identify_Cell_Types/'))
            DefaultAssay(sc_object) <- 'integrated'
        }else{
            sc_object <- mapDataToRef(ref_object = HematoMap.reference,
                                      ref_labels = HematoMap.reference@meta.data$CellType,
                                      query_object = sc_object,
                                      PCs = PCs,
                                      output.dir = paste0(output.dir, '/Step4.Identify_Cell_Types/'))
        }

      }

      # set the cell labels
      if(Step4_Use_Which_Labels == 'clustering'){
          sc_object@meta.data$selectLabels <- sc_object@meta.data$seurat_clusters
          Idents(sc_object) <- sc_object@meta.data$selectLabels
      }else if(Step4_Use_Which_Labels == 'abcCellmap.1'){
          sc_object@meta.data$selectLabels <- sc_object@meta.data$Seurat.RNACluster
          Idents(sc_object) <- sc_object@meta.data$selectLabels
      }else if(Step4_Use_Which_Labels == 'abcCellmap.2'){
          sc_object@meta.data$selectLabels <- sc_object@meta.data$scmap.RNACluster
          Idents(sc_object) <- sc_object@meta.data$selectLabels
      }else if(Step4_Use_Which_Labels == 'abcCellmap.3'){
          sc_object@meta.data$selectLabels <- sc_object@meta.data$Seurat.Immunophenotype
          Idents(sc_object) <- sc_object@meta.data$selectLabels
      }else if(Step4_Use_Which_Labels == 'abcCellmap.4'){
          sc_object@meta.data$selectLabels <- sc_object@meta.data$scmap.Immunophenotype
          Idents(sc_object) <- sc_object@meta.data$selectLabels
      }else if(Step4_Use_Which_Labels == 'HematoMap'){
          if(Org == 'hsa'){
            sc_object@meta.data$selectLabels <- sc_object@meta.data$predicted.id
            Idents(sc_object) <- sc_object@meta.data$selectLabels
          }else{print("'HematoMap' is only applicable to human data ('Org' = 'hsa').")}
      }else if(Step4_Use_Which_Labels == 'changeLabels'){
          if (!is.null(Step4_Cluster_Labels) && !is.null(Step4_Changed_Labels) && length(Step4_Cluster_Labels) == length(Step4_Changed_Labels)){
           sc_object@meta.data$selectLabels <- plyr::mapvalues(sc_object@meta.data$seurat_clusters,
                                                               from = as.character(Step4_Cluster_Labels),
                                                               to = as.character(Step4_Changed_Labels),
                                                               warn_missing = FALSE)
           Idents(sc_object) <- sc_object@meta.data$selectLabels
          }else{
           print("Please input the 'Step4_Cluster_Labels' parameter as Seurat clustering labels, and the 'Step4_Changed_Labels' parameter as new labels. Please note that these two parameters should be of equal length.")
          }
      }else{
          print('Please set the "Step4_Use_Which_Labels" parameter as "clustering", "abcCellmap.1", "abcCellmap.2", "HematoMap" or "changeLabels".')
      }

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step4. Identify cell types automatically.')}

  if(Step4_run_sc_CNV==TRUE){
          sc_CNV(sc_object=sc_object,
                 save_path=paste0(output.dir,'/Step4.Identify_Cell_Types/'),
                 assay = 'RNA',
                 LOW.DR = 0.05,
                 UP.DR = 0.1,
                 win.size = 25,
                 distance = "euclidean",
                 genome = NULL,
                 n.cores = ncores,
                 species = Org)
        }
  # Step5. Visualization-----------------------------------------------------------------------------------------------
  if(Step5_Visualization == TRUE){
      print('Step5. Visualization.')
      if (!file.exists(paste0(output.dir, '/Step5.Visualization/'))) {
        dir.create(paste0(output.dir, '/Step5.Visualization/'))
      }
      # statistical results
      cells_labels <- as.data.frame(cbind(rownames(sc_object@meta.data), as.character(sc_object@meta.data$selectLabels)))
      colnames(cells_labels) <- c('cell_id', 'cluster_id')
      cluster_counts <- cells_labels %>%
        group_by(cluster_id) %>%
        summarise(count = n())
      total_cells <- nrow(cells_labels)
      cluster_counts <- cluster_counts %>%
        mutate(proportion = count / total_cells)

      cluster_counts <- as.data.frame(cluster_counts)

      cluster_counts$percentages <- scales::percent(cluster_counts$proportion, accuracy = 0.1)
      cluster_counts <- cluster_counts[,-which(colnames(cluster_counts)=='proportion')]

      cluster_counts$cluster_id_count_percentages <- paste(cluster_counts$cluster_id, " (", cluster_counts$count, ' cells; ', cluster_counts$percentages, ")", sep='')

      cluster_counts <- cluster_counts[order(cluster_counts$count, decreasing = TRUE),]

      cluster_counts <- rbind(cluster_counts, c('Total', sum(cluster_counts$count), '100%', 'all cells'))

      sc_object@meta.data$cluster_id_count_percentages <- mapvalues(sc_object@meta.data$selectLabels,
                                                                    from=cluster_counts$cluster_id,
                                                                    to=cluster_counts$cluster_id_count_percentages,
                                                                    warn_missing=FALSE)

      colnames(sc_object@meta.data)[which(colnames(sc_object@meta.data) == 'cluster_id_count_percentages')] <- paste('Total ', nrow(sc_object@meta.data), ' cells', sep='')

      cluster_counts <- cluster_counts[,-which(colnames(cluster_counts)=='cluster_id_count_percentages')]

      colnames(cluster_counts) <- c('Cell types', 'Cell counts', 'Percentages')

      # names(colorvector) <- mapvalues(names(colorvector),
      #                                 from=cluster_counts$cluster_id,
      #                                 to=cluster_counts$cluster_id_count_percentages,
      #                                 warn_missing=FALSE)

      write.csv(cluster_counts, file=paste(paste0(output.dir, '/Step5.Visualization/'), '/cell types_cell counts_percentages.csv', sep=''), quote=FALSE, row.names=FALSE)

      pdf(paste(paste0(output.dir, '/Step5.Visualization/'), '/cell types_cell counts_percentages_umap.pdf', sep=''), width = 14, height = 6)
        print(DimPlot(sc_object,
                      reduction = "umap",
                      group.by = paste('Total ', nrow(sc_object@meta.data), ' cells', sep=''),
                      label = FALSE,
                      pt.size = 0.1,
                      raster = FALSE))
      dev.off()

      # run phateR
      if( (length(input.data.dirs) > 1) & Step2_Quality_Control.RemoveBatches ){
          DefaultAssay(sc_object) <- 'integrated'
      }else{
          DefaultAssay(sc_object) <- 'RNA'}

      if(!is.null(pythonPath)){
        run_phateR(sc_object = sc_object,
                   output.dir = paste0(output.dir,'/Step5.Visualization/'),
                   pythonPath = pythonPath,
                   phate.knn = phate.knn,
                   phate.npca = phate.npca,
                   phate.t = phate.t,
                   phate.ndim = phate.ndim)
      }

      # plot cell types
      pdf(paste0(paste0(output.dir,'/Step5.Visualization/'), '/sc_object ','tsne cell types.pdf'), width = 6, height = 6)
       print(DimPlot(sc_object, reduction = "tsne", group.by = "ident", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      pdf(paste0(paste0(output.dir,'/Step5.Visualization/'), '/sc_object ','umap cell types.pdf'), width = 6, height = 6)
       print(DimPlot(sc_object, reduction = "umap", group.by = "ident", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      png(paste0(paste0(output.dir,'/Step5.Visualization/'), '/sc_object ','tsne cell types.png'), width = 600, height = 600)
       print(DimPlot(sc_object, reduction = "tsne", group.by = "ident", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      png(paste0(paste0(output.dir,'/Step5.Visualization/'), '/sc_object ','umap cell types.png'), width = 600, height = 600)
       print(DimPlot(sc_object, reduction = "umap", group.by = "ident", label = FALSE, pt.size = 0.1, raster = FALSE))
      dev.off()

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step5. Visualization.')}

  # Step6. Find DEGs----------------------------------------------------------------------------------------------------
  if(Step6_Find_DEGs == TRUE){
      print('Step6. Find DEGs.')
      if (!file.exists(paste0(output.dir, '/Step6.Find_DEGs/'))) {
        dir.create(paste0(output.dir, '/Step6.Find_DEGs/'))
      }

     if (!file.exists(paste0(output.dir, '/Step6.Find_DEGs/OpenXGR/'))) {
        dir.create(paste0(output.dir, '/Step6.Find_DEGs/OpenXGR/'))
      }

      sc_object.markers <- FindAllMarkers(sc_object, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
      write.csv(sc_object.markers,
                file = paste0(paste0(output.dir, '/Step6.Find_DEGs/'),'sc_object.markerGenes.csv'),
                quote=FALSE)

      # GO enrichment
      if(Org=='mmu'){
        OrgDb <- 'org.Mm.eg.db'
      }else if(Org=='hsa'){
        OrgDb <- 'org.Hs.eg.db'
      }else{
        stop("Org should be 'mmu' or 'hsa'.")
      }

      GPT_annotation( marker.genes = sc_object.markers,
                      your_openai_API_key = your_openai_API_key,
                      tissuename = tissuename,
                      gptmodel = gptmodel,
                      output.dir = paste0(output.dir, '/Step6.Find_DEGs/'))

      HemaScopeREnrichment(DEGs=sc_object.markers,
                           OrgDb=OrgDb,
                           output.dir=paste0(output.dir, '/Step6.Find_DEGs/'))

      sc_object.markers.top5 <- sc_object.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

      pdf(paste0(paste0(output.dir, '/Step6.Find_DEGs/'), 'sc_object_markerGenesTop5.pdf'),
          width = 0.5*length(unique(sc_object.markers.top5$gene)),
          height = 0.5*length(unique(Idents(sc_object))))
          print(DotPlot(sc_object,
                  features = unique(sc_object.markers.top5$gene),
                  cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1)))
      dev.off()

      png(paste0(paste0(output.dir, '/Step6.Find_DEGs/'), 'sc_object_markerGenesTop5.png'),
          width = 20*length(unique(sc_object.markers.top5$gene)),
          height = 30*length(unique(Idents(sc_object))))
          print(DotPlot(sc_object,
                  features = unique(sc_object.markers.top5$gene),
                  cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1)))
      dev.off()

      OpenXGR_SAG(sc_object.markers = sc_object.markers,
                  output.dir = paste0(output.dir, '/Step6.Find_DEGs/OpenXGR/'),
                  subnet.size = 10)
      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step6. Find DEGs.')}

  # Step7. Assign Cell Cycles---------------------------------------------------------------------------------------------------------
  if(Step7_Assign_Cell_Cycle == TRUE){
      print('Step7. Assign cell cycles.')
      if (!file.exists(paste0(output.dir, '/Step7.Assign_cell_cycles/'))) {
        dir.create(paste0(output.dir, '/Step7.Assign_cell_cycles/'))
      }

      datasets.before.batch.removal <- readRDS(paste0(paste0(output.dir, '/RDSfiles/'),'datasets.before.batch.removal.rds'))
      sc_object <- cellCycle(sc_object=sc_object,
                             counts_matrix = GetAssayData(object = datasets.before.batch.removal, slot = "counts")%>%as.matrix(),
                             data_matrix = GetAssayData(object = datasets.before.batch.removal, slot = "data")%>%as.matrix(),
                             cellcycleCutoff = cellcycleCutoff,
                             cellTypeOrders = unique(sc_object@meta.data$selectLabels),
                             output.dir=paste0(output.dir, '/Step7.Assign_cell_cycles/'),
                             databasePath = databasePath,
                             Org = Org)

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step7. Assign cell cycles.')}

  # Step8. Calculate Heterogeneity--------------------------------------------------------------------------------------------------
  if(Step8_Calculate_Heterogeneity == TRUE){
      print('Step8. Calculate heterogeneity.')
      if (!file.exists(paste0(output.dir, '/Step8.Calculate_heterogeneity/'))) {
        dir.create(paste0(output.dir, '/Step8.Calculate_heterogeneity/'))
      }
      expression_matrix <- GetAssayData(object = datasets.before.batch.removal, slot = "data")%>%as.matrix()
      expression_matrix <- expression_matrix[,rownames(sc_object@meta.data)]
      cell_types_groups <- as.data.frame(cbind(sc_object@meta.data$selectLabels,
                                               sc_object@meta.data$datasetID))
      colnames(cell_types_groups) <- c('clusters', 'datasetID')

      if(is.null(ViolinPlot.cellTypeOrders)){
        cellTypes_orders <- unique(sc_object@meta.data$selectLabels)
      }else{
        cellTypes_orders <- ViolinPlot.cellTypeOrders
      }

      heterogeneity(expression_matrix = expression_matrix,
                    cell_types_groups = cell_types_groups,
                    cellTypeOrders = cellTypes_orders,
                    output.dir = paste0(output.dir, '/Step8.Calculate_heterogeneity/'))

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step8. Calculate heterogeneity.')}

  # Step9. Violin Plot for Marker Genes------------------------------------------------------------------------------------------------------------------
  if(Step9_Violin_Plot_for_Marker_Genes == TRUE){
      print('Step9. Violin plot for marker genes.')
      if (!file.exists(paste0(output.dir, '/Step9.Violin_plot_for_marker_genes/'))) {
        dir.create(paste0(output.dir, '/Step9.Violin_plot_for_marker_genes/'))
      }

      if( (length(input.data.dirs) > 1) & Step2_Quality_Control.RemoveBatches ){
          DefaultAssay(sc_object) <- 'integrated'
      }else{
          DefaultAssay(sc_object) <- 'RNA'}

      dataMatrix <- GetAssayData(object = sc_object, slot = "scale.data")

      if(is.null(marker.genes)&(Org == 'mmu')){
          # mpp genes are from 'The bone marrow microenvironment at single cell resolution'
          # the other genes are from 'single cell characterization of haematopoietic progenitors and their trajectories in homeostasis and perturbed haematopoiesis'
          # the aliases of these genes were changed in gecodeM16：Gpr64 -> Adgrg2, Sdpr -> Cavin2, Hbb-b1 -> Hbb-bs, Sfpi1 -> Spi1
          HSC_lineage_signatures <- c('Slamf1', 'Itga2b', 'Kit', 'Ly6a', 'Bmi1', 'Gata2', 'Hlf', 'Meis1', 'Mpl', 'Mcl1', 'Gfi1', 'Gfi1b', 'Hoxb5')
          Mpp_genes <- c('Mki67', 'Mpo', 'Elane', 'Ctsg', 'Calr')
          Erythroid_lineage_signatures <- c('Klf1', 'Gata1', 'Mpl', 'Epor', 'Vwf', 'Zfpm1', 'Fhl1', 'Adgrg2', 'Cavin2','Gypa', 'Tfrc', 'Hbb-bs', 'Hbb-y')
          Lymphoid_lineage_signatures <- c('Tcf3', 'Ikzf1', 'Notch1', 'Flt3', 'Dntt', 'Btg2', 'Tcf7', 'Rag1', 'Ptprc', 'Ly6a', 'Blnk')
          Myeloid_lineage_signatures <- c('Gfi1', 'Spi1', 'Mpo', 'Csf2rb', 'Csf1r', 'Gfi1b', 'Hk3', 'Csf2ra', 'Csf3r', 'Sp1', 'Fcgr3')
          marker.genes <- c(HSC_lineage_signatures, Mpp_genes, Erythroid_lineage_signatures, Lymphoid_lineage_signatures, Myeloid_lineage_signatures)
      }else if(is.null(marker.genes)&(Org == 'hsa')){
          HSPCs_lineage_signatures <- c('CD34','KIT','AVP','FLT3','MME','CD7','CD38','CSF1R','FCGR1A','MPO','ELANE','IL3RA')
          Myeloids_lineage_signatures <- c('LYZ','CD36','MPO','FCGR1A','CD4','CD14','CD300E','ITGAX','FCGR3A','FLT3','AXL',
                        'SIGLEC6','CLEC4C','IRF4','LILRA4','IL3RA','IRF8','IRF7','XCR1','CD1C','THBD',
                        'MRC1','CD34','KIT','ITGA2B','PF4','CD9','ENG','KLF','TFRC')
          B_cells_lineage_signatures <- c('CD79A','IGLL1','RAG1','RAG2','VPREB1','MME','IL7R','DNTT','MKI67','PCNA','TCL1A','MS4A1','IGHD','CD27','IGHG3')
          T_NK_cells_lineage_signatures <- c('CD3D','CD3E','CD8A','CCR7','IL7R','SELL','KLRG1','CD27','GNLY',
                          'NKG7','PDCD1','TNFRSF9','LAG3','CD160','CD4','CD40LG','IL2RA',
                          'FOXP3','DUSP4','IL2RB','KLRF1','FCGR3A','NCAM1','XCL1','MKI67','PCNA','KLRF')
          marker.genes <- c(HSPCs_lineage_signatures, Myeloids_lineage_signatures, B_cells_lineage_signatures, T_NK_cells_lineage_signatures)
      }

      if(is.null(ViolinPlot.cellTypeOrders)){
        ViolinPlot.cellTypeOrders <- unique(sc_object@meta.data$selectLabels)
      }

      if(is.null(ViolinPlot.cellTypeColors)){
        ViolinPlot.cellTypeColors <- viridis::viridis(length(unique(sc_object@meta.data$selectLabels)))
      }

      combinedViolinPlot(dataMatrix = dataMatrix,
                         features = marker.genes,
                         CellTypes = sc_object@meta.data$selectLabels,
                         cellTypeOrders = ViolinPlot.cellTypeOrders,
                         cellTypeColors = ViolinPlot.cellTypeColors,
                         Org = Org,
                         output.dir = paste0(output.dir, '/Step9.Violin_plot_for_marker_genes/'),
                         databasePath = databasePath)

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step9. Violin plot for marker genes.')}

  # Step10. Calculate Lineage Scores-----------------------------------------------------------
  if( Step10_Calculate_Lineage_Scores == TRUE){
      print('Step10. Calculate lineage scores.')
      # we use normalized data here
      if (!file.exists(paste0(output.dir, '/Step10.Calculate_lineage_scores/'))) {
        dir.create(paste0(output.dir, '/Step10.Calculate_lineage_scores/'))
      }

      if(is.null(lineage.genelist)&is.null(lineage.names)&(Org == 'mmu')){
          lineage.genelist <- c(list(HSC_lineage_signatures),
                                list(Mpp_genes),
                                list(Erythroid_lineage_signatures),
                                list(Lymphoid_lineage_signatures),
                                list(Myeloid_lineage_signatures))
          lineage.names <- c('HSC_lineage_signatures',
                             'Mpp_genes',
                             'Erythroid_lineage_signatures',
                             'Lymphoid_lineage_signatures',
                             'Myeloid_lineage_signatures')
      }else if(is.null(lineage.genelist)&is.null(lineage.names)&(Org == 'hsa')){
          lineage.genelist <- c(list(HSPCs_lineage_signatures),
                                list(Myeloids_lineage_signatures),
                                list(B_cells_lineage_signatures),
                                list(T_NK_cells_lineage_signatures))
          lineage.names <- c('HSPCs_lineage_signatures',
                             'Myeloids_lineage_signatures',
                             'B_cells_lineage_signatures',
                             'T_NK_cells_lineage_signatures')
      }

      if(is.null(ViolinPlot.cellTypeOrders)){
        cellTypes_orders <- unique(sc_object@meta.data$selectLabels)
      }else{
        cellTypes_orders <- ViolinPlot.cellTypeOrders
      }

      lineageScores(expression_matrix = expression_matrix,
                    cellTypes = sc_object@meta.data$selectLabels,
                    cellTypes_orders = cellTypes_orders,
                    cellTypes_colors = ViolinPlot.cellTypeColors,
                    groups = sc_object@meta.data$datasetID,
                    groups_orders = unique(sc_object@meta.data$datasetID),
                    groups_colors = groups_colors,
                    lineage.genelist = lineage.genelist,
                    lineage.names = lineage.names,
                    Org = Org,
                    output.dir = paste0(output.dir, '/Step10.Calculate_lineage_scores/'),
                    databasePath = databasePath)

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step10. Calculate lineage scores.')}

  # Step11. GSVA-------------------------------------------------------------------------------------------
  if(Step11_GSVA == TRUE){
      print('Step11. GSVA.')
      if (!file.exists(paste0(output.dir, '/Step11.GSVA/'))) {
        dir.create(paste0(output.dir, '/Step11.GSVA/'))
      }

      setwd(wdir)

      if(Org=='mmu'){
        load(paste0(databasePath,"/mouse_c2_v5p2.rdata"))
        GSVA.genelist <- Mm.c2
        assign('OrgDB', org.Mm.eg.db)
      }else if(Org=='hsa'){
        load(paste0(databasePath,"/human_c2_v5p2.rdata"))
        GSVA.genelist <- Hs.c2
        assign('OrgDB', org.Hs.eg.db)
      }else{
        stop("Org should be 'mmu' or 'hsa'.")
      }

      if(is.null(ViolinPlot.cellTypeOrders)){
        cellTypes_orders <- unique(sc_object@meta.data$selectLabels)
      }else{
        cellTypes_orders <- ViolinPlot.cellTypeOrders
      }
      run_GSVA(sc_object = sc_object,
               GSVA.genelist = GSVA.genelist,
               GSVA.cellTypes = sc_object@meta.data$selectLabels,
               GSVA.cellTypes.orders = cellTypes_orders,
               GSVA.cellGroups = sc_object@meta.data$datasetID,
               GSVA.identify.cellType.features = Step11_GSVA.identify.cellType.features,
               GSVA.identify.diff.features = Step11_GSVA.identify.diff.features,
               GSVA.comparison.design = Step11_GSVA.comparison.design,
               OrgDB = OrgDB,
               output.dir = paste0(output.dir, '/Step11.GSVA/'))

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step11. GSVA.')}

  # Step12. Construct Trajectories-------------------------------------------------------------
  #data("genecode_geneSymbolandEnsembleID")

  DefaultAssay(sc_object) <- 'RNA'
  countsSlot <- GetAssayData(object = sc_object, slot = "counts")
  gene_metadata <- as.data.frame(rownames(countsSlot))
  rownames(gene_metadata) <- gene_metadata[,1]
  if(Org == 'mmu'){
     load(paste0(databasePath,"/mouseGeneSymbolandEnsembleID.rdata"))
     gene_metadata $ ensembleID <- mapvalues(x = gene_metadata[,1],
                                             from = mouseGeneSymbolandEnsembleID$geneName,
                                             to = mouseGeneSymbolandEnsembleID$ensemblIDNoDot,
                                             warn_missing = FALSE)
  }else if(Org == 'hsa'){
     load(paste0(databasePath,"/humanGeneSymbolandEnsembleID.rdata"))
     gene_metadata $ ensembleID <- mapvalues(x = gene_metadata[,1],
                                             from = humanGeneSymbolandEnsembleID$geneName,
                                             to = humanGeneSymbolandEnsembleID$ensemblIDNoDot,
                                             warn_missing = FALSE)
  }

  colnames(gene_metadata) <- c('gene_short_name','ensembleID')

  if(Step12_Construct_Trajectories == TRUE){
      print('Step12. Construct trajectories.')
      if (!file.exists(paste0(output.dir, '/Step12.Construct_trajectories/'))) {
        dir.create(paste0(output.dir, '/Step12.Construct_trajectories/'))
      }

      if(is.null(Step12_Construct_Trajectories.clusters)){
          sc_object.subset <- sc_object
          countsSlot.subset <- GetAssayData(object = sc_object.subset, slot = "counts")
      }else{
          sc_object.subset <- subset(sc_object, subset = selectLabels %in% Step12_Construct_Trajectories.clusters)
          countsSlot.subset <- GetAssayData(object = sc_object.subset, slot = "counts")
      }

      if(Step12_Construct_Trajectories.monocle == TRUE){
          # monocle2
          if (!file.exists(paste0(output.dir, '/Step12.Construct_trajectories/monocle2/'))) {
            dir.create(paste0(output.dir, '/Step12.Construct_trajectories/monocle2/'))
          }

          phenoData <- sc_object.subset@meta.data
          featureData <- gene_metadata
          run_monocle(cellData = countsSlot.subset,
                      phenoData = phenoData,
                      featureData = featureData,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = VGAM::negbinomial.size(),
                      cellTypes='selectLabels',
                      monocle.orders=Step12_Construct_Trajectories.clusters,
                      monocle.colors = ViolinPlot.cellTypeColors,
                      output.dir = paste0(output.dir, '/Step12.Construct_trajectories/monocle2/'))
      }

      if(Step12_Construct_Trajectories.slingshot == TRUE){
          # slingshot
          if (!file.exists(paste0(output.dir, '/Step12.Construct_trajectories/slingshot/'))) {
            dir.create(paste0(output.dir, '/Step12.Construct_trajectories/slingshot/'))
          }

          if( (length(input.data.dirs) > 1) & Step2_Quality_Control.RemoveBatches ){
              DefaultAssay(sc_object.subset) <- 'integrated'
          }else{
              DefaultAssay(sc_object.subset) <- 'RNA'}
          run_slingshot(slingshot.PCAembeddings = Embeddings(sc_object.subset, reduction = "pca")[, PCs],
                        slingshot.cellTypes = sc_object.subset@meta.data$selectLabels,
                        slingshot.start.clus = slingshot.start.clus,
                        slingshot.end.clus = slingshot.end.clus,
                        slingshot.colors = slingshot.colors,
                        output.dir = paste0(output.dir, '/Step12.Construct_trajectories/slingshot/'))
      }

      if(Step12_Construct_Trajectories.scVelo == TRUE){
          # scVelo
          if((!is.null(loom.files.path))&(!is.null(pythonPath))){
              if (!file.exists(paste0(output.dir, '/Step12.Construct_trajectories/scVelo/'))) {
                dir.create(paste0(output.dir, '/Step12.Construct_trajectories/scVelo/'))
              }

              prepareDataForScvelo(sc_object = sc_object.subset,
                                   loom.files.path = loom.files.path,
                                   scvelo.reduction = 'pca',
                                   scvelo.column = 'selectLabels',
                                   output.dir = paste0(output.dir, '/Step12.Construct_trajectories/scVelo/'))

              reticulate::py_run_string(paste0("import os\noutputDir = '", output.dir, "'"))
              reticulate::py_run_file(file.path(system.file(package = "HemaScopeR"), "python/sc_run_scvelo.py"), convert = FALSE)
          }
      }
      # URD
      # sc_object.markers.top50 <- sc_object.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
      # sc_object.markers.top50 <- sc_object.markers.top50$gene%>%unique()
      # URD.meta <- sc_object@meta.data$selectLabels %>% as.data.frame()
      # URD.meta$sampleNames <- colnames(countsSlot)
      # rownames(URD.meta) <- URD.meta[,2]
      # colnames(URD.meta) <- c('label', 'sampleName')
      # save.image('./beforeURD.RData')
      # run_URD(URD.count.data = as.matrix(countsSlot),
      #         URD.meta = URD.meta,
      #         URD.select.deatures = sc_object.markers.top50,
      #         root.cellTypes=NULL,
      #         tips.cellTypes=NULL,
      #         output.dir = output.dir)

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step12. Construct trajectories.')}

  # Step13. TF Analysis------------------------------------------------------------------------------------
  if(Step13_TF_Analysis == TRUE){
      print('Step13. TF analysis.')
      if (!file.exists(paste0(output.dir, '/Step13.TF_analysis/'))) {
        dir.create(paste0(output.dir, '/Step13.TF_analysis/'))
      }

      run_SCENIC(countMatrix = countsSlot,
                 cellTypes = sc_object@meta.data$selectLabels,
                 datasetID = sc_object@meta.data$datasetID,
                 cellTypes_colors = Step13_TF_Analysis.cellTypes_colors,
                 cellTypes_orders = unique(sc_object@meta.data$selectLabels),
                 groups_colors = Step13_TF_Analysis.groups_colors,
                 groups_orders = unique(sc_object@meta.data$datasetID),
                 Org = Org,
                 output.dir = paste0(output.dir, '/Step13.TF_analysis/'),
                 pythonPath = pythonPath,
                 databasePath = databasePath)

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step13. TF analysis.')}

  # Step14. Cell-Cell Interaction----------------------------------------------------------------------------
  if( Step14_Cell_Cell_Interaction == TRUE){
      print('Step14. Cell-cell interaction.')
      if (!file.exists(paste0(output.dir, '/Step14.Cell_cell_interection/'))) {
        dir.create(paste0(output.dir, '/Step14.Cell_cell_interection/'))
      }
      tempwd <- getwd()
      run_CellChat(data.input=countsSlot,
                   labels = sc_object@meta.data$selectLabels,
                   cell.orders = ViolinPlot.cellTypeOrders,
                   cell.colors = ViolinPlot.cellTypeColors,
                   sample.names = rownames(sc_object@meta.data),
                   Org = Org,
                   sorting = sorting,
                   output.dir = paste0(output.dir, '/Step14.Cell_cell_interection/'))
      setwd(tempwd)
      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
  }else{print('Skip Step14. Cell-cell interaction.')}

  # Step15. Generate the Report-----------------------------------------------------------------------------
  if(Step15_Generate_the_Report == TRUE){
      print('Step15. Generate the report.')
  # generate the report
  }else{print('Skip Step15. Generate the report.')}
}

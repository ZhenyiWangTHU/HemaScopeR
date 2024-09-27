# server---------------------------------------------------------------------------------------------------------------------------------
server = function(input, output, session){
  output$logo <- renderImage({
    list(src = '../images/hemascoper_logo.png')
  }, deleteFile = FALSE)

  output$ui_styles <- renderUI({
    tags$style(HTML("
      #logo img {
        max-width: 500px; 
        height: auto;
        margin-left: -0px;
        margin-bottom: -350px;
      }
      .h1-font {
        font-size: 35px;
      }
      .main-panel {
        margin-top: 100px; /* 调整上边距 */
        margin-bottom: 100px;
      }
      .spacer {
        margin-bottom: 10px; /* 缩小空行的间距 */
      }
    "))
  })

  output$stepContent <- renderUI({
        if(step() == 'sc_step1') {
          step1_fluidRow
        }else if(step() == 'sc_step2') {
          step2_fluidRow
        }else if(step() == 'sc_step3') {
          step3_fluidRow
        }else if(step() == 'sc_step4') {
          step4_fluidRow
        }else if(step() == 'sc_step5') {
          step5_fluidRow
        }else if(step() == 'sc_step6') {
          step6_fluidRow
        }else if(step() == 'sc_step7') {
          step7_fluidRow
        }else if(step() == 'sc_step8') {
          step8_fluidRow
        }else if(step() == 'sc_step9') {
          step9_fluidRow
        }else if(step() == 'sc_step10') {
          step10_fluidRow
        }else if(step() == 'sc_step11') {
          step11_fluidRow
        }else if(step() == 'sc_step12') {
          step12_fluidRow
        }else if(step() == 'sc_step13') {
          step13_fluidRow
        }else if(step() == 'sc_step14') {
          step14_fluidRow
        }else if(step() == 'sc_step15') {
          step15_fluidRow
        }
  })

  output$stepContent_st <- renderUI({
        if(step() == 'st_step1') {
          step1_fluidRow_st
        }else if(step() == 'st_step2') {
          step2_fluidRow_st
        }else if(step() == 'st_step3') {
          step3_fluidRow_st
        }else if(step() == 'st_step4') {
          step4_fluidRow_st
        }else if(step() == 'st_step5') {
          step5_fluidRow_st
        }else if(step() == 'st_step6') {
          step6_fluidRow_st
        }else if(step() == 'st_step7') {
          step7_fluidRow_st
        }else if(step() == 'st_step8') {
          step8_fluidRow_st
        }else if(step() == 'st_step9') {
          step9_fluidRow_st
        }else if(step() == 'st_step10') {
          step10_fluidRow_st
        }else if(step() == 'st_step11') {
          step11_fluidRow_st
        }
  })

 ### scRNA-seq pipeline---------------------------------------------------------------------------------------------------------
  # Start----------------------------------------------------------------------------  
  observeEvent(input$start_button, {
    shinyjs::hide("ui1")
    shinyjs::show("ui2")
    #updateTabsetPanel(session, "mainTabs2", selected = "homeTab")  # ui2 tab
  })

  step <- reactiveVal('sc_step1')
  input.data.dirs.temp <- reactiveVal(NULL)
  project.names.temp <- reactiveVal(NULL)
  output.dir.temp <- reactiveVal(NULL)
  pythonPath.temp <- reactiveVal(NULL)
  databasePath.temp <- reactiveVal(NULL)
  # quality control and preprocessing
  gene.column.temp <- reactiveVal(NULL)
  min.cells.temp <- reactiveVal(NULL)
  min.feature.temp <- reactiveVal(NULL)
  mt.pattern.temp <- reactiveVal(NULL)
  nFeature_RNA.limit.temp <- reactiveVal(NULL)
  percent.mt.limit.temp <- reactiveVal(NULL)
  scale.factor.temp <- reactiveVal(NULL)
  nfeatures.temp <- reactiveVal(NULL)
  ndims.temp <- reactiveVal(NULL)
  vars.to.regress.temp <- reactiveVal(NULL)
  PCs.temp <- reactiveVal(NULL)
  PCs.clustering.temp <- reactiveVal(NULL)
  resolution.temp <- reactiveVal(NULL)
  n.neighbors.temp <- reactiveVal(NULL)
  # remove doublets
  doublet.percentage.temp <- reactiveVal(NULL)
  doublerFinderwraper.PCs.temp <- reactiveVal(NULL)
  doublerFinderwraper.pN.temp <- reactiveVal(NULL)
  doublerFinderwraper.pK.temp <- reactiveVal(NULL)
  # phateR
  phate.knn.temp <- reactiveVal(NULL)
  phate.npca.temp <- reactiveVal(NULL)
  phate.t.temp <- reactiveVal(NULL)
  phate.ndim.temp <- reactiveVal(NULL)
  min.pct.temp <- reactiveVal(NULL)
  logfc.threshold.temp <- reactiveVal(NULL)
  # visualization
  marker.genes.temp <- reactiveVal(NULL)
  ViolinPlot.cellTypeOrders.temp <- reactiveVal(NULL)
  ViolinPlot.cellTypeColors.temp <- reactiveVal(NULL)
  Org.temp <- reactiveVal(NULL)
  lineage.genelist.temp <- reactiveVal(NULL)
  lineage.names.temp <- reactiveVal(NULL)
  groups_colors.temp <- reactiveVal(NULL)
  slingshot.start.clus.temp <- reactiveVal(NULL)
  slingshot.end.clus.temp <- reactiveVal(NULL)
  slingshot.colors.temp <- reactiveVal(NULL)
  loom.files.path.temp <- reactiveVal(NULL)
  # cell cycle
  cellcycleCutoff.temp <- reactiveVal(NULL)
  # cell chat
  sorting.temp <- reactiveVal(NULL)
  ncores.temp <- reactiveVal(NULL)
  # # Verbose = FALSE,
  # # activeEachStep
  # Step1_Input_Data = TRUE,
  Step1_Input_Data.type.temp <- reactiveVal(NULL)
  # Step2_Quality_Control = TRUE,
  Step2_Quality_Control.RemoveBatches.temp <- reactiveVal(NULL)
  Step2_Quality_Control.RemoveDoublets.temp <- reactiveVal(NULL)
  # Step3_Clustering = TRUE,
  # Step4_Identify_Cell_Types = TRUE,
  Step4_Use_Which_Labels.temp <- reactiveVal(NULL)
  Step4_Cluster_Labels.temp <- reactiveVal(NULL)
  Step4_Changed_Labels.temp <- reactiveVal(NULL)
  Step4_run_sc_CNV.temp <- reactiveVal(NULL)
  # Step5_Visualization = TRUE,
  # Step6_Find_DEGs = TRUE,
  # Step7_Assign_Cell_Cycle = TRUE,
  # Step8_Calculate_Heterogeneity = TRUE,
  # Step9_Violin_Plot_for_Marker_Genes = TRUE,
  # Step10_Calculate_Lineage_Scores = TRUE,
  # Step11_GSVA = TRUE,
  # Step11_GSVA.identify.cellType.features=TRUE,
  # Step11_GSVA.identify.diff.features=FALSE,
  Step11_GSVA.identify.cellType.features.temp <- reactiveVal(NULL)
  Step11_GSVA.identify.diff.features.temp <- reactiveVal(NULL)
  Step11_GSVA.comparison.design.temp <- reactiveVal(NULL)
  # Step12_Construct_Trajectories = TRUE,
  Step12_Construct_Trajectories.clusters.temp <- reactiveVal(NULL)
  Step12_Construct_Trajectories.monocle.temp <- reactiveVal(NULL)
  Step12_Construct_Trajectories.slingshot.temp <- reactiveVal(NULL)
  Step12_Construct_Trajectories.scVelo.temp <- reactiveVal(NULL)
  # Step12_Construct_Trajectories.monocle = TRUE,
  # Step12_Construct_Trajectories.slingshot = TRUE,
  # Step12_Construct_Trajectories.scVelo = TRUE,
  # Step13_TF_Analysis = TRUE,
  Step13_TF_Analysis.cellTypes_colors.temp <- reactiveVal(NULL)
  Step13_TF_Analysis.groups_colors.temp <- reactiveVal(NULL)
  # Step14_Cell_Cell_Interaction = TRUE,
  # Step15_Generate_the_Report = TRUE

  # Step1. Input data-------------------------------------------------------------------------------------
  previous_results_path.temp <- reactiveVal(NULL)
  #sc_object.temp <- reactiveVal(NULL)
  observeEvent(input$step1, {
    step('sc_step1')
  })

  observeEvent(input$load_data_button, {
      # set parameters
      input.data.dirs <- input$input.data.dirs
      project.names <- input$project.names
      output.dir <- input$output.dir
      pythonPath <- input$pythonPath
      databasePath <- input$databasePath
      gene.column <- input$gene.column
      Step1_Input_Data.type <- input$Step1_Input_Data.type
      gene.column <- input$gene.column
      min.cells <- input$min.cells
      min.feature <- input$min.feature
      mt.pattern <- input$mt.pattern
      # convert some parameters
      if(pythonPath=="NULL"){pythonPath <- NULL}else{pythonPath <- pythonPath}
      if(input.data.dirs=='NULL'){input.data.dirs <- NULL}else{
        input.data.dirs <- unlist(strsplit(input.data.dirs, ";"))
        input.data.dirs <- gsub("^([^/].*)$", "/\\1", input.data.dirs)
      }

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
      
      shinyjs::runjs('$("#loadingData").text("Loading data...");')

      # output.dir.temp(output.dir)
      if (!is.null(input.data.dirs) && file.exists(input.data.dirs)){
        shinyjs::disable("input.data.dirs")
        shinyjs::disable("project.names")
        shinyjs::disable("output.dir")
        shinyjs::disable("pythonPath")
        shinyjs::disable("databasePath")
        shinyjs::disable("gene.column")
        shinyjs::disable("Step1_Input_Data.type")
        shinyjs::disable("gene.column")
        shinyjs::disable("min.cells")
        shinyjs::disable("min.feature")
        shinyjs::disable("mt.pattern")
        shinyjs::disable("load_data_button")

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

        # update variables
        previous_results_path.temp(previous_results_path)
        # update parameters
        input.data.dirs.temp(input.data.dirs)
        project.names.temp(project.names)
        output.dir.temp(output.dir)
        pythonPath.temp(pythonPath)
        databasePath.temp(databasePath)
        gene.column.temp(gene.column)
        Step1_Input_Data.type.temp(Step1_Input_Data.type)
        gene.column.temp(gene.column)
        min.cells.temp(min.cells)
        min.feature.temp(min.feature)
        mt.pattern.temp(mt.pattern)
        # Update output with dimension information
        output$data_dim_output <- renderText({
            paste0("OK! Data dimensions: ", paste0(dim(sc_object)[2], ' cells, ', dim(sc_object)[1], ' genes.'))
        })

        shinyjs::enable("input.data.dirs")
        shinyjs::enable("project.names")
        shinyjs::enable("output.dir")
        shinyjs::enable("pythonPath")
        shinyjs::enable("databasePath")
        shinyjs::enable("gene.column")
        shinyjs::enable("Step1_Input_Data.type")
        shinyjs::enable("gene.column")
        shinyjs::enable("min.cells")
        shinyjs::enable("min.feature")
        shinyjs::enable("mt.pattern")
        shinyjs::enable("load_data_button")
      }
  })

  # Step2. Quality control---------------------------------------------------------------------------------
  observeEvent(input$step2, {
    step('sc_step2')
  })

  observeEvent(input$RunStep2, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    input.data.dirs <- input.data.dirs.temp()
    output.dir <- output.dir.temp()
    project.names <- project.names.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    nFeature_RNA.limit <- input$nFeature_RNA.limit
    percent.mt.limit <- input$percent.mt.limit
    scale.factor <- input$scale.factor
    nfeatures <- input$nfeatures
    ndims <- input$ndims
    vars.to.regress <- input$vars.to.regress
    PCs <- input$PCs
    resolution <- input$resolution
    n.neighbors <- input$n.neighbors
    doublet.percentage <- input$doublet.percentage
    doublerFinderwraper.PCs <- input$doublerFinderwraper.PCs
    doublerFinderwraper.pN <- input$doublerFinderwraper.pN
    doublerFinderwraper.pK <- input$doublerFinderwraper.pK
    Step2_Quality_Control.RemoveBatches <- input$Step2_Quality_Control.RemoveBatches
    Step2_Quality_Control.RemoveDoublets <- input$Step2_Quality_Control.RemoveDoublets

    # convert some parameters
    PCs <- seq(from = as.numeric(unlist(strsplit(PCs, ":")))[1], to = as.numeric(unlist(strsplit(PCs, ":")))[2])
    doublerFinderwraper.PCs <- seq(from = as.numeric(unlist(strsplit(doublerFinderwraper.PCs, ":")))[1], to = as.numeric(unlist(strsplit(doublerFinderwraper.PCs, ":")))[2])
    if(vars.to.regress=='NULL'){vars.to.regress <- NULL}else{vars.to.regress <- vars.to.regress}

    shinyjs::runjs('$("#runningStep2").text("Running Step2...");')  
    shinyjs::disable('nFeature_RNA.limit')
    shinyjs::disable('percent.mt.limit')
    shinyjs::disable('scale.factor')
    shinyjs::disable('nfeatures')
    shinyjs::disable('ndims')
    shinyjs::disable('vars.to.regress')
    shinyjs::disable('PCs')
    shinyjs::disable('resolution')
    shinyjs::disable('n.neighbors')
    shinyjs::disable('doublet.percentage')
    shinyjs::disable('doublerFinderwraper.PCs')
    shinyjs::disable('doublerFinderwraper.pN')
    shinyjs::disable('doublerFinderwraper.pK')
    shinyjs::disable('Step2_Quality_Control.RemoveBatches')
    shinyjs::disable('Step2_Quality_Control.RemoveDoublets')
    shinyjs::disable('RunStep2')

    if (!file.exists(paste0(output.dir, '/Step2.Quality_control/'))) {
          dir.create(paste0(output.dir, '/Step2.Quality_control/'))
    }
    Step2_Quality_Control.RemoveBatches <- as.logical(Step2_Quality_Control.RemoveBatches)
    Step2_Quality_Control.RemoveDoublets <- as.logical(Step2_Quality_Control.RemoveDoublets)

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

    # update parameters
    nFeature_RNA.limit.temp(nFeature_RNA.limit)
    percent.mt.limit.temp(percent.mt.limit)
    scale.factor.temp(scale.factor)
    nfeatures.temp(nfeatures)
    ndims.temp(ndims)
    vars.to.regress.temp(vars.to.regress)
    PCs.temp(PCs)
    resolution.temp(resolution)
    n.neighbors.temp(n.neighbors)
    doublet.percentage.temp(doublet.percentage)
    doublerFinderwraper.PCs.temp(doublerFinderwraper.PCs)
    doublerFinderwraper.pN.temp(doublerFinderwraper.pN)
    doublerFinderwraper.pK.temp(doublerFinderwraper.pK)
    Step2_Quality_Control.RemoveBatches.temp(Step2_Quality_Control.RemoveBatches)
    Step2_Quality_Control.RemoveDoublets.temp(Step2_Quality_Control.RemoveDoublets)

    shinyjs::enable('nFeature_RNA.limit')
    shinyjs::enable('percent.mt.limit')
    shinyjs::enable('scale.factor')
    shinyjs::enable('nfeatures')
    shinyjs::enable('ndims')
    shinyjs::enable('vars.to.regress')
    shinyjs::enable('PCs')
    shinyjs::enable('resolution')
    shinyjs::enable('n.neighbors')
    shinyjs::enable('doublet.percentage')
    shinyjs::enable('doublerFinderwraper.PCs')
    shinyjs::enable('doublerFinderwraper.pN')
    shinyjs::enable('doublerFinderwraper.pK')
    shinyjs::enable('Step2_Quality_Control.RemoveBatches')
    shinyjs::enable('Step2_Quality_Control.RemoveDoublets')
    shinyjs::enable('RunStep2')
    output$step2_completed <- renderText({'Step2 completed.'})

    output$Step2.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step2.Quality_control/')  # 指定的文件夹路径
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = paste0(folder_path, file)
            )
          })
          
          tagList(
            p("Files in the folder:"),
            tags$ul(
              lapply(file_links, function(link) {
                tags$li(link)
              })
            )
          )
        } else {
          p("No files found in the folder")
        }
      } else {
        p("Folder does not exist.")
      }
    })

    selected_file <- reactiveVal(NULL)

    # Listen to the click event on the file links
    observe({
      for (file in list.files(paste0(output.dir, '/Step2.Quality_control/'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step2.Quality_control/', file))
        })
      }
    })

    output$download <- downloadHandler(
      filename = function() {
        if (!is.null(selected_file())) {
          basename(selected_file())
        }
      },
      content = function(file) {
        if (!is.null(selected_file())) {
          file_content <- readBin(selected_file(), "raw", n=file.size(selected_file()))
          return(file_content)
        }
      }
    )

  })

  # Step3. Clustering-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step3, {
    step('sc_step3')
  })

  observeEvent(input$RunStep3, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()
      input.data.dirs <- input.data.dirs.temp()
      Step2_Quality_Control.RemoveBatches <- Step2_Quality_Control.RemoveBatches.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      PCs.clustering <- input$PCs.clustering
      n.neighbors <- input$n.neighbors
      resolution <- input$resolution
      # convert some parameters
      PCs.clustering <- seq(from = as.numeric(unlist(strsplit(PCs.clustering, ":")))[1], to = as.numeric(unlist(strsplit(PCs.clustering, ":")))[2])

      shinyjs::runjs('$("#runningStep3").text("Running Step3...");')  
      shinyjs::disable('PCs.clustering')
      shinyjs::disable('n.neighbors')
      shinyjs::disable('resolution')
      shinyjs::disable('RunStep3')

      if (!file.exists(paste0(output.dir, '/Step3.Clustering/'))) {
        dir.create(paste0(output.dir, '/Step3.Clustering/'))
      }

      if( (length(input.data.dirs) > 1) & Step2_Quality_Control.RemoveBatches ){graph.name <- 'integrated_snn'}else{graph.name <- 'RNA_snn'}
      sc_object <- FindNeighbors(sc_object, dims = PCs.clustering, k.param = n.neighbors, force.recalc = TRUE)
      sc_object <- FindClusters(sc_object, resolution = resolution, graph.name = graph.name)
      sc_object@meta.data$seurat_clusters <- as.character(as.numeric(sc_object@meta.data$seurat_clusters))

      # plot clustering
      pdf(paste0(paste0(output.dir,'/Step3.Clustering/'), '/sc_object ','tsne_cluster.pdf'), width = 6, height = 6)
       print(DimPlot(sc_object, reduction = "tsne", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1))
      dev.off()

      pdf(paste0(paste0(output.dir,'/Step3.Clustering/'), '/sc_object ','umap_cluster.pdf'), width = 6, height = 6)
       print(DimPlot(sc_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1))
      dev.off()

      png(paste0(paste0(output.dir,'/Step3.Clustering/'), '/sc_object ','tsne_cluster.png'), width = 600, height = 600)
       print(DimPlot(sc_object, reduction = "tsne", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1))
      dev.off()

      png(paste0(paste0(output.dir,'/Step3.Clustering/'), '/sc_object ','umap_cluster.png'), width = 600, height = 600)
       print(DimPlot(sc_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1))
      dev.off()    
      
      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }  
      
      #update parameters
      PCs.clustering.temp(PCs.clustering)
      n.neighbors.temp(n.neighbors)
      resolution.temp(resolution)

      shinyjs::enable('PCs.clustering')
      shinyjs::enable('n.neighbors')
      shinyjs::enable('resolution')
      shinyjs::enable('RunStep3')
      output$step3_completed <- renderText({'Step3 completed.'})
  })
    
  # Step4. Identify cell types automatically-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step4, {
    step('sc_step4')
  })
    
  observeEvent(input$RunStep4, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()
      PCs <- PCs.temp()
      databasePath <- databasePath.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      Org <- input$Org
      Step4_Use_Which_Labels <- input$Step4_Use_Which_Labels
      Step4_run_sc_CNV <- as.logical(input$Step4_run_sc_CNV)
      ncores <- input$ncores

      shinyjs::runjs('$("#runningStep4").text("Running Step4...");')  
      shinyjs::disable('Org')
      shinyjs::disable('Step4_Use_Which_Labels')
      shinyjs::disable('Step4_run_sc_CNV')
      shinyjs::disable('ncores')
      shinyjs::disable('RunStep4')

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

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }  

      #update parameters
      Org.temp(Org)
      Step4_Use_Which_Labels.temp(Step4_Use_Which_Labels)
      Step4_run_sc_CNV.temp(Step4_run_sc_CNV)
      ncores.temp(ncores)

      shinyjs::enable('Org')
      shinyjs::enable('Step4_Use_Which_Labels')
      shinyjs::enable('Step4_run_sc_CNV')
      shinyjs::enable('ncores')
      shinyjs::enable('RunStep4')
      output$step4_completed <- renderText({'Step4 completed.'})
  })
  
  # Step5. Visualization-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step5, {
    step('sc_step5')
  })

  observeEvent(input$RunStep5, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()
      input.data.dirs <- input.data.dirs.temp()
      Step2_Quality_Control.RemoveBatches <- Step2_Quality_Control.RemoveBatches.temp()
      pythonPath <- pythonPath.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      phate.knn <- input$phate.knn
      phate.npca <- input$phate.npca
      phate.t <- input$phate.t
      phate.ndim <- input$phate.ndim

      shinyjs::runjs('$("#runningStep5").text("Running Step5...");')  
      shinyjs::disable('phate.knn')
      shinyjs::disable('phate.npca')
      shinyjs::disable('phate.t')
      shinyjs::disable('phate.ndim')
      shinyjs::disable('RunStep5')

      if (!file.exists(paste0(output.dir, '/Step5.Visualization/'))) {
        dir.create(paste0(output.dir, '/Step5.Visualization/'))
      }    
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
       print(DimPlot(sc_object, reduction = "tsne", group.by = "ident", label = FALSE, pt.size = 0.1))
      dev.off()

      pdf(paste0(paste0(output.dir,'/Step5.Visualization/'), '/sc_object ','umap cell types.pdf'), width = 6, height = 6)
       print(DimPlot(sc_object, reduction = "umap", group.by = "ident", label = FALSE, pt.size = 0.1))
      dev.off()

      png(paste0(paste0(output.dir,'/Step5.Visualization/'), '/sc_object ','tsne cell types.png'), width = 600, height = 600)
       print(DimPlot(sc_object, reduction = "tsne", group.by = "ident", label = FALSE, pt.size = 0.1))
      dev.off()

      png(paste0(paste0(output.dir,'/Step5.Visualization/'), '/sc_object ','umap cell types.png'), width = 600, height = 600)
       print(DimPlot(sc_object, reduction = "umap", group.by = "ident", label = FALSE, pt.size = 0.1))
      dev.off()    

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }  
      
      #update parameters
      phate.knn.temp(phate.knn)
      phate.npca.temp(phate.npca)
      phate.t.temp(phate.t)
      phate.ndim.temp(phate.ndim)

      shinyjs::enable('phate.knn')
      shinyjs::enable('phate.npca')
      shinyjs::enable('phate.t')
      shinyjs::enable('phate.ndim')
      shinyjs::enable('RunStep5')
      output$step5_completed <- renderText({'Step5 completed.'})
  })
    
  # Step6. Find DEGs-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step6, {
    step('sc_step6')
  })
    
  observeEvent(input$RunStep6, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()
      Org <- Org.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      min.pct <- input$min.pct
      logfc.threshold <- input$logfc.threshold

      shinyjs::runjs('$("#runningStep6").text("Running Step6...");')  
      shinyjs::disable('min.pct')
      shinyjs::disable('logfc.threshold')
      shinyjs::disable('RunStep6')

      if (!file.exists(paste0(output.dir, '/Step6.Find_DEGs/'))) {
        dir.create(paste0(output.dir, '/Step6.Find_DEGs/'))
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

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }  
      
      #update parameters
      min.pct.temp(min.pct)
      logfc.threshold.temp(logfc.threshold)

      shinyjs::enable('min.pct')
      shinyjs::enable('logfc.threshold')
      shinyjs::enable('RunStep6')
      output$step6_completed <- renderText({'Step6 completed.'})
  })
    
  # Step7. Assign Cell Cycles-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step7, {
    step('sc_step7')
  })
      # log_file <- file(paste0(output.dir, '/log.txt'), open = "a")
      # sink(log_file, type = "output")
      # sink(NULL)
  observeEvent(input$RunStep7, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      databasePath <- databasePath.temp()
      output.dir <- output.dir.temp()
      Org <- Org.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      cellcycleCutoff <- input$cellcycleCutoff

      # convert some parameters
      if(cellcycleCutoff=='NULL'){cellcycleCutoff <- NULL}else{cellcycleCutoff <- as.numeric(cellcycleCutoff)}

      shinyjs::runjs('$("#runningStep7").text("Running Step7...");')  
      shinyjs::disable('cellcycleCutoff')
      shinyjs::disable('RunStep7')

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
      
      #update parameters
      cellcycleCutoff.temp(cellcycleCutoff)

      shinyjs::enable('cellcycleCutoff')
      shinyjs::enable('RunStep7')
      output$step7_completed <- renderText({'Step7 completed.'})
  })  
  
  # Step8. Calculate Heterogeneity-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step8, {
    step('sc_step8')
  })
    
  observeEvent(input$RunStep8, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      ViolinPlot.cellTypeOrders <- input$ViolinPlot.cellTypeOrders

      # convert some parameters
      if(ViolinPlot.cellTypeOrders=='NULL'){ViolinPlot.cellTypeOrders <- NULL}else{ViolinPlot.cellTypeOrders <- unlist(strsplit(ViolinPlot.cellTypeOrders, ","))}

      shinyjs::runjs('$("#runningStep8").text("Running Step8...");')
      shinyjs::disable('ViolinPlot.cellTypeOrders')
      shinyjs::disable('RunStep8')

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
      
      #update parameters
      ViolinPlot.cellTypeOrders.temp(ViolinPlot.cellTypeOrders)

      shinyjs::enable('ViolinPlot.cellTypeOrders')
      shinyjs::enable('RunStep8')
      output$step8_completed <- renderText({'Step8 completed.'})
  })    
    
  # Step9. Violin Plot for Marker Genes-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step9, {
    step('sc_step9')
  })
    
  observeEvent(input$RunStep9, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()
      input.data.dirs <- input.data.dirs.temp()
      databasePath <- databasePath.temp()
      Step2_Quality_Control.RemoveBatches <- Step2_Quality_Control.RemoveBatches.temp()
      Org <- Org.temp()
      ViolinPlot.cellTypeOrders <- ViolinPlot.cellTypeOrders.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      marker.genes <- input$marker.genes
      ViolinPlot.cellTypeColors <- input$ViolinPlot.cellTypeColors

      # convert some parameters
      if(marker.genes=='NULL'){marker.genes <- NULL}else{marker.genes <- unlist(strsplit(marker.genes, ","))}
      if(ViolinPlot.cellTypeColors=='NULL'){ViolinPlot.cellTypeColors <- NULL}else{ViolinPlot.cellTypeColors <- unlist(strsplit(ViolinPlot.cellTypeColors, ","))}

      shinyjs::runjs('$("#runningStep9").text("Running Step9...");')  
      shinyjs::disable('marker.genes')
      shinyjs::disable('ViolinPlot.cellTypeColors')
      shinyjs::disable('RunStep9')

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
      
      #update parameters
      marker.genes.temp(marker.genes)
      ViolinPlot.cellTypeColors.temp(ViolinPlot.cellTypeColors)

      shinyjs::enable('marker.genes')
      shinyjs::enable('ViolinPlot.cellTypeColors')
      shinyjs::enable('RunStep9')
      output$step9_completed <- renderText({'Step9 completed.'})
  })      

  # Step10. Calculate Lineage Scores-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step10, {
    step('sc_step10')
  })
    
  observeEvent(input$RunStep10, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()
      databasePath <- databasePath.temp()
      ViolinPlot.cellTypeOrders <- ViolinPlot.cellTypeOrders.temp()
      ViolinPlot.cellTypeColors <- ViolinPlot.cellTypeColors.temp()
      Org <- Org.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      lineage.genelist <- input$lineage.genelist
      lineage.names <- input$lineage.names
      groups_colors <- input$groups_colors
      # convert some parameters
      if(lineage.genelist=='NULL'){lineage.genelist <- NULL}else{
        lineage.genelist <- strsplit(lineage.genelist, ";")[[1]]
        lineage.genelist <- lapply(lineage.genelist, function(x) unlist(strsplit(gsub("\"", "", x), ",")))
      }
      if(lineage.names=='NULL'){lineage.names <- NULL}else{lineage.names <- unlist(strsplit(lineage.names, ","))}
      if(groups_colors=='NULL'){groups_colors <- NULL}else{groups_colors <- unlist(strsplit(groups_colors, ","))}

      shinyjs::runjs('$("#runningStep3").text("Running Step3...");')  
      shinyjs::disable('lineage.genelist')
      shinyjs::disable('lineage.names')
      shinyjs::disable('groups_colors')
      shinyjs::disable('RunStep10')

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
      
      #update parameters
      lineage.genelist.temp(lineage.genelist)
      lineage.names.temp(lineage.names)
      groups_colors.temp(groups_colors)

      shinyjs::enable('lineage.genelist')
      shinyjs::enable('lineage.names')
      shinyjs::enable('groups_colors')
      shinyjs::enable('RunStep10')
      output$step10_completed <- renderText({'Step10 completed.'})
  })
  
  # Step11. GSVA-----------------------------------------------------------------------------------------------------------------------------   
  observeEvent(input$step11, {
    step('sc_step11')
  })
    
  observeEvent(input$RunStep11, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      databasePath <- databasePath.temp()
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      ViolinPlot.cellTypeOrders <- ViolinPlot.cellTypeOrders.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      Step11_GSVA.identify.cellType.features <- as.logical(input$Step11_GSVA.identify.cellType.features)
      Step11_GSVA.identify.diff.features <- as.logical(input$Step11_GSVA.identify.diff.features)
      Step11_GSVA.comparison.design <- input$Step11_GSVA.comparison.design
      # convert some parameters
      if(Step11_GSVA.comparison.design=='NULL'){Step11_GSVA.comparison.design <- NULL}else{
        Step11_GSVA.comparison.design <- eval(parse(text = Step11_GSVA.comparison.design))
      }

      shinyjs::runjs('$("#runningStep11").text("Running Step11...");')
      shinyjs::disable('Step11_GSVA.identify.cellType.features')
      shinyjs::disable('Step11_GSVA.identify.diff.features')
      shinyjs::disable('Step11_GSVA.comparison.design')
      shinyjs::disable('RunStep11')

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
      
      #update parameters
      Step11_GSVA.identify.cellType.features.temp(Step11_GSVA.identify.cellType.features)
      Step11_GSVA.identify.diff.features.temp(Step11_GSVA.identify.diff.features)
      Step11_GSVA.comparison.design.temp(Step11_GSVA.comparison.design)

      shinyjs::enable('Step11_GSVA.identify.cellType.features')
      shinyjs::enable('Step11_GSVA.identify.diff.features')
      shinyjs::enable('Step11_GSVA.comparison.design')
      shinyjs::enable('RunStep11')
      output$step11_completed <- renderText({'Step11 completed.'})
  })

  # Step12. Construct Trajectories-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step12, {
    step('sc_step12')
  })
    
  observeEvent(input$RunStep12, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      input.data.dirs <- input.data.dirs.temp()
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      pythonPath <- pythonPath.temp()
      databasePath <- databasePath.temp()
      ViolinPlot.cellTypeColors <- ViolinPlot.cellTypeColors.temp()
      Step2_Quality_Control.RemoveBatches <- Step2_Quality_Control.RemoveBatches.temp()
      PCs <- PCs.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      Step12_Construct_Trajectories.clusters <- input$Step12_Construct_Trajectories.clusters
      Step12_Construct_Trajectories.monocle <- as.logical(input$Step12_Construct_Trajectories.monocle)
      Step12_Construct_Trajectories.slingshot <- as.logical(input$Step12_Construct_Trajectories.slingshot)
      slingshot.start.clus <- input$slingshot.start.clus
      slingshot.end.clus <- input$slingshot.end.clus
      slingshot.colors <- input$slingshot.colors
      Step12_Construct_Trajectories.scVelo <- as.logical(input$Step12_Construct_Trajectories.scVelo)
      loom.files.path <- input$loom.files.path
      # convert some parameters
      if(Step12_Construct_Trajectories.clusters=='all'){Step12_Construct_Trajectories.clusters <- NULL}else{Step12_Construct_Trajectories.clusters <- unlist(strsplit(Step12_Construct_Trajectories.clusters, ","))}
      if(slingshot.start.clus=='NULL'){slingshot.start.clus <- NULL}else{slingshot.start.clus <- unlist(strsplit(slingshot.start.clus, ","))}
      if(slingshot.end.clus=='NULL'){slingshot.end.clus <- NULL}else{slingshot.end.clus <- unlist(strsplit(slingshot.end.clus, ","))}
      if(slingshot.colors=='NULL'){slingshot.colors <- NULL}else{slingshot.colors <- unlist(strsplit(slingshot.colors, ","))}
      if(loom.files.path=='NULL'){loom.files.path <- NULL}else{
        loom.files.path <- unlist(strsplit(loom.files.path, ";"))
        loom.files.path <- gsub("^([^/].*)$", "./\\1", loom.files.path)
      }

      shinyjs::runjs('$("#runningStep12").text("Running Step12...");')
      shinyjs::disable('Step12_Construct_Trajectories.clusters')
      shinyjs::disable('Step12_Construct_Trajectories.monocle')
      shinyjs::disable('Step12_Construct_Trajectories.slingshot')
      shinyjs::disable('slingshot.start.clus')
      shinyjs::disable('slingshot.end.clus')
      shinyjs::disable('slingshot.colors')
      shinyjs::disable('Step12_Construct_Trajectories.scVelo')
      shinyjs::disable('loom.files.path')

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
              reticulate::py_run_file("~/HemaScopeR/R/sc_run_scvelo.py", convert = FALSE)
          }
      }

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }  
      
      #update parameters
      Step12_Construct_Trajectories.clusters.temp(Step12_Construct_Trajectories.clusters)
      Step12_Construct_Trajectories.monocle.temp(Step12_Construct_Trajectories.monocle) 
      Step12_Construct_Trajectories.slingshot.temp(Step12_Construct_Trajectories.slingshot)
      slingshot.start.clus.temp(slingshot.start.clus)
      slingshot.end.clus.temp(slingshot.end.clus)
      slingshot.colors.temp(slingshot.colors)
      Step12_Construct_Trajectories.scVelo.temp(Step12_Construct_Trajectories.scVelo)
      loom.files.path.temp(loom.files.path)

      shinyjs::enable('Step12_Construct_Trajectories.clusters')
      shinyjs::enable('Step12_Construct_Trajectories.monocle')
      shinyjs::enable('Step12_Construct_Trajectories.slingshot')
      shinyjs::enable('slingshot.start.clus')
      shinyjs::enable('slingshot.end.clus')
      shinyjs::enable('slingshot.colors')
      shinyjs::enable('Step12_Construct_Trajectories.scVelo')
      shinyjs::enable('loom.files.path')
      output$step12_completed <- renderText({'Step12 completed.'})
  })
    
  # Step13. TF Analysis-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step13, {
    step('sc_step13')
  })
    
  observeEvent(input$RunStep13, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      pythonPath <- pythonPath.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      Step13_TF_Analysis.cellTypes_colors <- input$Step13_TF_Analysis.cellTypes_colors
      Step13_TF_Analysis.groups_colors <- input$Step13_TF_Analysis.groups_colors

      shinyjs::runjs('$("#runningStep13").text("Running Step13...");')  
      shinyjs::disable('Step13_TF_Analysis.cellTypes_colors')
      shinyjs::disable('Step13_TF_Analysis.groups_colors')
      shinyjs::disable('RunStep13')

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
                 pythonPath = pythonPath)

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }  
      
      #update parameters
      Step13_TF_Analysis.cellTypes_colors.temp(Step13_TF_Analysis.cellTypes_colors)
      Step13_TF_Analysis.groups_colors.temp(Step13_TF_Analysis.groups_colors)

      shinyjs::enable('Step13_TF_Analysis.cellTypes_colors')
      shinyjs::enable('Step13_TF_Analysis.groups_colors')
      shinyjs::enable('RunStep13')
      output$step13_completed <- renderText({'Step13 completed.'})
  })  
  
  # Step14. -----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step14, {
    step('sc_step14')
  })
    
  observeEvent(input$RunStep14, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      ViolinPlot.cellTypeOrders <- ViolinPlot.cellTypeOrders.temp()
      ViolinPlot.cellTypeColors <- ViolinPlot.cellTypeColors.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
      sorting <- as.logical(input$sorting)

      shinyjs::runjs('$("#runningStep14").text("Running Step14...");')  
      shinyjs::disable('sorting')
      shinyjs::disable('RunStep14')

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
      
      #update parameters
      sorting.temp(sorting)

      shinyjs::enable('sorting')
      shinyjs::enable('RunStep14')
      output$step14_completed <- renderText({'Step14 completed.'})

      output$Step14.file_list <- renderUI({
        # dirs
        if (dir.exists(paste0(output.dir, '/Step14.Cell_cell_interection/'))) {
          # files
          files <- list.files(paste0(output.dir, '/Step14.Cell_cell_interection/'))
          
          if (length(files) > 0) {
            file_links <- lapply(files, function(file) {
              a(href = paste0(output.dir, '/Step14.Cell_cell_interection/', file), file) # create links
            })
            
            tagList(
              p("Files in Step14.Cell_cell_interection:"),
              tags$ul(
                lapply(file_links, function(link) {
                  tags$li(link)
                })
              )
            )
          } else {
            p("No files found in Step14.Cell_cell_interection.")
          }
        } else {
          p("Directory Step14.Cell_cell_interection does not exist.")
        }
      })
  })

  # Step15. -----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$step15, {
    step('sc_step15')
  })
    
  observeEvent(input$RunStep15, {
      # load previous variables
      previous_results_path <- previous_results_path.temp()

      # load previous parameters
      output.dir <- output.dir.temp()

      # load previous results
      Load_previous_results(previous_results_path = previous_results_path)

      # set parameters
  
      shinyjs::runjs('$("#runningStep15").text("Running Step15...");')  
      #shinyjs::disable('')
      shinyjs::disable('RunStep15')



      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }  
      
      #update parameters

      shinyjs::enable('RunStep15')
      output$step15_completed <- renderText({'Step15 completed.'})
  })

  observeEvent(input$back_button, {
    shinyjs::hide("ui2")
    shinyjs::show("ui1")
  })
###-------------------------------------------------------------------------------------------------------------------------------

### st-seq pipeline---------------------------------------------------------------------------------------------------------------
  observeEvent(input$start_button_st, {
    shinyjs::hide("ui1")
    shinyjs::show("ui3")
  })

  step <- reactiveVal('st_step1')
  input.data.dir.temp <- reactiveVal(NULL)
  output.dir.temp <- reactiveVal(NULL)
  sampleName.temp <- reactiveVal(NULL)

  # For Step2 QC
  #Step2_QC.temp <- reactiveVal(NULL)
  min.gene.temp <- reactiveVal(NULL)
  min.nUMI.temp <- reactiveVal(NULL)
  max.gene.temp <- reactiveVal(NULL)
  max.nUMI.temp <- reactiveVal(NULL)
  min.spot.temp <- reactiveVal(NULL)
  bool.remove.mito.temp <- reactiveVal(NULL)
  species.temp <- reactiveVal(NULL)

  # For Step3 Clustering
  #Step3_Clustering.temp <- reactiveVal(NULL)
  normalization.method.temp <- reactiveVal(NULL)
  npcs.temp <- reactiveVal(NULL)
  pcs.used.temp <- reactiveVal(NULL)
  resolution.temp <- reactiveVal(NULL)

  # For Step4 Find DEGs
  #Step4_Find_DEGs.temp <- reactiveVal(NULL)
  only.pos.temp <- reactiveVal(NULL)
  min.pct.temp <- reactiveVal(NULL)
  logfc.threshold.temp <- reactiveVal(NULL)
  test.use.temp <- reactiveVal(NULL)

  # For Step5 SVF
  #Step5_SVFs.temp <- reactiveVal(NULL)
  selection.method.temp <- reactiveVal(NULL)
  n.top.show.temp <- reactiveVal(NULL)
  n.col.show.temp <- reactiveVal(NULL)

  # For Step6 Interaction
  #Step6_Interaction.temp <- reactiveVal(NULL)
  commot.signaling_type.temp <- reactiveVal(NULL)
  commot.database.temp <- reactiveVal(NULL)
  commot.min_cell_pct.temp <- reactiveVal(NULL)
  commot.dis_thr.temp <- reactiveVal(NULL)
  commot.n_permutations.temp <- reactiveVal(NULL)

  # For Step7 CNV analysis
  #Step7_CNV.temp <- reactiveVal(NULL)
  copykat.genome.temp <- reactiveVal(NULL)
  copykat.LOW.DR.temp <- reactiveVal(NULL)
  copykat.UP.DR.temp <- reactiveVal(NULL)
  copykat.win.size.temp <- reactiveVal(NULL)
  copykat.distance.temp <- reactiveVal(NULL)
  copykat.n.cores.temp <- reactiveVal(NULL)

  # For Step8 Deconvolution
  #Step8_Deconvolution.temp <- reactiveVal(NULL)
  cell2loc.sc.h5ad.dir.temp <- reactiveVal(NULL)
  cell2loc.sc.max.epoch.temp <- reactiveVal(NULL)
  cell2loc.st.max.epoch.temp <- reactiveVal(NULL)
  cell2loc.use.gpu.temp <- reactiveVal(NULL)

  # For Step9 Cellcycle
  #Step9_Cellcycle.temp <- reactiveVal(NULL)
  s.features.temp <- reactiveVal(NULL)
  g2m.features.temp <- reactiveVal(NULL)

  # settings
  #condaenv = 'r-reticulate',
  #verbose.temp <- reactiveVal(NULL)
  #genReport.temp <- reactiveVal(NULL)

  previous_results_path.temp <- reactiveVal(NULL)
  pythonPath.temp <- reactiveVal(NULL)

  # Step1. Loading data---------------------------------------------------------------------------------
  observeEvent(input$step1_st, {
    step('st_step1')
  })

  observeEvent(input$RunStep1_st, {

    # set parameters
    input.data.dir <- input$input.data.dir
    output.dir <- input$output.dir
    sampleName <- input$sampleName
    pythonPath <- input$pythonPath

    if(!dir.exists(file.path(output.dir, sampleName))){
        dir.create(file.path(output.dir, sampleName))
    }else{
        warning(paste0('The new results will overwrite the file under ',
                       file.path(output.dir, sampleName)))
    }
    output.dir <- file.path(output.dir, sampleName)
    #print(paste0('The results will be saved in ', output.dir))
    if (!file.exists(paste0(output.dir, '/RDSfiles/'))) {
        dir.create(paste0(output.dir, '/RDSfiles/'))
      }
 
    previous_results_path <- paste0(output.dir, '/RDSfiles/')

    shinyjs::runjs('$("#loadingData").text("Loading data...");')
    shinyjs::disable('input.data.dir')
    shinyjs::disable('output.dir')
    shinyjs::disable('sampleName')
    shinyjs::disable('pythonPath')
    shinyjs::disable('load_data_button')

    #### Step1: Loading data ####
    #print('Loading data...')
    st_obj <- Load10X_Spatial(input.data.dir,
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "slice1",
                              filter.matrix = TRUE,
                              to.upper = FALSE,
                              image = NULL)
    st_obj@project.name <- sampleName
    st_obj$orig.ident <- factor(sampleName)
    st_obj@active.ident <- st_obj$orig.ident

   
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
    }

    # update parameters
    input.data.dir.temp(input.data.dir)
    output.dir.temp(output.dir)
    sampleName.temp(sampleName)
    pythonPath.temp(pythonPath)
    previous_results_path.temp(previous_results_path)

    shinyjs::enable('input.data.dir')
    shinyjs::enable('output.dir')
    shinyjs::enable('sampleName')
    shinyjs::enable('pythonPath')
    shinyjs::enable('load_data_button')
    output$step1_completed <- renderText({'Step1 completed.'})

  })

  # Step2. QC---------------------------------------------------------------------------------
  observeEvent(input$step2_st, {
    step('st_step2')
  })

  observeEvent(input$RunStep2_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    output.dir <- output.dir.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    min.gene <- input$min.gene
    min.nUMI <- input$min.nUMI
    max.gene <- input$max.gene
    max.nUMI <- input$max.nUMI
    min.spot <- input$min.spot
    bool.remove.mito <- as.logical(input$bool.remove.mito)
    species <- input$species

    # convert some parameters
    #

    shinyjs::runjs('$("#runningStep2").text("Running Step2...");')  
    shinyjs::disable('min.gene')
    shinyjs::disable('min.nUMI')
    shinyjs::disable('max.gene')
    shinyjs::disable('max.nUMI')
    shinyjs::disable('min.spot')
    shinyjs::disable('bool.remove.mito')
    shinyjs::disable('species')
    shinyjs::disable('RunStep2')

    SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
    st_obj <- QC_Spatial(
        st_obj = st_obj,
        output.dir = file.path(output.dir, 'Step2_QC'),
        min.gene = min.gene,
        min.nUMI = min.nUMI,
        max.gene = max.gene,
        max.nUMI = max.nUMI,
        min.spot = min.spot,
        species = species,
        bool.remove.mito = bool.remove.mito,
        SpatialColors = SpatialColors
    )


   
      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    min.gene.temp(min.gene)
    min.nUMI.temp(min.nUMI)
    max.gene.temp(max.gene)
    max.nUMI.temp(max.nUMI)
    min.spot.temp(min.spot)
    bool.remove.mito.temp(bool.remove.mito)
    species.temp(species)
    RunStep2.temp(RunStep2)

    shinyjs::enable('min.gene')
    shinyjs::enable('min.nUMI')
    shinyjs::enable('max.gene')
    shinyjs::enable('max.nUMI')
    shinyjs::enable('min.spot')
    shinyjs::enable('bool.remove.mito')
    shinyjs::enable('species')
    shinyjs::enable('RunStep2')
    output$step2_completed <- renderText({'Step2 completed.'})

  })

  # Step3. Normalization, PCA and Clustering---------------------------------------------------------------------------------
  observeEvent(input$step3_st, {
    step('st_step3')
  })

  observeEvent(input$RunStep3_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    output.dir <- output.dir.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    normalization.method <- input$normalization.method
    npcs <- input$npcs
    pcs.used <- input$pcs.used
    resolution <- input$resolution

    # convert some parameters
    pcs.used <- seq(from = as.numeric(unlist(strsplit(pcs.used, ":")))[1], to = as.numeric(unlist(strsplit(pcs.used, ":")))[2])

    shinyjs::runjs('$("#runningStep3").text("Running Step3...");')
    shinyjs::disable('normalization.method')
    shinyjs::disable('npcs')
    shinyjs::disable('pcs.used')
    shinyjs::disable('resolution')
    shinyjs::disable('RunStep3')

    st_obj <- st_Clustering(
        st_obj = st_obj,
        output.dir = file.path(output.dir, 'Step3_Clustering'),
        normalization.method = normalization.method,
        npcs = npcs,
        pcs.used = pcs.used,
        resolution = resolution,
        verbose = FALSE
    )
      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    normalization.method.temp(normalization.method)
    npcs.temp(npcs)
    pcs.used.temp(pcs.used)
    resolution.temp(resolution)

    shinyjs::enable('normalization.method')
    shinyjs::enable('npcs')
    shinyjs::enable('pcs.used')
    shinyjs::enable('resolution')
    shinyjs::enable('RunStep3')
    output$step3_completed <- renderText({'Step3 completed.'})

  })

  # Step4. Differential expressed genes---------------------------------------------------------------------------------
  observeEvent(input$step4_st, {
    step('st_step4')
  })

  observeEvent(input$RunStep4_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    output.dir <- output.dir.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    only.pos <- as.logical(input$only.pos)
    min.pct <- input$min.pct
    logfc.threshold <- input$logfc.threshold
    test.use <- input$test.use

    # convert some parameters
    #

    shinyjs::runjs('$("#runningStep4").text("Running Step4...");')  
    shinyjs::disable('only.pos')
    shinyjs::disable('min.pct')
    shinyjs::disable('logfc.threshold')
    shinyjs::disable('test.use')
    shinyjs::disable('RunStep4')

    st.markers <- st_Find_DEGs(
        st_obj = st_obj,
        output.dir = file.path(output.dir, 'Step4_Find_DEGs'),
        ident.label = 'seurat_clusters',
        only.pos = only.pos,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        test.use = test.use,
        verbose = verbose
    )

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    only.pos.temp(only.pos)
    min.pct.temp(min.pct)
    logfc.threshold.temp(logfc.threshold)
    test.use.temp(test.use)

    shinyjs::enable('only.pos')
    shinyjs::enable('min.pct')
    shinyjs::enable('logfc.threshold')
    shinyjs::enable('test.use')
    shinyjs::enable('RunStep4')
    output$step4_completed <- renderText({'Step4 completed.'})

  })

  # Step5. Spatially variable features---------------------------------------------------------------------------------
  observeEvent(input$step5_st, {
    step('st_step5')
  })

  observeEvent(input$RunStep5_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    output.dir <- output.dir.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    selection.method <- input$selection.method
    n.top.show <- input$n.top.show
    n.col.show <- input$n.col.show

    # convert some parameters
    #

    shinyjs::runjs('$("#runningStep5").text("Running Step5...");')  
    shinyjs::disable('selection.method')
    shinyjs::disable('n.top.show')
    shinyjs::disable('n.col.show')
    shinyjs::disable('RunStep5')

    st_obj <- st_SpatiallyVariableFeatures(
        st_obj = st_obj,
        output.dir = file.path(output.dir, 'Step5_SpatiallyVariableFeatures'),
        assay = st_obj@active.assay,
        selection.method = selection.method,
        n.top.show = n.top.show,
        n.col = n.col.show,
        verbose = FALSE
    )
      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    selection.method.temp(selection.method)
    n.top.show.temp(n.top.show)
    n.col.show.temp(n.col.show)

    shinyjs::enable('selection.method')
    shinyjs::enable('n.top.show')
    shinyjs::enable('n.col.show')
    shinyjs::enable('RunStep5')
    output$step5_completed <- renderText({'Step5 completed.'})

  })

  # Step6. Spatial interaction---------------------------------------------------------------------------------
  observeEvent(input$step6_st, {
    step('st_step6')
  })

  observeEvent(input$RunStep6_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    input.data.dir <- input.data.dir.temp()
    species <- species.temp()
    output.dir <- output.dir.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    commot.signaling_type <- input$commot.signaling_type
    commot.database <- input$commot.database
    commot.min_cell_pct <- input$commot.min_cell_pct
    commot.dis_thr <- input$commot.dis_thr
    commot.n_permutations <- input$commot.n_permutations

    # convert some parameters
    #

    shinyjs::runjs('$("#runningStep6").text("Running Step6...");')  
    shinyjs::disable('commot.signaling_type')
    shinyjs::disable('commot.database')
    shinyjs::disable('commot.min_cell_pct')
    shinyjs::disable('commot.dis_thr')
    shinyjs::disable('commot.n_permutations')
    shinyjs::disable('RunStep6')

    interaction_path = file.path(output.dir, 'Step6_Interaction')
    if(!dir.exists(interaction_path)){
        dir.create(interaction_path)
    }
    write.csv(st_obj@meta.data,
              file.path(interaction_path, 'metadata.csv'),
              row.names = TRUE)
    st_interaction(
        st_data_path = input.data.dir,
        metadata_path = file.path(interaction_path, 'metadata.csv'),
        label_key = 'seurat_clusters',
        save_path = interaction_path,
        species = species,
        signaling_type = commot.signaling_type,
        database = commot.database,
        min_cell_pct = commot.min_cell_pct,
        dis_thr = commot.dis_thr,
        n_permutations = commot.n_permutations
        #condaenv = condaenv
    )

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    commot.signaling_type.temp(commot.signaling_type)
    commot.database.temp(commot.database)
    commot.min_cell_pct.temp(commot.min_cell_pct)
    commot.dis_thr.temp(commot.dis_thr)
    commot.n_permutations.temp(commot.n_permutations)

    shinyjs::enable('commot.signaling_type')
    shinyjs::enable('commot.database')
    shinyjs::enable('commot.min_cell_pct')
    shinyjs::enable('commot.dis_thr')
    shinyjs::enable('commot.n_permutations')
    shinyjs::enable('RunStep6')
    output$step6_completed <- renderText({'Step6 completed.'})

  })

  # Step7. Loading data---------------------------------------------------------------------------------
  observeEvent(input$step7_st, {
    step('st_step7')
  })

  observeEvent(input$RunStep7_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    species <- species.temp()
    output.dir <- output.dir.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    copykat.LOW.DR <- input$copykat.LOW.DR
    copykat.UP.DR <- input$copykat.UP.DR
    copykat.win.size <- input$copykat.win.size
    copykat.distance <- input$copykat.distance
    copykat.genome <- input$copykat.genome
    copykat.n.cores <- input$copykat.n.cores
    # convert some parameters
    if(copykat.genome=='NULL'){copykat.genome <- NULL}else{copykat.genome <- copykat.genome}

    shinyjs::runjs('$("#runningStep7").text("Running Step7...");')  
    shinyjs::disable('copykat.LOW.DR')
    shinyjs::disable('copykat.UP.DR')
    shinyjs::disable('copykat.win.size')
    shinyjs::disable('copykat.distance')
    shinyjs::disable('copykat.genome')
    shinyjs::disable('copykat.n.cores')
    shinyjs::disable('RunStep7')

    st_obj <- st_CNV(
        st_obj = st_obj,
        save_path = file.path(output.dir, 'Step7_CNV_analysis'),
        assay = 'Spatial',
        LOW.DR = copykat.LOW.DR,
        UP.DR = copykat.UP.DR,
        win.size = copykat.win.size,
        distance = copykat.distance,
        genome = copykat.genome,
        n.cores = copykat.n.cores,
        species = species
    )

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    copykat.LOW.DR.temp(copykat.LOW.DR)
    copykat.UP.DR.temp(copykat.UP.DR)
    copykat.win.size.temp(copykat.win.size)
    copykat.distance.temp(copykat.distance)
    copykat.genome.temp(copykat.genome)
    copykat.n.cores.temp(copykat.n.cores)

    shinyjs::enable('copykat.LOW.DR')
    shinyjs::enable('copykat.UP.DR')
    shinyjs::enable('copykat.win.size')
    shinyjs::enable('copykat.distance')
    shinyjs::enable('copykat.genome')
    shinyjs::enable('copykat.n.cores')
    shinyjs::enable('RunStep7')
    output$step7_completed <- renderText({'Step7 completed.'})

  })

  # Step8. Deconvolution---------------------------------------------------------------------------------
  observeEvent(input$step8_st, {
    step('st_step8')
  })

  observeEvent(input$RunStep8_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    input.data.dir <- input.data.dir.temp()
    species <- species.temp()
    output.dir <- output.dir.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    cell2loc.sc.h5ad.dir <- input$cell2loc.sc.h5ad.dir
    cell2loc.sc.max.epoch <- input$cell2loc.sc.max.epoch
    cell2loc.st.max.epoch <- input$cell2loc.st.max.epoch
    cell2loc.use.gpu <- as.logical(input$cell2loc.use.gpu)

    # convert some parameters
    if(cell2loc.sc.h5ad.dir=='NULL'){cell2loc.sc.h5ad.dir <- NULL}else{cell2loc.sc.h5ad.dir <- cell2loc.sc.h5ad.dir}

    shinyjs::runjs('$("#runningStep8").text("Running Step8...");')  
    shinyjs::disable('cell2loc.sc.h5ad.dir')
    shinyjs::disable('cell2loc.sc.max.epoch')
    shinyjs::disable('cell2loc.st.max.epoch')
    shinyjs::disable('cell2loc.use.gpu')
    shinyjs::disable('RunStep8')

    st_obj <- st_Deconvolution(
        st.data.dir = input.data.dir,
        sc.h5ad.dir = cell2loc.sc.h5ad.dir,
        st_obj = st_obj,
        save_path = file.path(output.dir, 'Step8_Deconvolution'),
        sc.labels.key = 'seurat_clusters',
        species = species,
        sc.max.epoch = cell2loc.sc.max.epoch,
        st.max.epoch = cell2loc.st.max.epoch,
        use.gpu = cell2loc.use.gpu
        # condaenv = condaenv
    )

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    cell2loc.sc.h5ad.dir.temp(cell2loc.sc.h5ad.dir)
    cell2loc.sc.max.epoch.temp(cell2loc.sc.max.epoch)
    cell2loc.st.max.epoch.temp(cell2loc.st.max.epoch)
    cell2loc.use.gpu.temp(cell2loc.use.gpu)

    shinyjs::enable('cell2loc.sc.h5ad.dir')
    shinyjs::enable('cell2loc.sc.max.epoch')
    shinyjs::enable('cell2loc.st.max.epoch')
    shinyjs::enable('cell2loc.use.gpu')
    shinyjs::enable('RunStep8')
    output$step8_completed <- renderText({'Step8 completed.'})

  })

  # Step9. Cellcycle---------------------------------------------------------------------------------
  observeEvent(input$step9_st, {
    step('st_step9')
  })

  observeEvent(input$RunStep9_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    output.dir <- output.dir.temp()
    species <- species.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    s.features <- input$s.features
    g2m.features <- input$g2m.features
    FeatureColors.bi <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu')))

    # convert some parameters
    s.features <- unlist(strsplit(s.features, ","))
    g2m.features <- unlist(strsplit(g2m.features, ","))

    shinyjs::runjs('$("#runningStep9").text("Running Step9...");')  
    shinyjs::disable('s.features')
    shinyjs::disable('g2m.features')
    shinyjs::disable('RunStep9')

    st_obj <- st_Cell_cycle(
        st_obj = st_obj,
        save_path = file.path(output.dir, 'Step9_Cellcycle'),
        s.features = s.features,
        g2m.features = g2m.features,
        species = species,
        FeatureColors.bi = FeatureColors.bi
    )

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    s.features.temp(s.features)
    g2m.features.temp(g2m.features)

    shinyjs::enable('s.features')
    shinyjs::enable('g2m.features')
    shinyjs::enable('RunStep9')
    output$step9_completed <- renderText({'Step9 completed.'})

  })

  # Step10. Niche analysis---------------------------------------------------------------------------------
  observeEvent(input$step10_st, {
    step('st_step10')
  })

  observeEvent(input$RunStep10_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    output.dir <- output.dir.temp()
    slice <- slice.temp()
    species <- species.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    Nich.cluster.n <- input$Nich.cluster.n

    # convert some parameters
    #

    shinyjs::runjs('$("#runningStep10").text("Running Step10...");')  
    shinyjs::disable('Nich.cluster.n')
    shinyjs::disable('RunStep10')

    if(!file.exists(file.path(output.dir, 'Step8_Deconvolution', 'cell2loc_res.csv'))){
        stop('Please run step 8  deconvolution first.')
    }
    tmp <- read.csv(file.path(output.dir, 'Step8_Deconvolution', 'cell2loc_res.csv'),
                    row.names = 1)
    features <- colnames(tmp)

    if(!all(features %in% names(st_obj@meta.data))){
        common.barcodes <- intersect(colnames(st_obj), rownames(tmp))
        tmp <- tmp[common.barcodes, ]
        st_obj <- st_obj[, common.barcodes]
        st_obj <- AddMetaData(st_obj,
                              metadata = tmp)
    }

    rm(tmp)
    gc()

    st_obj <- st_NicheAnalysis(
        st_obj,
        features = features,
        save_path = file.path(output.dir, 'Step10_NicheAnalysis'),
        kmeans.n = Nich.cluster.n,
        st_data_path = file.path(output.dir, 'Step1_Loading_Data'),
        slice = slice,
        species = species,
        condaenv = condaenv
    )

      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    Nich.cluster.n.temp(Nich.cluster.n)

    shinyjs::enable('Nich.cluster.n')
    shinyjs::enable('RunStep10')
    output$step10_completed <- renderText({'Step10 completed.'})
  })

  # Step11. Generate the Report---------------------------------------------------------------------------------
  observeEvent(input$step11_st, {
    step('st_step11')
  })

  observeEvent(input$RunStep11_st, {
    # load previous variables
    previous_results_path <- previous_results_path.temp()

    # load previous parameters
    output.dir <- output.dir.temp()

    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)

    # set parameters
    #

    # convert some parameters
    #

    shinyjs::runjs('$("#runningStep11").text("Running Step11...");')  
    shinyjs::disable('RunStep11')

    #### Save data ####
    output.dir.final <- file.path(output.dir, 'Data_and_report')
    if(!dir.exists(output.dir.final)){
        dir.create(output.dir.final)
    }
    saveRDS(st_obj, file.path(output.dir.final, 'st_object.rds'))
    suppressMessages(suppressWarnings(
        SaveH5Seurat(st_obj,
                     filename = file.path(output.dir.final, 'st_object.h5Seurat'),
                     overwrite = TRUE)
    ))
    suppressMessages(suppressWarnings(
        Convert(file.path(output.dir.final, 'st_object.h5Seurat'), dest = "h5ad",
                overwrite = TRUE)
    ))

    #### Generate the report ####
    print('Generating the report...')
    if(genReport){
        knitr::knit(file.path(system.file(package = "HemaScopeR"), "rmd/st_base.Rmd"),
                    file.path(output.dir.final, 'st_pipeline.md'))
        markdown::markdownToHTML(file.path(output.dir.final, 'st_pipeline.md'),
                                 file.path(output.dir.final, 'st_pipeline.html'))
    }
      # Get the names of all variables in the current environment
      variable_names <- ls()
      # Loop through the variable names and save them as RDS files
      for (var_name in variable_names) {
        var <- get(var_name)  # Get the variable by its name
        saveRDS(var, file = paste0(output.dir, '/RDSfiles/', var_name, ".rds"))  # Save as RDS with the variable's name
      }

    # update parameters
    #

    shinyjs::enable('RunStep11')
    output$step10_completed <- renderText({'Step11 completed.'})
  })

  observeEvent(input$back_button_st, {
    shinyjs::hide("ui3")
    shinyjs::show("ui1")
  })

}
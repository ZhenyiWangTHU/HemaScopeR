library(shiny) 
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
library(slickR)
# sc libraries
library(Seurat)
library(phateR)
library(DoubletFinder)
library(monocle)
library(slingshot)
#library(URD) #没有这个包
library(GSVA)
library(limma)
library(plyr)
library(dplyr)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(CellChat)
library(velocyto.R)
library(SeuratWrappers)
library(stringr)
library(scran)
library(ggpubr)
library(viridis)
library(pheatmap)
library(parallel)
library(reticulate)
library(SCENIC)
library(feather)
library(AUCell)
library(RcisTarget)
library(Matrix)
library(foreach)
library(doParallel)
library(clusterProfiler)
library(tools)
#st libraries
library(RColorBrewer)
library(Rfast2) 
library(SeuratDisk)
library(abcCellmap)
library(biomaRt)
library(copykat)
library(gelnet)
library(ggplot2)
library(parallelDist)
library(patchwork)
library(markdown)
# getpot
#library(getopt) #没有这个包
library(tools)
# HemaScopeR
library(HemaScopeR) 

Load_previous_results = function(previous_results_path=NULL){
  # Get a list of all .RDS files in the specified path
  rds_files_list <- list.files(previous_results_path, pattern = "\\.rds$", full.names = TRUE)
  # Loop through the list of RDS files and load them into variables
  for (rds_file in rds_files_list) {
    # Extract the variable name from the file name
    var_name <- tools::file_path_sans_ext(basename(rds_file))
    if(!(var_name %in% c(   # input and output
      #'input.data.dirs',
      #'project.names', 
      #'output.dir',
      #'pythonPath',
      # quality control and preprocessing
      'gene.column',
      'min.cells',
      'min.feature',
      'mt.pattern',
      'nFeature_RNA.limit',
      'percent.mt.limit',
      'scale.factor',
      'nfeatures',
      'ndims',
      'vars.to.regress',
      #'PCs',
      'resolution',
      'n.neighbors',
      # remove doublets
      'doublet.percentage',
      'doublerFinderwraper.PCs',
      'doublerFinderwraper.pN',
      'doublerFinderwraper.pK',
      # phateR
      'phate.knn',
      'phate.npca',
      'phate.t',
      'phate.ndim',
      'min.pct',
      'logfc.threshold',
      # visualization
      'marker.genes',
      #'ViolinPlot.cellTypeOrders',
      #'ViolinPlot.cellTypeColors',
      #'Org',
      'lineage.genelist',
      'lineage.names',
      'groups_colors',
      #slingshot
      'slingshot.start.clus',
      'slingshot.end.clus',
      'slingshot.colors',
      'loom.files.path',
      # cell cycle
      'cellcycleCutoff',
      # cell chat
      'sorting',
      'ncores',
      # Verbose = FALSE,
      # activeEachStep
      'Whether_load_previous_results',
      'Step1_Input_Data',
      'Step1_Input_Data.type',
      'Step2_Quality_Control',
      #'Step2_Quality_Control.RemoveBatches',
      #'Step2_Quality_Control.RemoveDoublets',
      'Step3_Clustering',
      'Step4_Identify_Cell_Types',
      'Step4_Use_Which_Labels',
      'Step4_Cluster_Labels',
      'Step4_Changed_Labels',
      'Step4_run_sc_CNV',
      'Step5_Visualization',
      'Step6_Find_DEGs',
      'Step7_Assign_Cell_Cycle',
      'Step8_Calculate_Heterogeneity',
      'Step9_Violin_Plot_for_Marker_Genes',
      'Step10_Calculate_Lineage_Scores',
      'Step11_GSVA',
      'Step11_GSVA.identify.cellType.features',
      'Step11_GSVA.identify.diff.features',
      'Step11_GSVA.comparison.design',
      'Step12_Construct_Trajectories',
      'Step12_Construct_Trajectories.monocle',
      'Step12_Construct_Trajectories.slingshot',
      'Step12_Construct_Trajectories.scVelo',
      'Step13_TF_Analysis',
      'Step13_TF_Analysis.cellTypes_colors',
      'Step13_TF_Analysis.groups_colors',
      'Step14_Cell_Cell_Interaction',
      'Step15_Generate_the_Report'
    ))){  
      # Load the RDS file and assign it to the variable with the same name 
      rds_file.temp <- base::readRDS(rds_file)
      assign(var_name, rds_file.temp, envir = .GlobalEnv)
    }
  }
}
# ui--------------------------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(
  shinyjs::useShinyjs(),  
  # ui1
  div(id = "ui1", style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 70vh;",
      fluidRow(),
      fluidRow(
        column(3, align = "center", imageOutput('logo'))
      ),
      fluidRow(
        column(12, align = "center", h1("HemaScopeR: A Specialized Bioinformatics Toolkit Designed for Analyzing both Single-cell and Spatial Transcriptome Sequencing Data from Hematopoietic Cells", 
                                        class = "h1-font",style = "font-family: 'arial'; font-size: 28pt;font-weight: bold;"))
      ),
      fluidRow(div(class = "spacer")),  # empty line
      fluidRow(div(class = "spacer")),  # empty line
      div(style = "display: flex; justify-content: center; width: 100%;",
          actionButton("start_button", "Start scRNA-seq pipeline",style = "font-family: 'arial'; font-size: 18pt;background-color: #dae3f5;padding: 20px 20px; margin-right: 20px;"),
          actionButton("start_button_st", "Start ST-seq pipeline",style = "font-family: 'arial'; font-size: 18pt;background-color: #dae3f5;padding: 20px 20px; margin-right: 20px;")
      ),
      uiOutput("ui_styles")
  ),
  #ui2.1 ,zyt add this code
  div(id = "ui2.1",style = "display: none;",#设置为隐藏不显示
      fluidRow(
        box(width = 4, status = "primary", solidHeader = TRUE,
            actionButton("new_analysize_btn", "Begin New Analysis", 
                         style = "width:100%; height:100px; font-size:20px;background-color: white; color: black;")),
        box(width = 4, status = "warning", solidHeader = TRUE,
            actionButton("continue_analysize_btn", "Continue Previous Analysis", 
                         style = "width:100%; height:100px; font-size:20px;background-color: white; color: black;")),
        box(width = 4, status = "danger", solidHeader = TRUE,
            actionButton("sc_return_home", "Back to Home", 
                         style = "width:100%; height:100px; font-size:20px;background-color: white; color: black;"))
      )
  ),
  
  
  div(id = "ui2.2",style = "display: none;",
      uiOutput("dynamic_ui"),
      uiOutput("continue_stepContent")
  ),
  
  # ui2
  div(id = "ui2", style = "display: none;",
      dashboardPage(
        skin = c("blue"),
        dashboardHeader(title = "scRNA-seq pipeline"),
        dashboardSidebar(
          disable = FALSE,
          sidebarMenu(
            id="sc_step_tab",selected=TRUE,
            menuItem("Step 1. Input Data", tabName = "step_1"),
            menuItem("Step 2. Quality Control", tabName = "step_2"),
            menuItem("Step 3. Clustering", tabName = "step_3"),
            menuItem("Step 4. Identify Cell Types", tabName = "step_4"),
            menuItem("Step 5. Visualization", tabName = "step_5"),
            menuItem("Step 6. Find DEGs", tabName = "step_6"),
            #menuItem("Step 6. Find DEGs", tabName = "Step 6"),
            menuItem("Step 7. Assign Cell Cycles", tabName = "step_7"),
            menuItem("Step 8. Calculate Heterogeneity", tabName = "step_8"),
            menuItem("Step 9. Violin Plot for Marker Genes", tabName = "step_9"),
            menuItem("Step 10. Calculate Lineage Scores", tabName = "step_10"),
            menuItem("Step 11. GSVA", tabName = "step_11"),
            menuItem("Step 12. Construct Trajectories", tabName = "step_12"),
            menuItem("Step 13. TF Analysis" , tabName = "step_13"),
            menuItem("Step 14. Cell-Cell Interaction", tabName = "step_14"),
            menuItem("Step 15. Generate the Report" , tabName = "step_15"),
            menuItem("Back to Prior Page",tabName = "back_button")
          )
        ),
        dashboardBody(
          uiOutput("stepContent")
        )
      )
  ),
  
  # ui3
  div(id = "ui3", style = "display: none;",
      dashboardPage(
        skin = c("blue"),
        dashboardHeader(title = "ST-seq pipeline"),
        dashboardSidebar(
          disable = FALSE,
          sidebarMenu(
            id="st_step_tab",selected=TRUE,
            menuItem("Step 1. Input Data", tabName = "step1_st"),
            menuItem("Step 2. Quality Control", tabName = "step2_st"),
            menuItem("Step 3. Clustering", tabName = "step3_st"),
            menuItem("Step 4. Find Differential Genes", tabName = "step4_st"),
            menuItem("Step 5. Spatially Variable Features", tabName = "step5_st"),
            menuItem("Step 6. Spatial Interaction", tabName = "step6_st"),
            menuItem("Step 7. CNV Analysis", tabName = "step7_st"),
            menuItem("Step 8. Deconvolution", tabName = "step8_st"),
            menuItem("Step 9. Cell Cycle Analysis", tabName = "step9_st"),
            menuItem("Step 10. Niche Analysis", tabName = "step10_st"),
            menuItem("Step 11. Generate the Report", tabName = "step11_st"),
            menuItem("Back to Prior Page",tabName = "back_button_st")
          )
        ),
        dashboardBody(
          uiOutput("stepContent_st")
        )
      )
  ),
  div(id = "ui3.1",style = "display: none;",
      fluidRow(
        box(width = 4, status = "primary", solidHeader = TRUE,
            actionButton("st_new_analysize_btn", "Begin New Analysis", 
                         style = "width:100%; height:100px; font-size:20px;background-color: white; color: black;")),
        box(width = 4, status = "warning", solidHeader = TRUE,
            actionButton("st_continue_analysize_btn", "Continue Previous Analysis", 
                         style = "width:100%; height:100px; font-size:20px;background-color: white; color: black;")),
        box(width = 4, status = "danger", solidHeader = TRUE,
            actionButton("st_return_home", "Back to Home", 
                         style = "width:100%; height:100px; font-size:20px;background-color: white; color: black;"))
      )
  ),
  
  div(id = "ui3.2",style = "display: none;",
      uiOutput("dynamic_st_ui"),
      uiOutput("st_continue_step")
  )
)
#scRNA UI step1-15
step1_fluidRow<-fluidRow(
  style = "margin-left: 10px;",
  column(
    6, align = "left", h3("Step 1. Input Data"),
    #p("Please select the input data."),
    textInput("input.data.dirs", 
              HTML("Enter data path and use ';' to separate files<br>(e.g. /path1/file1/data1;/path2/file2/data2)"), 
              value = "NULL"),
    textInput("project.names", "Enter project name:", value = "project.names"),
    textInput("output.dir", "Enter output path:", value = ""),
    textInput("pythonPath", "Enter the path of Python:", value = "NULL"),
    textInput("databasePath", "Enter the path of database:", value = "NULL"),
    pickerInput("Step1_Input_Data.type", "Select Data Type:", choices = c("cellranger-count", "Seurat", "Matrix")),
    numericInput("gene.column", "Gene Column (default: 2):", value = 2),
    numericInput("min.cells", "Minimum Cells (default: 10):", value = 10),
    numericInput("min.feature", "Minimum Features (default: 200):", value = 200),
    textInput("mt.pattern", "Mt Pattern (default: '^MT-'):", value = '^MT-'),
    actionBttn("load_data_button", "Load Data", style = "unite", color = "primary"),
    div(class = "spacer"),
    uiOutput("loadingData"),
    div(class = "spacer"),
    uiOutput("data_dim_output")
  ),
  column(3,align="middle",
         h3("Please record your job ID for use in future analyses"),
         #h5("Please copy your jobid and you will use it in subsequent analysis"),
         textOutput("jobid")
  )
)

step2_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 2. Quality Control"),
    #p("Please input the parameters for quality control and preprocessing."),
    #textInput("jobid","Please Paste the token generated in step 1"), #zyt add this code
    numericInput("nFeature_RNA.limit", "nFeature_RNA.limit (default: 200):", value = 200),
    numericInput("percent.mt.limit", "percent.mt.limit (default: 20):", value = 20),
    numericInput("scale.factor", "scale factor (default: 10,000):", value = 10000),
    numericInput("nfeatures", "nfeatures (default: 3,000):", value = 3000),
    numericInput("ndims", "ndims (default: 50):", value = 50),
    textInput("vars.to.regress", "vars.to.regress (default: NULL):", value = 'NULL'),
    textInput("PCs", "PCs (default: 1:35):", value = "1:35"),
    numericInput("resolution", "resolution (default: 0.4):", value = 0.4),
    numericInput("n.neighbors", "n.neighbors (default: 50):", value = 50),
    numericInput("doublet.percentage", "doublet percentage (default: 0.04):", value = 0.04),
    textInput("doublerFinderwraper.PCs", "doubletFinder Wrapper PCs (default: 1:20):", value = "1:20"),
    numericInput("doublerFinderwraper.pN", "doubletFinder Wrapper pN (default: 0.25):", value = 0.25),
    numericInput("doublerFinderwraper.pK", "doubletFinder Wrapper pK (default: 0.1):", value = 0.1),
    selectInput("Step2_Quality_Control.RemoveBatches", "Step2_Quality_Control.RemoveBatches:", choices = c("FALSE" = FALSE, "TRUE" = TRUE)),
    selectInput("Step2_Quality_Control.RemoveDoublets", "Step2_Quality_Control.RemoveDoublets:", choices = c("FALSE" = FALSE, "TRUE" = TRUE)),
    actionBttn("RunStep2", "Run Step2",style = "unite",color = "primary"),
    #actionButton("RunStep2", "Run Step2"),
    div(class = "spacer"), 
    uiOutput("runningStep2"),
    div(class = "spacer"),  
    uiOutput("step2_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 2. Quality Control:"),
    uiOutput("Step2.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step2_plot", width='auto',height = "500px"),
         textOutput("step2_text"))
)


step3_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 3. Clustering"),
    #p("Please input the parameters for clustering."),
    textInput("PCs.clustering", "PCs for clustering (default: 1:20):", value = "1:20"),
    numericInput("n.neighbors", "n.neighbors for clustering (default: 50):", value = 50),
    numericInput("resolution", "resolution for clustering (default: 0.4):", value = 0.4),
    actionBttn("RunStep3", "Run Step3",style = "unite",color = "primary"),
    #actionButton("RunStep3", "Run Step3"),
    div(class = "spacer"),
    uiOutput("runningStep3"),
    div(class = "spacer"),
    uiOutput("step3_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 3. Clustering:"),
    uiOutput("Step3.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step3_plot", width='auto',height = "500px"),
         textOutput("step3_text"))
)

step4_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 4. Identify Cell Types"),
    #p("Please input the parameters for identifying cell types automatically."),
    selectInput("Org", "Choose organism (e.g., 'mmu' for mouse, 'hsa' for human):", choices = c("mmu" = "mmu", "hsa" = "hsa")),
    selectInput("Step4_Use_Which_Labels", "Choose Labels:", choices = c('clustering' = 'clustering',
                                                                        'abcCellmap.1' = 'abcCellmap.1',
                                                                        'abcCellmap.2' = 'abcCellmap.2',
                                                                        'abcCellmap.3' = 'abcCellmap.3',
                                                                        'abcCellmap.4' = 'abcCellmap.4',
                                                                        'HematoMap' = 'HematoMap')),
    selectInput("Step4_run_sc_CNV", "Run CNV:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
    numericInput("ncores", "CPU cores for parallel processing (default: 1):", value = 1),
    actionBttn("RunStep4", "Run Step4",style = "unite",color = "primary"),
    #actionButton("RunStep4", "Run Step4"),
    div(class = "spacer"),
    uiOutput("runningStep4"),
    div(class = "spacer"),
    uiOutput("step4_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 4. Identify Cell Types:"),
    uiOutput("Step4.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step4_plot", width='auto',height = "500px"),
         textOutput("step4_text"))
)

step5_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 5. Visualization"),
    #p("Please input the parameters for visualization."),
    numericInput("phate.knn", "Nearest neighbors for PhateR analysis (default: 50):", value = 50),
    numericInput("phate.npca", "Principal components for PhateR (default: 20):", value = 20),
    numericInput("phate.t", "t parameter for PhateR (default: 10):", value = 10),
    numericInput("phate.ndim", "Dimensions for PhateR (default: 2):", value = 2),  
    actionBttn("RunStep5", "Run Step5",style = "unite",color = "primary"),
    #actionButton("RunStep5", "Run Step5"),
    div(class = "spacer"), 
    uiOutput("runningStep5"),
    div(class = "spacer"),  
    uiOutput("step5_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 5. Visualization:"),
    uiOutput("Step5.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step5_plot", width='auto',height = "500px"),
         textOutput("step5_text"))
)

step6_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 6. Find Differential Genes"),
    #p("Please input the parameters for finding DEGs."),
    
    numericInput("min.pct", "Minimum gene percentage for differential detection (default: 0.25):", value = 0.25),
    numericInput("logfc.threshold", "Log-fold threshold for gene analysis(Default: 0.25):", value = 0.25),  
    actionBttn("RunStep6", "Run Step6",style = "unite",color = "primary"),
    #actionButton("RunStep6", "Run Step6"),
    div(class = "spacer"), 
    uiOutput("runningStep6"),
    div(class = "spacer"),  
    uiOutput("step6_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 6. Find Differential Genes:"),
    uiOutput("Step6.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step6_plot", width='auto',height = "500px"),
         textOutput("step6_text"))
)

step7_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 7. Assign Cell Cycles"),
    #p("Please input the parameters for assigning cell cycles."),
    textInput("cellcycleCutoff", "Define cell cycle cutoff (default: NULL):", value = "NULL"),
    actionBttn("RunStep7", "Run Step7",style = "unite",color = "primary"),
    #actionButton("RunStep7", "Run Step7"),
    div(class = "spacer"), 
    uiOutput("runningStep7"),
    div(class = "spacer"),  
    uiOutput("step7_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 7. Assign Cell Cycles:"),
    uiOutput("Step7.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step7_plot", width='auto',height = "500px"),
         textOutput("step7_text"))
)

step8_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 8. Calculate Heterogeneity"),
    #p("Please set the parameter needs for calculating heterogeneity."),
    textInput("ViolinPlot.cellTypeOrders", "Order cell types (seperate by ','):", value = "NULL"),
    actionBttn("RunStep8", "Run Step8",style = "unite",color = "primary"),
    #actionButton("RunStep8", "Run Step8"),
    div(class = "spacer"), 
    uiOutput("runningStep8"),
    div(class = "spacer"),  
    uiOutput("step8_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 8. Calculate Heterogeneity:"),
    uiOutput("Step8.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step8_plot", width='auto',height = "500px"),
         textOutput("step8_text"))
)

step9_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 9. Violin Plot for Marker Genes"),
    #p("Set parameters for violin plot of marker genes."),
    textInput("marker.genes", "Enter marker genes for violin plot (seperate by ','):", value = "NULL"),
    textInput("ViolinPlot.cellTypeColors", "Set the hexadecimal codes of colors for cell types (seperate by ','):", value = "NULL"),
    actionBttn("RunStep9", "Run Step9",style = "unite",color = "primary"),
    #actionButton("RunStep9", "Run Step9"),
    div(class = "spacer"), 
    uiOutput("runningStep9"),
    div(class = "spacer"),  
    uiOutput("step9_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 9. Violin Plot for Marker Genes:"),
    uiOutput("Step9.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step9_plot", width='auto',height = "500px"),
         textOutput("step9_text"))
)

step10_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 10. Calculate Lineage Scores"),
    #p("Please input the parameters for calculating lineage scores."),
    textInput("lineage.genelist", HTML('The gene sets for calculating lineage scores.<br>Please enclose the gene sets with double quotation marks and seperate them by ";". <br>(e.g. "gene1,gene2,gene3";"gene4,gene5,gene6"):'), value = "NULL"),
    textInput("lineage.names", HTML('The names for the lineages. Please use "," to seperate them.<br>(e.g. lineage1,lineage2):'), value = "NULL"),
    textInput("groups_colors", HTML('The hexadecimal codes of colors for groups. Please use "," to seperate them.<br>(e.g. #FF0000,#0000FF):'), value = "NULL"),
    actionBttn("RunStep10", "Run Step10",style = "unite",color = "primary"),
    #actionButton("RunStep10", "Run Step10"),
    div(class = "spacer"), 
    uiOutput("runningStep10"),
    div(class = "spacer"),  
    uiOutput("step10_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 10. Calculate Lineage Scores:"),
    uiOutput("Step10.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step10_plot", width='auto',height = "500px"),
         textOutput("step10_text"))
)

step11_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 11. GSVA"),
    #p("Please input the parameters for GSVA."),
    selectInput("Step11_GSVA.identify.cellType.features", "Option to identify cell type-specific GSVA terms:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
    selectInput("Step11_GSVA.identify.diff.features", "Option to identify differential GSVA terms:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
    textInput("Step11_GSVA.comparison.design", 
              HTML('The comparison design for GSVA<br>
                       e.g. list(
                            list(c("sample1","sample2"),
                                 c("sample3","sample4")),<br>
                            list(c("sample5","sample6"),
                                 c("sample7","sample8")))<br>
                            which indicates "sample1 and sample2 vs sample 3 and sample4",<br>
                             and "sample5 and sample6 vs sample 7 and sample8" :'), 
              value = "NULL"),
    actionBttn("RunStep11", "Run Step11",style = "unite",color = "primary"),
    #actionButton("RunStep11", "Run Step11"),
    div(class = "spacer"), 
    uiOutput("runningStep11"),
    div(class = "spacer"),  
    uiOutput("step11_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 11. GSVA:"),
    uiOutput("Step11.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step11_plot", width='auto',height = "500px"),
         textOutput("step11_text"))
)

step12_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 12. Construct Trajectories"),
    textInput("Step12_Construct_Trajectories.clusters",
              "Set the cell types for constructing trajectories(seperate by ',' default is 'all'):", value = "all"),
    # monocle
    h4("Set the parameters for monocle2",style = "font-weight: bold; margin-top: 15px; margin-bottom: 15px;"),
    selectInput("Step12_Construct_Trajectories.monocle", "Option to run monocle2:", 
                choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
    # slingshot
    h4("Set the parameters for slingshot",style = "font-weight: bold; margin-top: 15px; margin-bottom: 15px;"),
    selectInput("Step12_Construct_Trajectories.slingshot", "Option to run slingshot:", 
                choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
    textInput("slingshot.start.clus", "Set the root clusters:", value = "NULL"),
    textInput("slingshot.end.clus", "Set the tip clusters:", value = "NULL"),
    textInput("slingshot.colors", "Set the hexadecimal codes of colors for clusters (seperate by ','):", 
              value = "NULL"),
    # scVelo
    h4("Set the parameters for scVelo",style = "font-weight: bold; margin-top: 15px; margin-bottom: 15px;"),
    selectInput("Step12_Construct_Trajectories.scVelo", "Option to run scVelo:", 
                choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
    textInput("loom.files.path", "Enter the paths of loom files (seperate by ';'):", value = "NULL"),
    actionBttn("RunStep12", "Run Step12",style = "unite",color = "primary"),
    #actionButton("RunStep12", "Run Step12"),
    div(class = "spacer"),
    uiOutput("runningStep12"),
    div(class = "spacer"),
    uiOutput("Step12_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 12. Construct Trajectories:"),
    uiOutput("Step12.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step12_plot", width='auto',height = "500px"),
         textOutput("step12_text"))
)

step13_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 13. Transcription Factors Analysis"),
    #p("Please input the parameters for TF analysis."),
    textInput("Step13_TF_Analysis.cellTypes_colors", "Set the hexadecimal codes of colors for cell types (seperate by ','):", value = "NULL"),
    textInput("Step13_TF_Analysis.groups_colors", "Set the hexadecimal codes of colors for groups (seperate by ','):", value = "NULL"),
    actionBttn("RunStep13", "Run Step13",style = "unite",color = "primary"),
    #actionButton("RunStep13", "Run Step13"),
    div(class = "spacer"), 
    uiOutput("runningStep13"),
    div(class = "spacer"),  
    uiOutput("step13_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 13. Transcription Factors Analysis:"),
    uiOutput("Step13.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step13_plot", width='auto',height = "500px"),
         textOutput("step13_text"))
)

step14_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 14. Cell-Cell Interaction"),
    #p("Please input the parameters for cell-cell interaction analysis."),
    selectInput("sorting", "The cell groups were sorted?:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
    actionBttn("RunStep14", "Run Step14",style = "unite",color = "primary"),
    #actionButton("RunStep14", "Run Step14"),
    div(class = "spacer"), 
    uiOutput("runningStep14"),
    div(class = "spacer"),  
    uiOutput("step14_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 14. Cell-Cell Interection:"),
    uiOutput("Step14.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("step14_plot", width='auto',height = "500px"),
         textOutput("step14_text"))
)

step15_fluidRow <- fluidRow(
  style = "margin-left: 10px;",
  column(
    6, align = "left", h3("Step 15. Generate the Report"),
    #p("Please input the parameters for generating the report."),
    actionBttn("RunStep15", "Run Step15",style = "unite",color = "primary"),
    #actionButton("RunStep15", "Run Step15"),
    div(class = "spacer"), 
    uiOutput("runningStep15"),
    div(class = "spacer"),  
    uiOutput("step15_completed"))
)

step_sc_continue_fluidRow<-fluidRow(
  style = "margin-left: 10px;",
  column(8,align = "left",
         h3("Continue Previous Analysis"),
         actionButton("sc_continue_back_btn","Back to Prior Page",style = "width: 25%;",
                      class = "btn-primary btn-lg"),
         numericInput("jobid","Enter your Job ID",value = NULL),#要是数字输入不然会报错
         selectInput("continue_step","Choose a step you want analysize",
                     choices = c(
                       #"Step 1. Input Data",
                       "Step 2. Quality Control" = "step 2",
                       "Step 3. Clustering" = "step 3",
                       "Step 4. Identify Cell Types" = "step 4",
                       "Step 5. Visualization" = "step 5",
                       "Step 6. Find DEGs" = "step 6",
                       "Step 7. Assign Cell Cycles" = "step 7",
                       "Step 8. Calculate Heterogeneity" = "step 8",
                       "Step 9. Violin Plot for Marker Genes" = "step 9",
                       "Step 10. Calculate Lineage Scores" = "step 10",
                       "Step 11. GSVA" = "step 11",
                       "Step 12. Construct Trajectories" = "step 12",
                       "Step 13. TF Analysis" = "step 13",
                       "Step 14. Cell-Cell Interaction" = "step 14",
                       "Step 15. Generate the Report" = "step 15"
                     ),selected = "step 2")
  ))

#ST UI step1-11
step1_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    6, align = "left", h3("Step 1. Input Data"),
    #p("Please select the input data."),
    textInput("input.data.dir", HTML("Enter data path:"), value = "NULL"),
    textInput("sampleName", "Enter sample name:", value = "Hema_ST"),
    textInput("output.dir", "Enter output path:", value = "NULL"),
    textInput("pythonPath", "Enter the path of Python:", value = "NULL"),
    actionBttn("RunStep1_st", "Load Data",style = "unite",color = "primary"),
    #actionButton("load_data_button", "Load Data"),
    div(class = "spacer"),
    uiOutput("loadingData"),
    uiOutput("step1_completed")),
  column(3,align="middle",
         h3("Please record your job ID for use in future analyses"),
         #p("Please copy your jobid and you will use it in subsequent analysis"),
         textOutput("st_jobid_1"))
)

step2_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 2. Quality Control"),
    #p("Please input the parameters for quality control and preprocessing."),
    numericInput("min.gene", "min.gene (default: 200):", value = 200),
    numericInput("min.nUMI", "min.nUMI (default: 500):", value = 500),
    numericInput("max.gene", "max.gene (default: Inf):", value = Inf),
    numericInput("max.nUMI", "max.nUMI (default: Inf):", value = Inf),
    numericInput("min.spot", "min.spot (default: 0):", value = 0),
    selectInput("bool.remove.mito", "bool.remove.mito:", choices = c("FALSE" = FALSE, "TRUE" = TRUE)),
    selectInput("species", "species:", choices = c("mouse" = "mouse", "human" = "human")),
    actionBttn("RunStep2_st", "Run Step2",style = "unite",color = "primary"),
    #actionButton("RunStep2", "Run Step2"),
    div(class = "spacer"), 
    uiOutput("runningStep2"),
    div(class = "spacer"),  
    uiOutput("step2_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 2. Quality Control:"),
    uiOutput("Step2.st.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("st_step2_plot", width='auto',height = "500px"),
         textOutput("st_step2_text"))
)


step3_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 3. Clustering"),
    #p("Please input the parameters for normalization, PCA and clustering."),
    textInput("normalization.method", "normalization.method (default: 'SCTransform'):", value = "SCTransform"),
    numericInput("npcs", "npcs (default: 50):", value = 50),
    textInput("pcs.used", "pcs.used (default: 1:10):", value = "1:10"),
    numericInput("resolution", "resolution (default: 0.8):", value = 0.8),
    actionBttn("RunStep3_st", "Run Step3",style = "unite",color = "primary"),
    #actionButton("RunStep3", "Run Step3"),
    div(class = "spacer"), 
    uiOutput("runningStep3"),
    div(class = "spacer"),  
    uiOutput("step3_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 3. Clustering:"),
    uiOutput("Step3.st.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("st_step3_plot", width='auto',height = "500px"),
         textOutput("st_step3_text"))
)


step4_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 4. Find Differential Genes"),
    #p("Please input the parameters for finding differentially expressed genes in each cluster."),
    selectInput("only.pos", "only.pos:", choices = c("FALSE" = FALSE, "TRUE" = TRUE)),
    numericInput("min.pct", "min.pct (default: 0.25):", value = 0.25),
    numericInput("logfc.threshold", "logfc.threshold (default: 0.25):", value = 0.25),
    textInput("test.use", "test.use (default: 'wilcox'):", value = 'wilcox'),
    actionBttn("RunStep4_st", "Run Step4",style = "unite",color = "primary"),
    #actionButton("RunStep4", "Run Step4"),
    div(class = "spacer"), 
    uiOutput("runningStep4"),
    div(class = "spacer"),  
    uiOutput("step4_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 4. Find Differential Genes:"),
    uiOutput("Step4.st.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("st_step4_plot", width='auto',height = "500px"),
         textOutput("st_step4_text"))
)

step5_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 5. Spatially variable features"),
    #p("Please input the parameters for finding spatially variable features."),
    textInput("selection.method", "selection.method (default: 'moransi'):", value = 'moransi'),
    numericInput("n.top.show", "n.top.show (default: 10):", value = 10),
    numericInput("n.col.show", "n.col.show (default: 5):", value = 5),
    actionBttn("RunStep5_st", "Run Step5",style = "unite",color = "primary"),
    #actionButton("RunStep5", "Run Step5"),
    div(class = "spacer"), 
    uiOutput("runningStep5"),
    div(class = "spacer"),  
    uiOutput("step5_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 5. Spatially variable features:"),
    uiOutput("Step5.st.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("st_step5_plot", width='auto',height = "500px"),
         textOutput("st_step5_text"))
)

step6_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 6. Spatial interaction"),
    #p("Please input the parameters for analyzing spatial interaction."),
    textInput("commot.signaling_type", "commot.signaling_type (default: 'Secreted Signaling'):", value = 'Secreted Signaling'),
    textInput("commot.database", "commot.database (default: 'CellChat'):", value = 'CellChat'),
    numericInput("commot.min_cell_pct", "commot.min_cell_pct (default: 0.05):", value = 0.05),
    numericInput("commot.dis_thr", "commot.dis_thr (default: 500):", value = 500),
    numericInput("commot.n_permutations", "commot.n_permutations (default: 100):", value = 100),
    actionBttn("RunStep6_st", "Run Step6",style = "unite",color = "primary"),
    #actionButton("RunStep6", "Run Step6"),
    div(class = "spacer"), 
    uiOutput("runningStep6"),
    div(class = "spacer"),  
    uiOutput("step6_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 6. Spatial interaction:"),
    uiOutput("Step6.st.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("st_step6_plot", width='auto',height = "500px"),
         textOutput("st_step6_text"))
)

step7_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 7. CNV analysis"),
    #p("Please input the parameters for CNV analysis."),
    textInput("copykat.genome", "copykat.genome (default: 'hg20'):", value = 'hg20'),
    numericInput("copykat.LOW.DR", "copykat.LOW.DR (default: 0.05):", value = 0.05),
    numericInput("copykat.UP.DR", "copykat.UP.DR (default: 0.1):", value = 0.1),
    numericInput("copykat.win.size", "copykat.win.size (default: 25):", value = 25),
    textInput("copykat.distance", "copykat.distance (default: 'euclidean'):", value = 'euclidean'),
    numericInput("copykat.n.cores", "copykat.n.cores (default: 1):", value = 1),
    actionBttn("RunStep7_st", "Run Step7",style = "unite",color = "primary"),
    #actionButton("RunStep7", "Run Step7"),
    div(class = "spacer"), 
    uiOutput("runningStep7"),
    div(class = "spacer"),  
    uiOutput("step7_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 7. CNV analysis:"),
    uiOutput("Step7.st.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("st_step7_plot", width='auto',height = "500px"),
         textOutput("st_step7_text"))
)

step8_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 8. Deconvolution"),
    #p("Please input the parameters for deconvolution."),
    textInput("cell2loc.sc.h5ad.dir", "cell2loc.sc.h5ad.dir (default: 'NULL'):", value = 'NULL'),
    numericInput("cell2loc.sc.max.epoch", "cell2loc.sc.max.epoch (default: 1000):", value = 1000),
    numericInput("cell2loc.st.max.epoch", "cell2loc.st.max.epoch (default: 10000):", value = 10000),
    selectInput("cell2loc.use.gpu", "cell2loc.use.gpu (default: FALSE):", choices = c("FALSE" = FALSE, "TRUE" = TRUE)),
    actionBttn("RunStep8_st", "Run Step8",style = "unite",color = "primary"),
    #actionButton("RunStep8", "Run Step8"),
    div(class = "spacer"), 
    uiOutput("runningStep8"),
    div(class = "spacer"),  
    uiOutput("step8_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 8. Deconvolution:"),
    uiOutput("Step8.st.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("st_step8_plot", width='auto',height = "500px"),
         textOutput("st_step8_text"))
)

step9_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 9. Cell cycle analysis"),
    #p("Please input the parameters for cell cycle analysis."),
    textInput("s.features", HTML('The gene sets for calculating S phase scores(e.g. "gene1,gene2,gene3"):'), value = "NULL"),
    textInput("g2m.features", HTML('The gene sets for calculating G2M phase scores(e.g. "gene1,gene2,gene3"):'), value = "NULL"),
    actionBttn("RunStep9_st", "Run Step9",style = "unite",color = "primary"),
    #actionButton("RunStep9", "Run Step9"),
    div(class = "spacer"), 
    uiOutput("runningStep9"),
    div(class = "spacer"),  
    uiOutput("step9_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 9. Cell cycle analysis:"),
    uiOutput("Step9.st.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("st_step9_plot", width='auto',height = "500px"),
         textOutput("st_step9_text"))
)

step10_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    4, align = "left", h3("Step 10. Niche analysis"),
    #p("Please input the parameters for niche analysis (Please run step 8  deconvolution first)."),
    numericInput("Nich.cluster.n", "Nich.cluster.n (default: 4):", value = 4),
    actionBttn("RunStep10_st", "Run Step10",style = "unite",color = "primary"),
    #actionButton("RunStep10", "Run Step10"),
    div(class = "spacer"), 
    uiOutput("runningStep10"),
    div(class = "spacer"),  
    uiOutput("step10_completed")),
  column(
    4, align = "left",
    h3("Browse files in Step 10. Niche analysis:"),
    uiOutput("Step10.st.file_list")),
  column(4, align = "left",
         h3("Important Results Figure"),
         slickROutput("st_step10_plot", width='auto',height = "500px"),
         textOutput("st_step10_text"))
)

step11_fluidRow_st <- fluidRow(
  style = "margin-left: 10px;",
  column(
    6, align = "left", h3("Step 11. Generate the Report"),
    actionBttn("RunStep11_st", "Run Step11",style = "unite",color = "primary"),
    #actionButton("RunStep11", "Run Step11"),
    div(class = "spacer"), 
    uiOutput("runningStep11"),
    div(class = "spacer"),  
    uiOutput("step11_completed"))
)

step_st_continue_fluidRow<-fluidRow(
  style = "margin-left: 10px;",
  column(8,align = "left",
         h3("Continue Previous Analysis"),
         actionButton("st_continue_back_btn","Back to Prior Page",style = "width: 25%;",
                      class = "btn-primary btn-lg"),
         numericInput("st_jobid","Enter your Job ID",value = 0),#要是数字输入不然会报错
         selectInput("st_continue_step","Choose a step you want analysize",
                     choices = c(
                       #"Step 1. Input Data",
                       "Step 2. Quality Control" = "step 2",
                       "Step3. Normalization, PCA and Clustering" = "step 3",
                       "Step4. Differential expressed genes" = "step 4",
                       "Step5. Spatially variable features" = "step 5",
                       "Step6. Spatial interaction" = "step 6",
                       "Step 7. CNV analysis" = "step 7",
                       "Step8. Deconvolution" = "step 8",
                       "Step9. Cellcycle" = "step 9",
                       "Step10. Niche analysis" = "step 10",
                       "Step11. Generate the Report" = "step 11"
                     ),selected = "step 2")
  ))

# server---------------------------------------------------------------------------------------------------------------------------------
server = function(input, output, session){
  output$logo <- renderImage({
    list(src = '../images/hemascoper_logo.png') #加载特定位置下的图片
  }, deleteFile = FALSE) #加载后不删除
  
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
    switch(req(input$sc_step_tab),
           "step_1" = step1_fluidRow,
           "step_2" = step2_fluidRow,
           "step_3" = step3_fluidRow,
           "step_4" = step4_fluidRow,
           "step_5" = step5_fluidRow,
           "step_6" = step6_fluidRow,
           "step_7" = step7_fluidRow,
           "step_8" = step8_fluidRow,
           "step_9" = step9_fluidRow,
           "step_10" = step10_fluidRow,
           "step_11" = step11_fluidRow,
           "step_12" = step12_fluidRow,
           "step_13" = step13_fluidRow,
           "step_14" = step14_fluidRow,
           "step_15" = step15_fluidRow
    )
  })
  
  output$continue_stepContent <- renderUI({
    switch(req(input$continue_step),
           "step 2"=step2_fluidRow,
           "step 3"=step3_fluidRow,
           "step 4"=step4_fluidRow,
           "step 5"=step5_fluidRow,
           "step 6"=step6_fluidRow,
           "step 7"=step7_fluidRow,
           "step 8"=step8_fluidRow,
           "step 9"=step9_fluidRow,
           "step 10"=step10_fluidRow,
           "step 11"=step11_fluidRow,
           "step 12"=step12_fluidRow,
           "step 13"=step13_fluidRow,
           "step 14"=step14_fluidRow,
           "step 15"=step15_fluidRow)
    
  })
  
  output$stepContent_st<-renderUI({
    switch(req(input$st_step_tab),
           "step1_st" = step1_fluidRow_st,
           "step2_st" = step2_fluidRow_st,
           "step3_st" = step3_fluidRow_st,
           "step4_st" = step4_fluidRow_st,
           "step5_st" = step5_fluidRow_st,
           "step6_st" = step6_fluidRow_st,
           "step7_st" = step7_fluidRow_st,
           "step8_st" = step8_fluidRow_st,
           "step9_st" = step9_fluidRow_st,
           "step10_st" = step10_fluidRow_st,
           "step11_st" = step11_fluidRow_st)
  })
  
  output$st_continue_step <- renderUI({
    switch(req(input$st_continue_step),
           "step 2"=step2_fluidRow_st,
           "step 3"=step3_fluidRow_st,
           "step 4"=step4_fluidRow_st,
           "step 5"=step5_fluidRow_st,
           "step 6"=step6_fluidRow_st,
           "step 7"=step7_fluidRow_st,
           "step 8"=step8_fluidRow_st,
           "step 9"=step9_fluidRow_st,
           "step 10"=step10_fluidRow_st,
           "step 11"=step11_fluidRow_st)
  })
  
  ### scRNA-seq pipeline---------------------------------------------------------------------------------------------------------
  # Start----------------------------------------------------------------------------  
  observeEvent(input$start_button, {
    shinyjs::hide("ui1") #隐藏ui1 
    shinyjs::show("ui2.1") #展示ui2.1
  })
  
  observeEvent(input$sc_return_home,{
    shinyjs::hide("ui2.1")
    shinyjs::show("ui1")
  })
  
  observeEvent(input$new_analysize_btn,{
    shinyjs::hide("ui1") #隐藏ui1 
    shinyjs::show("ui2") 
    shinyjs::hide("ui2.1")
  })
  
  observeEvent(input$continue_analysize_btn,{
    shinyjs::hide("ui2")
    shinyjs::hide("ui2.1")
    shinyjs::show("ui2.2")
    output$dynamic_ui<-renderUI(
      {
        step_sc_continue_fluidRow
      }
    )
  })
  observeEvent(input$sc_continue_back_btn,{
    shinyjs::show("ui2.1")
    shinyjs::hide("ui2.2")
  })
  observeEvent(input$sc_step_tab,{
    if(input$sc_step_tab=="back_button"){
      shinyjs::hide("ui2")
      shinyjs::show("ui2.1")
    }
  })
  
  observeEvent(input$st_return_home,{
    shinyjs::hide("ui3.1")
    shinyjs::show("ui1")
  })
  observeEvent(input$st_step_tab, {
    if(input$st_step_tab=="back_button_st"){
      shinyjs::hide("ui3")
      shinyjs::show("ui3.1")
    }
  })
  
  
  #step <- reactiveVal('sc_step1')
  input.data.dirs.temp <- reactiveVal(NULL) #初始值设置为 NULL
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
  Step11_GSVA.identify.cellType.features.temp <- reactiveVal(NULL)
  Step11_GSVA.identify.diff.features.temp <- reactiveVal(NULL)
  Step11_GSVA.comparison.design.temp <- reactiveVal(NULL)
  # Step12_Construct_Trajectories = TRUE,
  Step12_Construct_Trajectories.clusters.temp <- reactiveVal(NULL)
  Step12_Construct_Trajectories.monocle.temp <- reactiveVal(NULL)
  Step12_Construct_Trajectories.slingshot.temp <- reactiveVal(NULL)
  Step12_Construct_Trajectories.scVelo.temp <- reactiveVal(NULL)
  Step13_TF_Analysis.cellTypes_colors.temp <- reactiveVal(NULL)
  Step13_TF_Analysis.groups_colors.temp <- reactiveVal(NULL)
  
  # Step1. Input data-----------------------------------------------------------------------------------------
  previous_results_path.temp <- reactiveVal(NULL)
  
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
    }
    
    if(project.names=='NULL'){project.names=NULL}else{
      project.names<-unlist(strsplit(project.names, ";"))
    }
    
    wdir <- getwd() 
    
    if(is.null(pythonPath)==FALSE){ reticulate::use_python(pythonPath) }else{print('Please set the path of Python.')}
    
    if (!file.exists(paste0(output.dir, '/HemaScopeR_results'))) {
      dir.create(paste0(output.dir, '/HemaScopeR_results'))
    }
    
    output.dir <- paste0(output.dir,'/HemaScopeR_results')
    current_time<-as.character(round(as.numeric(Sys.time()), 0)) #将jobid设置为当前的时间 2024/07/30
    
    if (!file.exists(file.path("/home/zyt/HemaScopeR_R",current_time))) {
      dir.create(file.path("/home/zyt/HemaScopeR_R",current_time)) #home/HemaScopeR_R  
    }
    
    previous_results_path <- file.path("/home/zyt/HemaScopeR_R",current_time) #每一次结果的输出路径
    
    shinyjs::runjs('$("#loadingData").text("Loading data...");')
    
    # output.dir.temp(output.dir)在 input.data.dirs 不为空且存在时，禁用多个输入控件，防止用户在数据已经存在时更改这些控件的值
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
      
      if (!file.exists(paste0(output.dir, '/Step1.Input_data'))) {
        dir.create(paste0(output.dir, '/Step1.Input_data'))
      }  
      
      file.copy(from = input.data.dirs, to = paste0(output.dir,'/Step1.Input_data'), recursive = TRUE)
      
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
        saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
      }
      
      # update variables
      previous_results_path.temp(previous_results_path)
      # update parameters 更新临时存储的路径和参数变量
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
      output$jobid<-renderText({
        paste0("Your Job ID:",current_time)
      })
      output$data_dim_output <- renderText({
        paste0("OK! Data dimensions: ", paste0(dim(sc_object)[2], ' cells, ', dim(sc_object)[1], ' genes.'))
      })
      
      
      
      shinyjs::enable("input.data.dirs") #重新启用先前禁用的输入控件
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
  observeEvent(input$jobid,{
    if(file.exists(file.path("/home/HemaScopeR_R",input$jobid))){
      
      previous_results_path<-file.path("/home/HemaScopeR_R",input$jobid)
      previous_results_path.temp(previous_results_path)
      
    }else if(!file.exists(file.path("/home/HemaScopeR_R",input$jobid))){
      showNotification("Job id does not exist",type = "warning")
    }
    
  })
  observeEvent(input$RunStep2, {
    
    if(is.null(input$jobid)){
      # load previous parameters
      input.data.dirs <- input.data.dirs.temp() 
      output.dir <- output.dir.temp() 
      project.names <- project.names.temp() 
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp() 
    Load_previous_results(previous_results_path)
    
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
    Step2_Quality_Control.RemoveBatches <- as.logical(input$Step2_Quality_Control.RemoveBatches)
    Step2_Quality_Control.RemoveDoublets <- as.logical(input$Step2_Quality_Control.RemoveDoublets)
    
    # convert some parameters
    PCs <- seq(from = as.numeric(unlist(strsplit(PCs, ":")))[1], to = as.numeric(unlist(strsplit(PCs, ":")))[2]) 
    doublerFinderwraper.PCs <- seq(from = as.numeric(unlist(strsplit(doublerFinderwraper.PCs, ":")))[1], 
                                   to = as.numeric(unlist(strsplit(doublerFinderwraper.PCs, ":")))[2])
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
    
    #Load_previous_results(previous_results_path)
    if (!file.exists(paste0(output.dir, '/Step2.Quality_control'))) {
      dir.create(paste0(output.dir, '/Step2.Quality_control'))
    }
    
    if(length(input.data.dirs) > 1){
      # preprocess and quality control for multiple scRNA-Seq data sets
      sc_object <- QC_multiple_scRNASeq(seuratObjects = input.data.list,
                                        datasetID = project.names,
                                        output.dir = paste0(output.dir,'/Step2.Quality_control'),
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
                                      output.dir = paste0(output.dir,'/Step2.Quality_control'),
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    output$step2_completed <- renderText({'Step2 completed'})
    
    output$Step2.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step2.Quality_control')  
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
      for (file in list.files(paste0(output.dir, '/Step2.Quality_control'))) {
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
    # Set the path to the folder containing the images
    img_dir <- file.path(output.dir,'Step2.Quality_control')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    
    if(length(images)!=0){
      output$step2_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step2_text<-renderText({
        print("NO Figure !")
      })
    }
    
  })
  
  # Step3. Clustering-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$RunStep3, {
    if(is.null(input$jobid)){
      # load previous parameters
      output.dir <- output.dir.temp()
      input.data.dirs <- input.data.dirs.temp()
      Step2_Quality_Control.RemoveBatches <- Step2_Quality_Control.RemoveBatches.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    
    if (!file.exists(paste0(output.dir, '/Step3.Clustering'))) {
      dir.create(paste0(output.dir, '/Step3.Clustering'))
    }
    
    if( (length(input.data.dirs) > 1) & Step2_Quality_Control.RemoveBatches ){graph.name <- 'integrated_snn'}else{graph.name <- 'RNA_snn'}
    sc_object <- FindNeighbors(sc_object, dims = PCs.clustering, k.param = n.neighbors, force.recalc = TRUE)
    sc_object <- FindClusters(sc_object, resolution = resolution, graph.name = graph.name)
    sc_object@meta.data$seurat_clusters <- as.character(as.numeric(sc_object@meta.data$seurat_clusters))
    
    # plot clustering
    pdf(paste0(paste0(output.dir,'/Step3.Clustering'), '/sc_object','tsne_cluster.pdf'), width = 6, height = 6)
    print(DimPlot(sc_object, reduction = "tsne", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1))
    dev.off()
    
    pdf(paste0(paste0(output.dir,'/Step3.Clustering'), '/sc_object','umap_cluster.pdf'), width = 6, height = 6)
    print(DimPlot(sc_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1))
    dev.off()
    
    png(paste0(paste0(output.dir,'/Step3.Clustering'), '/sc_object','tsne_cluster.png'), width = 600, height = 600)
    print(DimPlot(sc_object, reduction = "tsne", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1))
    dev.off()
    
    png(paste0(paste0(output.dir,'/Step3.Clustering'), '/sc_object','umap_cluster.png'), width = 600, height = 600)
    print(DimPlot(sc_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.1))
    dev.off()    
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }  
    
    #update parameters
    PCs.clustering.temp(PCs.clustering)
    n.neighbors.temp(n.neighbors)
    resolution.temp(resolution)
    
    shinyjs::enable('PCs.clustering')
    shinyjs::enable('n.neighbors')
    shinyjs::enable('resolution')
    shinyjs::enable('RunStep3')
    output$step3_completed <- renderText({'Step3 completed'})
    output$Step3.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step3.Clustering')  # 指定的文件夹路径
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
      for (file in list.files(paste0(output.dir, '/Step3.Clustering'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step3.Clustering/', file))
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
    #show figures
    img_dir <- file.path(output.dir,'Step3.Clustering')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      
      output$step3_plot<-renderSlickR({
        
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step3_text<-renderText({
        print("NO Figure !")
      })
    }
  })
  
  # Step4. Identify cell types automatically-----------------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$RunStep4, {
    
    if(is.null(input$jobid)){
      output.dir <- output.dir.temp()
      PCs <- PCs.temp()
      databasePath <- databasePath.temp()
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    
    if (!file.exists(paste0(output.dir, '/Step4.Identify_Cell_Types'))) {
      dir.create(paste0(output.dir, '/Step4.Identify_Cell_Types'))
    } 
    
    sc_object <- run_cell_annotation(object = sc_object, 
                                     assay = 'RNA', 
                                     species = Org,
                                     output.dir = paste0(output.dir,'/Step4.Identify_Cell_Types'))
    
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
                                  output.dir = paste0(output.dir, '/Step4.Identify_Cell_Types'))
        DefaultAssay(sc_object) <- 'integrated'
      }else{
        sc_object <- mapDataToRef(ref_object = HematoMap.reference,
                                  ref_labels = HematoMap.reference@meta.data$CellType,
                                  query_object = sc_object,
                                  PCs = PCs,
                                  output.dir = paste0(output.dir, '/Step4.Identify_Cell_Types'))
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
             save_path=paste0(output.dir,'/Step4.Identify_Cell_Types'),
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    output$step4_completed <- renderText({'Step4 completed'})
    
    output$Step4.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step4.Identify_Cell_Types')  # 指定的文件夹路径
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
      for (file in list.files(paste0(output.dir, '/Step4.Identify_Cell_Types'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step4.Identify_Cell_Types/', file))
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
    #show figures
    img_dir <- file.path(output.dir,'Step4.Identify_Cell_Types')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step4_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step4_text<-renderText({
        print("NO Figure !")
      })
    }
  })
  
  # Step5. Visualization-----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep5, {
    
    if(is.null(input$jobid)){
      
      # load previous parameters
      output.dir <- output.dir.temp()
      input.data.dirs <- input.data.dirs.temp()
      Step2_Quality_Control.RemoveBatches <- Step2_Quality_Control.RemoveBatches.temp()
      pythonPath <- pythonPath.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    
    if (!file.exists(paste0(output.dir, '/Step5.Visualization'))) {
      dir.create(paste0(output.dir, '/Step5.Visualization'))
    }    
    # run phateR
    if( (length(input.data.dirs) > 1) & Step2_Quality_Control.RemoveBatches ){
      DefaultAssay(sc_object) <- 'integrated'
    }else{
      DefaultAssay(sc_object) <- 'RNA'}
    
    if(!is.null(pythonPath)){
      run_phateR(sc_object = sc_object,
                 output.dir = paste0(output.dir,'/Step5.Visualization'),
                 pythonPath = pythonPath,
                 phate.knn = phate.knn,
                 phate.npca = phate.npca,
                 phate.t = phate.t,
                 phate.ndim = phate.ndim)     
    }
    
    # plot cell types
    pdf(paste0(paste0(output.dir,'/Step5.Visualization'), '/sc_object ','tsne cell types.pdf'), width = 6, height = 6)
    print(DimPlot(sc_object, reduction = "tsne", group.by = "ident", label = FALSE, pt.size = 0.1))
    dev.off()
    
    pdf(paste0(paste0(output.dir,'/Step5.Visualization'), '/sc_object ','umap cell types.pdf'), width = 6, height = 6)
    print(DimPlot(sc_object, reduction = "umap", group.by = "ident", label = FALSE, pt.size = 0.1))
    dev.off()
    
    png(paste0(paste0(output.dir,'/Step5.Visualization'), '/sc_object ','tsne cell types.png'), width = 600, height = 600)
    print(DimPlot(sc_object, reduction = "tsne", group.by = "ident", label = FALSE, pt.size = 0.1))
    dev.off()
    
    png(paste0(paste0(output.dir,'/Step5.Visualization'), '/sc_object ','umap cell types.png'), width = 600, height = 600)
    print(DimPlot(sc_object, reduction = "umap", group.by = "ident", label = FALSE, pt.size = 0.1))
    dev.off()    
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    output$step5_completed <- renderText({'Step5 completed'})
    
    output$Step5.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step5.Visualization')  # 指定的文件夹路径  '/Step5.Visualization'
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(paste0(output.dir, '/Step5.Visualization'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step5.Visualization/', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step5.Visualization')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step5_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step5_text<-renderText({
        print("NO Figure !")
      })
    }
  })
  
  # Step6. Find DEGs-----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep6, {
    
    if(is.null(input$jobid)){
      
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    
    # set parameters
    min.pct <- input$min.pct
    logfc.threshold <- input$logfc.threshold
    
    shinyjs::runjs('$("#runningStep6").text("Running Step6...");')  
    shinyjs::disable('min.pct')
    shinyjs::disable('logfc.threshold')
    shinyjs::disable('RunStep6')
    
    if (!file.exists(paste0(output.dir, '/Step6.Find_DEGs'))) {
      dir.create(paste0(output.dir, '/Step6.Find_DEGs'))
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
                         output.dir=paste0(output.dir, '/Step6.Find_DEGs'))
    
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
      saveRDS(var, file = paste0(previous_results_path,'/',var_name,".rds"))  # Save as RDS with the variable's name
    }  
    
    #update parameters
    min.pct.temp(min.pct)
    logfc.threshold.temp(logfc.threshold)
    
    shinyjs::enable('min.pct')
    shinyjs::enable('logfc.threshold')
    shinyjs::enable('RunStep6')
    output$step6_completed <- renderText({'Step6 completed'})
    output$Step6.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step6.Find_DEGs')  # 指定的文件夹路径  '/Step5.Visualization'
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(paste0(output.dir, '/Step6.Find_DEGs'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step6.Find_DEGs/', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step6.Find_DEGs')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step6_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step6_text<-renderText({
        print("NO Figure !")
      })
    }
  })
  
  # Step7. Assign Cell Cycles-----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep7, {
    
    if(is.null(input$jobid)){
      
      databasePath <- databasePath.temp()
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    
    # set parameters
    cellcycleCutoff <- input$cellcycleCutoff
    
    # convert some parameters
    if(cellcycleCutoff=='NULL'){cellcycleCutoff <- NULL}else{cellcycleCutoff <- as.numeric(cellcycleCutoff)}
    
    shinyjs::runjs('$("#runningStep7").text("Running Step7...");')  
    shinyjs::disable('cellcycleCutoff')
    shinyjs::disable('RunStep7')
    
    if (!file.exists(paste0(output.dir, '/Step7.Assign_cell_cycles'))) {
      dir.create(paste0(output.dir, '/Step7.Assign_cell_cycles'))
    }
    datasets.before.batch.removal <- readRDS(file.path(previous_results_path,'sc_object.rds'))  
    sc_object <- cellCycle(sc_object=sc_object,
                           counts_matrix = GetAssayData(object = datasets.before.batch.removal, slot = "counts")%>%as.matrix(),
                           data_matrix = GetAssayData(object = datasets.before.batch.removal, slot = "data")%>%as.matrix(),
                           cellcycleCutoff = cellcycleCutoff,
                           cellTypeOrders = unique(sc_object@meta.data$selectLabels),
                           output.dir=paste0(output.dir, '/Step7.Assign_cell_cycles'),
                           databasePath = databasePath,
                           Org = Org)
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path,'/',var_name,".rds"))  # Save as RDS with the variable's name
    }  
    
    #update parameters
    cellcycleCutoff.temp(cellcycleCutoff)
    
    shinyjs::enable('cellcycleCutoff')
    shinyjs::enable('RunStep7')
    output$step7_completed <- renderText({'Step7 completed'})
    
    #download files
    output$Step7.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step7.Assign_cell_cycles')  # 指定的文件夹路径  '/Step5.Visualization'
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(paste0(output.dir, '/Step7.Assign_cell_cycles'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step7.Assign_cell_cycles/', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step7.Assign_cell_cycles')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step7_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step7_text<-renderText({
        print("NO Figure !")
      })
    }
  })  
  
  # Step8. Calculate Heterogeneity-----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep8, {
    
    if(is.null(input$jobid)){
      
      output.dir <- output.dir.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    
    # set parameters
    ViolinPlot.cellTypeOrders <- input$ViolinPlot.cellTypeOrders
    
    # convert some parameters
    if(ViolinPlot.cellTypeOrders=='NULL'){ViolinPlot.cellTypeOrders <- NULL}else{ViolinPlot.cellTypeOrders <- unlist(strsplit(ViolinPlot.cellTypeOrders, ","))}
    
    shinyjs::runjs('$("#runningStep8").text("Running Step8...");')
    shinyjs::disable('ViolinPlot.cellTypeOrders')
    shinyjs::disable('RunStep8')
    
    if (!file.exists(paste0(output.dir, '/Step8.Calculate_heterogeneity'))) {
      dir.create(paste0(output.dir, '/Step8.Calculate_heterogeneity'))
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
                  output.dir = paste0(output.dir, '/Step8.Calculate_heterogeneity'))
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }  
    0  
    #update parameters
    ViolinPlot.cellTypeOrders.temp(ViolinPlot.cellTypeOrders)
    
    shinyjs::enable('ViolinPlot.cellTypeOrders')
    shinyjs::enable('RunStep8')
    output$step8_completed <- renderText({'Step8 completed'})
    #download files
    output$Step8.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step8.Calculate_heterogeneity')  # 指定的文件夹路径  '/Step5.Visualization'
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(paste0(output.dir, '/Step8.Calculate_heterogeneity'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step8.Calculate_heterogeneity/', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step8.Calculate_heterogeneity')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step8_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step8_text<-renderText({
        print("NO Figure !")
      })
    }
  })    
  
  # Step9. Violin Plot for Marker Genes-----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep9, {
    
    if(is.null(input$jobid)){
      
      # load previous parameters
      output.dir <- output.dir.temp()
      input.data.dirs <- input.data.dirs.temp()
      databasePath <- databasePath.temp()
      Step2_Quality_Control.RemoveBatches <- Step2_Quality_Control.RemoveBatches.temp()
      Org <- Org.temp()
      ViolinPlot.cellTypeOrders <- ViolinPlot.cellTypeOrders.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
                       output.dir = paste0(output.dir, '/Step9.Violin_plot_for_marker_genes'),
                       databasePath = databasePath)
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }  
    
    #update parameters
    marker.genes.temp(marker.genes)
    ViolinPlot.cellTypeColors.temp(ViolinPlot.cellTypeColors)
    
    shinyjs::enable('marker.genes')
    shinyjs::enable('ViolinPlot.cellTypeColors')
    shinyjs::enable('RunStep9')
    output$step9_completed <- renderText({'Step9 completed'})
    #download files
    output$Step9.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step9.Violin_plot_for_marker_genes')  # 指定的文件夹路径  '/Step5.Visualization'
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(paste0(output.dir, '/Step9.Violin_plot_for_marker_genes'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step9.Violin_plot_for_marker_genes/', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step9.Violin_plot_for_marker_genes')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step9_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step9_text<-renderText({
        print("NO Figure !")
      })
    }
  })      
  
  # Step10. Calculate Lineage Scores-----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep10, {
    
    if(is.null(input$jobid)){
      
      # load previous parameters
      output.dir <- output.dir.temp()
      databasePath <- databasePath.temp()
      ViolinPlot.cellTypeOrders <- ViolinPlot.cellTypeOrders.temp()
      ViolinPlot.cellTypeColors <- ViolinPlot.cellTypeColors.temp()
      Org <- Org.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    if (!file.exists(paste0(output.dir, '/Step10.Calculate_lineage_scores'))) {
      dir.create(paste0(output.dir, '/Step10.Calculate_lineage_scores'))
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
                  output.dir = paste0(output.dir, '/Step10.Calculate_lineage_scores'),
                  databasePath = databasePath)
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }  
    
    #update parameters
    lineage.genelist.temp(lineage.genelist)
    lineage.names.temp(lineage.names)
    groups_colors.temp(groups_colors)
    
    shinyjs::enable('lineage.genelist')
    shinyjs::enable('lineage.names')
    shinyjs::enable('groups_colors')
    shinyjs::enable('RunStep10')
    output$step10_completed <- renderText({'Step10 completed'})
    #download files
    output$Step10.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step10.Calculate_lineage_scores')  # 指定的文件夹路径  '/Step5.Visualization'
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(paste0(output.dir, '/Step10.Calculate_lineage_scores'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step10.Calculate_lineage_scores/', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step10.Calculate_lineage_scores')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step10_plot<-renderSlickR({
        slickR(images,width = "80%",height = 160) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step10_text<-renderText({
        print("NO Figure !")
      })
    }
  })
  
  # Step11. GSVA-----------------------------------------------------------------------------------------------------------------------------   
  
  observeEvent(input$RunStep11, {
    
    if(is.null(input$jobid)){
      
      # load previous parameters
      databasePath <- databasePath.temp()
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      ViolinPlot.cellTypeOrders <- ViolinPlot.cellTypeOrders.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    
    if (!file.exists(paste0(output.dir, '/Step11.GSVA'))) {
      dir.create(paste0(output.dir, '/Step11.GSVA'))
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
             output.dir = paste0(output.dir, '/Step11.GSVA'))
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    #download files
    output$Step11.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step11.GSVA')  # 指定的文件夹路径  '/Step5.Visualization'
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(paste0(output.dir, '/Step11.GSVA'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step11.GSVA/', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step11.GSVA')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step11_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step11_text<-renderText({
        print("NO Figure !")
      })
    }
  })
  
  # Step12. Construct Trajectories-----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep12, {
    
    if(is.null(input$jobid)){
      
      # load previous parameters
      input.data.dirs <- input.data.dirs.temp()
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      pythonPath <- pythonPath.temp()
      databasePath <- databasePath.temp()
      ViolinPlot.cellTypeColors <- ViolinPlot.cellTypeColors.temp()
      Step2_Quality_Control.RemoveBatches <- Step2_Quality_Control.RemoveBatches.temp()
      PCs <- PCs.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
      loom.files.path <- gsub("^([^/].*)$", "/\\1", loom.files.path)
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
    
    if (!file.exists(paste0(output.dir, '/Step12.Construct_trajectories'))) {
      dir.create(paste0(output.dir, '/Step12.Construct_trajectories'))
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
      if (!file.exists(paste0(output.dir, '/Step12.Construct_trajectories/monocle2'))) {
        dir.create(paste0(output.dir, '/Step12.Construct_trajectories/monocle2'))
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
                  output.dir = paste0(output.dir, '/Step12.Construct_trajectories/monocle2'))
    }
    
    if(Step12_Construct_Trajectories.slingshot == TRUE){
      # slingshot
      if (!file.exists(paste0(output.dir, '/Step12.Construct_trajectories/slingshot'))) {
        dir.create(paste0(output.dir, '/Step12.Construct_trajectories/slingshot'))
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
                    output.dir = paste0(output.dir, '/Step12.Construct_trajectories/slingshot'))
    }
    
    if(Step12_Construct_Trajectories.scVelo == TRUE){
      # scVelo
      if((!is.null(loom.files.path))&(!is.null(pythonPath))){
        if (!file.exists(paste0(output.dir, '/Step12.Construct_trajectories/scVelo'))) {
          dir.create(paste0(output.dir, '/Step12.Construct_trajectories/scVelo'))
        } 
        
        prepareDataForScvelo(sc_object = sc_object.subset,
                             loom.files.path = loom.files.path,
                             scvelo.reduction = 'pca',
                             scvelo.column = 'selectLabels',
                             output.dir = paste0(output.dir, '/Step12.Construct_trajectories/scVelo'))
        
        reticulate::py_run_string(paste0("import os\noutputDir = '", output.dir, "'"))
        reticulate::py_run_file("~/HemaScopeR/R/sc_run_scvelo.py", convert = FALSE)
      }
    }
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    output$step12_completed <- renderText({'Step12 completed'})
    #download files
    output$Step12.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step12.Construct_trajectories')  # 指定的文件夹路径  '/Step5.Visualization'
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(paste0(output.dir, '/Step12.Construct_trajectories'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step12.Construct_trajectories/', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step12.Construct_trajectories')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step12_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step12_text<-renderText({
        print("NO Figure !")
      })
    }
  })
  
  # Step13. TF Analysis-----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep13, {
    
    if(is.null(input$jobid)){
      
      # load previous parameters
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      pythonPath <- pythonPath.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    
    # set parameters
    Step13_TF_Analysis.cellTypes_colors <- input$Step13_TF_Analysis.cellTypes_colors
    Step13_TF_Analysis.groups_colors <- input$Step13_TF_Analysis.groups_colors
    
    shinyjs::runjs('$("#runningStep13").text("Running Step13...");')  
    shinyjs::disable('Step13_TF_Analysis.cellTypes_colors')
    shinyjs::disable('Step13_TF_Analysis.groups_colors')
    shinyjs::disable('RunStep13')
    
    if (!file.exists(paste0(output.dir, '/Step13.TF_analysis'))) {
      dir.create(paste0(output.dir, '/Step13.TF_analysis'))
    }
    
    run_SCENIC(countMatrix = countsSlot,
               cellTypes = sc_object@meta.data$selectLabels,
               datasetID = sc_object@meta.data$datasetID,
               cellTypes_colors = Step13_TF_Analysis.cellTypes_colors,
               cellTypes_orders = unique(sc_object@meta.data$selectLabels),
               groups_colors = Step13_TF_Analysis.groups_colors,
               groups_orders = unique(sc_object@meta.data$datasetID),
               Org = Org,
               output.dir = paste0(output.dir, '/Step13.TF_analysis'),
               pythonPath = pythonPath,
               databasePath = databasePath)
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }  
    
    #update parameters
    Step13_TF_Analysis.cellTypes_colors.temp(Step13_TF_Analysis.cellTypes_colors)
    Step13_TF_Analysis.groups_colors.temp(Step13_TF_Analysis.groups_colors)
    
    shinyjs::enable('Step13_TF_Analysis.cellTypes_colors')
    shinyjs::enable('Step13_TF_Analysis.groups_colors')
    shinyjs::enable('RunStep13')
    output$step13_completed <- renderText({'Step13 completed'})
    #download files
    output$Step13.file_list <- renderUI({
      folder_path <- paste0(output.dir, '/Step13.TF_analysis')  # 指定的文件夹路径  '/Step5.Visualization'
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(paste0(output.dir, '/Step13.TF_analysis'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(paste0(output.dir, '/Step13.TF_analysis/', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step13.TF_analysis')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step13_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step13_text<-renderText({
        print("NO Figure !")
      })
    }
  })  
  
  # Step14. -----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep14, {
    
    if(is.null(input$jobid)){
      
      # load previous parameters
      output.dir <- output.dir.temp()
      Org <- Org.temp()
      ViolinPlot.cellTypeOrders <- ViolinPlot.cellTypeOrders.temp()
      ViolinPlot.cellTypeColors <- ViolinPlot.cellTypeColors.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    
    # set parameters
    sorting <- as.logical(input$sorting)
    
    shinyjs::runjs('$("#runningStep14").text("Running Step14...");')  
    shinyjs::disable('sorting')
    shinyjs::disable('RunStep14')
    
    if (!file.exists(paste0(output.dir, '/Step14.Cell_cell_interection'))) {
      dir.create(paste0(output.dir, '/Step14.Cell_cell_interection'))
    }
    tempwd <- getwd()
    run_CellChat(data.input=countsSlot,
                 labels = sc_object@meta.data$selectLabels,
                 cell.orders = ViolinPlot.cellTypeOrders,
                 cell.colors = ViolinPlot.cellTypeColors,
                 sample.names = rownames(sc_object@meta.data),
                 Org = Org,
                 sorting = sorting,
                 output.dir = paste0(output.dir, '/Step14.Cell_cell_interection'))
    setwd(tempwd)
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }  
    
    #update parameters
    sorting.temp(sorting)
    
    shinyjs::enable('sorting')
    shinyjs::enable('RunStep14')
    output$step14_completed <- renderText({'Step14 completed'})
    
    output$Step14.file_list <- renderUI({
      # dirs
      if (dir.exists(paste0(output.dir, '/Step14.Cell_cell_interection'))) {
        # files
        files <- list.files(paste0(output.dir, '/Step14.Cell_cell_interection'))
        
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
    #show figures
    img_dir <- file.path(output.dir,'Step14.Cell_cell_interection')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    
    if(length(images)!=0){
      image_names <- basename(images)
      
      output$step14_plot<-renderSlickR({
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$step14_text<-renderText({
        print("NO Figure !")
      })
    }
  })
  
  # Step15. -----------------------------------------------------------------------------------------------------------------------------  
  
  observeEvent(input$RunStep15, {
    
    if(is.null(input$jobid)){
      
      # load previous parameters
      output.dir <- output.dir.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }  
    
    #update parameters
    
    shinyjs::enable('RunStep15')
    output$step15_completed <- renderText({'Step15 completed'})
  })
  
  ### st-seq pipeline---------------------------------------------------------------------------------------------------------------
  observeEvent(input$start_button_st, {
    shinyjs::hide("ui1")
    shinyjs::show("ui3.1")
  })
  observeEvent(input$st_new_analysize_btn,{
    shinyjs::show("ui3")
    shinyjs::hide("ui3.1")
  })
  observeEvent(input$st_continue_analysize_btn,{
    shinyjs::hide("ui3")
    shinyjs::hide("ui3.1")
    shinyjs::show("ui3.2")
    output$dynamic_st_ui<-renderUI(
      {
        step_st_continue_fluidRow
      }
      
    )
  })
  observeEvent(input$st_continue_back_btn,{
    shinyjs::hide("ui3.2")
    shinyjs::show("ui3.1")
  })
  
  input.data.dir.temp <- reactiveVal(NULL)
  #output.dir.temp <- reactiveVal(NULL)
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
  
  # Step1. Loading data---------------------------------------------------------------------------------
  observeEvent(input$RunStep1_st, {
    
    # set parameters
    input.data.dir <- input$input.data.dir
    output.dir <- input$output.dir
    sampleName <- input$sampleName
    pythonPath <- input$pythonPath
    
    if(!file.exists(file.path(output.dir, sampleName))){
      dir.create(file.path(output.dir, sampleName))
    }else{
      warning(paste0('The new results will overwrite the file under ',
                     file.path(output.dir, sampleName)))
    }
    output.dir <- file.path(output.dir, sampleName)
    current_time<-as.character(round(as.numeric(Sys.time()),0))
    #print(paste0('The results will be saved in ', output.dir))
    if (!file.exists(file.path("/home/zyt/HemaScopeR_R",current_time))) {
      dir.create(file.path("/home/zyt/HemaScopeR_R",current_time))
    }
    
    previous_results_path <- file.path("/home/zyt/HemaScopeR_R",current_time)
    
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }
    
    output$st_jobid_1<-renderText({
      paste0("Your Job id:",current_time)
    })
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
  observeEvent(input$st_jobid,{
    
    if(file.exists(file.path("/home/HemaScopeR_R",input$st_jobid))){
      
      previous_results_path<-file.path("/home/HemaScopeR_R",input$st_jobid)
      
      print(paste0("previous_results_path:",previous_results_path))
      previous_results_path.temp(previous_results_path)
      print(paste0("previous_results_path.temp:",previous_results_path.temp()))
      
    }else if(!file.exists(file.path("/home/HemaScopeR_R",input$st_jobid))){
      showNotification("Job id does not exist",type = "warning")
    }
    
  })
  
  observeEvent(input$RunStep2_st, {
    
    if(is.null(input$st_jobid)){
      # load previous parameters
      output.dir <- output.dir.temp()
    }
    
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    
    shinyjs::runjs('$("#runningStep2").text("Running Step2...");')  
    shinyjs::disable('min.gene')
    shinyjs::disable('min.nUMI')
    shinyjs::disable('max.gene')
    shinyjs::disable('max.nUMI')
    shinyjs::disable('min.spot')
    shinyjs::disable('bool.remove.mito')
    shinyjs::disable('species')
    shinyjs::disable('RunStep2')
    
    if(!file.exists(file.path(output.dir, 'Step2_QC'))){
      dir.create(file.path(output.dir, 'Step2_QC'))
    }
    
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }
    
    # update parameters
    min.gene.temp(min.gene)
    min.nUMI.temp(min.nUMI)
    max.gene.temp(max.gene)
    max.nUMI.temp(max.nUMI)
    min.spot.temp(min.spot)
    bool.remove.mito.temp(bool.remove.mito)
    species.temp(species)
    #RunStep2.temp(RunStep2)
    
    shinyjs::enable('min.gene')
    shinyjs::enable('min.nUMI')
    shinyjs::enable('max.gene')
    shinyjs::enable('max.nUMI')
    shinyjs::enable('min.spot')
    shinyjs::enable('bool.remove.mito')
    shinyjs::enable('species')
    shinyjs::enable('RunStep2')
    output$step2_completed <- renderText({'Step2 completed.'})
    
    #download files
    output$Step2.st.file_list <- renderUI({
      folder_path <- file.path(output.dir,'Step2_QC')  # 
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(file.path(output.dir, 'Step2_QC'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(file.path(output.dir, 'Step2_QC', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step2_QC/png')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    if(length(images)!=0){
      output$st_step2_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$st_step2_text<-renderText({
        print("NO Figure !")
      })
    }
    
  })
  
  # Step3. Normalization, PCA and Clustering---------------------------------------------------------------------------------
  observeEvent(input$RunStep3_st, {
    if(is.null(input$st_jobid)){
      # load previous parameters
      output.dir <- output.dir.temp()
    }
    
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    
    if(!file.exists(file.path(output.dir, 'Step3_Clustering'))){
      dir.create(file.path(output.dir, 'Step3_Clustering'))
    }
    
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    output$step3_completed <- renderText({'Step3 completed'})
    #download files
    output$Step3.st.file_list <- renderUI({
      folder_path <- file.path(output.dir,'Step3_Clustering')  
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
      for (file in list.files(file.path(output.dir,'Step3_Clustering'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(file.path(output.dir,'Step3_Clustering', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step3_Clustering/png')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    if(length(images)!=0){
      output$st_step3_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$st_step3_text<-renderText({
        print("NO Figure !")
      })
    }
    
  })
  
  # Step4. Differential expressed genes---------------------------------------------------------------------------------
  observeEvent(input$RunStep4_st, {
    if(is.null(input$st_jobid)){
      # load previous parameters
      output.dir <- output.dir.temp()
    }
    
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    # set parameters
    only.pos <- as.logical(input$only.pos)
    min.pct <- input$min.pct
    logfc.threshold <- input$logfc.threshold
    test.use <- input$test.use
    
    # convert some parameters
    shinyjs::runjs('$("#runningStep4").text("Running Step4...");')  
    shinyjs::disable('only.pos')
    shinyjs::disable('min.pct')
    shinyjs::disable('logfc.threshold')
    shinyjs::disable('test.use')
    shinyjs::disable('RunStep4')
    
    if(!file.exists(file.path(output.dir, 'Step4_Find_DEGs'))){
      dir.create(file.path(output.dir, 'Step4_Find_DEGs'))
    }
    
    st.markers <- st_Find_DEGs(
      st_obj = st_obj,
      output.dir = file.path(output.dir, 'Step4_Find_DEGs'),
      ident.label = 'seurat_clusters',
      only.pos = only.pos,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      test.use = test.use,
      verbose = FALSE
    )
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    output$step4_completed <- renderText({'Step4 completed'})
    #download files
    output$Step3.st.file_list <- renderUI({
      folder_path <- file.path(output.dir,'Step4_Find_DEGs')  
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
      for (file in list.files(file.path(output.dir,'Step4_Find_DEGs'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(file.path(output.dir,'Step4_Find_DEGs', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step4_Find_DEGs/png')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    if(length(images)!=0){
      output$st_step4_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$st_step4_text<-renderText({
        print("NO Figure !")
      })
    }
    
  })
  
  # Step5. Spatially variable features---------------------------------------------------------------------------------
  observeEvent(input$RunStep5_st, {
    if(is.null(input$st_jobid)){
      # load previous parameters
      output.dir <- output.dir.temp()
    }
    
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    
    # set parameters
    selection.method <- input$selection.method
    n.top.show <- input$n.top.show
    n.col.show <- input$n.col.show
    
    # convert some parameters
    shinyjs::runjs('$("#runningStep5").text("Running Step5...");')  
    shinyjs::disable('selection.method')
    shinyjs::disable('n.top.show')
    shinyjs::disable('n.col.show')
    shinyjs::disable('RunStep5')
    
    if(!file.exists(file.path(output.dir, 'Step5_SpatiallyVariableFeatures'))){
      dir.create(file.path(output.dir, 'Step5_SpatiallyVariableFeatures'))
    }
    
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }
    
    # update parameters
    selection.method.temp(selection.method)
    n.top.show.temp(n.top.show)
    n.col.show.temp(n.col.show)
    
    shinyjs::enable('selection.method')
    shinyjs::enable('n.top.show')
    shinyjs::enable('n.col.show')
    shinyjs::enable('RunStep5')
    output$step5_completed <- renderText({'Step5 completed'})
    
    #download files
    output$Step5.st.file_list <- renderUI({
      folder_path <- file.path(output.dir,'Step5_SpatiallyVariableFeatures')  # 
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(file.path(output.dir,'Step5_SpatiallyVariableFeatures'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(file.path(output.dir,'Step5_SpatiallyVariableFeatures', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step5_SpatiallyVariableFeatures/png')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    if(length(images)!=0){
      output$st_step5_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$st_step5_text<-renderText({
        print("NO Figure !")
      })
    }
    
  })
  
  # Step6. Spatial interaction---------------------------------------------------------------------------------
  observeEvent(input$RunStep6_st, {
    if(is.null(input$st_jobid)){
      
      # load previous parameters
      input.data.dir <- input.data.dir.temp()
      species <- species.temp()
      output.dir <- output.dir.temp()
      
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    
    # set parameters
    commot.signaling_type <- input$commot.signaling_type
    commot.database <- input$commot.database
    commot.min_cell_pct <- input$commot.min_cell_pct
    commot.dis_thr <- input$commot.dis_thr
    commot.n_permutations <- input$commot.n_permutations
    
    # convert some parameters
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
    st_Interaction(
      st_data_path = input.data.dir,
      metadata_path = file.path(interaction_path, 'metadata.csv'),
      label_key = 'seurat_clusters',
      save_path = interaction_path,
      species = species,
      signaling_type = commot.signaling_type,
      database = commot.database,
      min_cell_pct = commot.min_cell_pct,
      dis_thr = commot.dis_thr,
      n_permutations = commot.n_permutations,
      pythonPath = pythonPath
      #condaenv = condaenv
    )
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    output$step6_completed <- renderText({'Step6 completed'})
    #download files
    output$Step6.st.file_list <- renderUI({
      folder_path <- file.path(output.dir,'Step6_Interaction')  # 
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(file.path(output.dir,'Step6_Interaction'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(file.path(output.dir,'Step6_Interaction', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step6_Interaction/png')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    if(length(images)!=0){
      output$st_step6_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$st_step6_text<-renderText({
        print("NO Figure !")
      })
    }
    
  })
  
  # Step7. Loading data---------------------------------------------------------------------------------
  observeEvent(input$RunStep7_st, {
    if(is.null(input$st_jobid)){
      # load previous parameters
      species <- species.temp()
      output.dir <- output.dir.temp()
      
    }
    
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    
    if(!file.exists(file.path(output.dir, 'Step7_CNV_analysis'))){
      dir.create(file.path(output.dir, 'Step7_CNV_analysis'))
    }
    
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    output$step7_completed <- renderText({'Step7 completed'})
    
    #download files
    output$Step7.st.file_list <- renderUI({
      folder_path <- file.path(output.dir,'Step7_CNV_analysis')  # 
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(file.path(output.dir,'Step7_CNV_analysis'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(file.path(output.dir,'Step7_CNV_analysis', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step7_CNV_analysis/png')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    if(length(images)!=0){
      output$st_step7_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$st_step7_text<-renderText({
        print("NO Figure !")
      })
    }
    
  })
  
  # Step8. Deconvolution---------------------------------------------------------------------------------
  observeEvent(input$RunStep8_st, {
    
    if(is.null(input$st_jobid)){
      # load previous parameters
      input.data.dir <- input.data.dir.temp()
      species <- species.temp()
      output.dir <- output.dir.temp()
    }
    
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    
    if(!file.exists(file.path(output.dir, 'Step8_Deconvolution'))){
      dir.create(file.path(output.dir, 'Step8_Deconvolution'))
    }
    
    st_obj <- st_Deconvolution(
      st.data.dir = input.data.dir,
      sc.h5ad.dir = cell2loc.sc.h5ad.dir,
      st_obj = st_obj,
      save_path = file.path(output.dir, 'Step8_Deconvolution'),
      sc.labels.key = 'seurat_clusters',
      species = species,
      sc.max.epoch = cell2loc.sc.max.epoch,
      st.max.epoch = cell2loc.st.max.epoch,
      use.gpu = cell2loc.use.gpu,
      pythonPath = pythonPath,
      use.Dataset = 'LymphNode'
      # condaenv = condaenv
    )
    
    # Get the names of all variables in the current environment
    variable_names <- ls()
    # Loop through the variable names and save them as RDS files
    for (var_name in variable_names) {
      var <- get(var_name)  # Get the variable by its name
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
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
    output$step8_completed <- renderText({'Step8 completed'})
    #download files
    output$Step8.st.file_list <- renderUI({
      folder_path <- file.path(output.dir,'Step8_Deconvolution')   
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(file.path(output.dir,'Step8_Deconvolution'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(file.path(output.dir,'Step8_Deconvolution', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step8_Deconvolution')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    if(length(images)!=0){
      output$st_step8_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$st_step8_text<-renderText({
        print("NO Figure !")
      })
    }
    
  })
  
  # Step9. Cellcycle---------------------------------------------------------------------------------
  observeEvent(input$RunStep9_st, {
    
    if(is.null(input$st_jobid)){
      
      # load previous parameters
      output.dir <- output.dir.temp()
      species <- species.temp()
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
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
    
    if(!file.exists(file.path(output.dir, 'Step9_Cellcycle'))){
      dir.create(file.path(output.dir, 'Step9_Cellcycle'))
    }
    
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }
    
    # update parameters
    s.features.temp(s.features)
    g2m.features.temp(g2m.features)
    
    shinyjs::enable('s.features')
    shinyjs::enable('g2m.features')
    shinyjs::enable('RunStep9')
    output$step9_completed <- renderText({'Step9 completed'})
    #download files
    output$Step9.st.file_list <- renderUI({
      folder_path <- file.path(output.dir,'Step9_Cellcycle')  # 
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(file.path(output.dir,'Step9_Cellcycle'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(file.path(output.dir,'Step9_Cellcycle', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step9_Cellcycle')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    if(length(images)!=0){
      output$st_step9_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$st_step9_text<-renderText({
        print("NO Figure !")
      })
    }
    
  })
  
  # Step10. Niche analysis---------------------------------------------------------------------------------
  observeEvent(input$RunStep10_st, {
    
    if(is.null(input$st_jobid)){
      
      # load previous parameters
      output.dir <- output.dir.temp()
      slice <- slice.temp()
      species <- species.temp()
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    
    # set parameters
    Nich.cluster.n <- input$Nich.cluster.n
    
    # convert some parameters
    shinyjs::runjs('$("#runningStep10").text("Running Step10...");')  
    shinyjs::disable('Nich.cluster.n')
    shinyjs::disable('RunStep10')
    
    if(!file.exists(file.path(output.dir, 'Step10_NicheAnalysis'))){
      dir.create(file.path(output.dir, 'Step10_NicheAnalysis'))
    }
    
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }
    
    # update parameters
    Nich.cluster.n.temp(Nich.cluster.n)
    
    shinyjs::enable('Nich.cluster.n')
    shinyjs::enable('RunStep10')
    output$step10_completed <- renderText({'Step10 completed'})
    #download files
    output$Step10.st.file_list <- renderUI({
      folder_path <- file.path(output.dir,'Step10_NicheAnalysis')  # 
      
      if (dir.exists(folder_path)) {
        files <- list.files(folder_path)
        
        if (length(files) > 0) {
          file_links <- lapply(files, function(file) {
            downloadLink(
              outputId = paste0("download_", file),
              label = file,
              href = file.path(folder_path, file)
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
    
    # Listen to the click event on the file links  # '/Step5.Visualization'
    observe({
      for (file in list.files(file.path(output.dir,'Step10_NicheAnalysis'))) {
        link_name <- paste0("download_", file)
        observeEvent(input[[link_name]], {
          selected_file(file.path(output.dir,'Step10_NicheAnalysis', file)) # '/Step5.Visualization/'
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
    #show figures
    img_dir <- file.path(output.dir,'Step10_NicheAnalysis')
    images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
    if(length(images)!=0){
      output$st_step10_plot<-renderSlickR({
        image_names <- basename(images)
        slickR(images,width = "80%",height = 100) %synch%
          (slickR(image_names, slideType = 'p') + settings(arrows = FALSE))
      })
    }else{
      output$st_step10_text<-renderText({
        print("NO Figure !")
      })
    }
  })
  
  # Step11. Generate the Report---------------------------------------------------------------------------------
  
  observeEvent(input$RunStep11_st, {
    
    if(is.null(input$st_jobid)){
      
      # load previous parameters
      output.dir <- output.dir.temp()
    }
    # load previous variables
    previous_results_path <- previous_results_path.temp()
    # load previous results
    Load_previous_results(previous_results_path = previous_results_path)
    
    # set parameters
    
    # convert some parameters
    
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
      saveRDS(var, file = paste0(previous_results_path, '/', var_name, ".rds"))  # Save as RDS with the variable's name
    }
    
    # update parameters
    
    shinyjs::enable('RunStep11')
    output$step10_completed <- renderText({'Step11 completed'})
  })
  observeEvent(input$back_button_st, {
    shinyjs::hide("ui3")
    shinyjs::show("ui1")
  })
  
}

  
  
  
  



















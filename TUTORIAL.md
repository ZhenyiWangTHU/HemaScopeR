# Usage instructions:

## 1. scRNA-seq pipeline

### A. Run the complete scRNA-seq pipeline automatically
```R
library(HemaScopeR)
scRNASeq_10x_pipeline(
                     # input and output
                     input.data.dirs = c('./HemaScopeR/test/caseData/GSE120221/1399/SRR7881399/SRR7881399/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1400/SRR7881400/SRR7881400/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1401/SRR7881401/SRR7881401/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1402/SRR7881402/SRR7881402/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1403/SRR7881403/SRR7881403/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1404/SRR7881404/SRR7881404/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1405/SRR7881405/SRR7881405/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1406/SRR7881406/SRR7881406/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1407/SRR7881407/SRR7881407/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1408/SRR7881408/SRR7881408/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1409/SRR7881409/SRR7881409/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1410/SRR7881410/SRR7881410/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1411/SRR7881411/SRR7881411/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1412/SRR7881412/SRR7881412/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1413/SRR7881413/SRR7881413/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1414/SRR7881414/SRR7881414/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1415/SRR7881415/SRR7881415/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1416/SRR7881416/SRR7881416/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1417/SRR7881417/SRR7881417/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1418/SRR7881418/SRR7881418/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1419/SRR7881419/SRR7881419/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1420/SRR7881420/SRR7881420/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1421/SRR7881421/SRR7881421/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1422/SRR7881422/SRR7881422/outs/filtered_feature_bc_matrix',
                                         './HemaScopeR/test/caseData/GSE120221/1423/SRR7881423/SRR7881423/outs/filtered_feature_bc_matrix'),
                     project.names = c(  'SRR7881399',
                                         'SRR7881400',
                                         'SRR7881401',
                                         'SRR7881402',
                                         'SRR7881403',
                                         'SRR7881404',
                                         'SRR7881405',
                                         'SRR7881406',
                                         'SRR7881407',
                                         'SRR7881408',
                                         'SRR7881409',
                                         'SRR7881410',
                                         'SRR7881411',
                                         'SRR7881412',
                                         'SRR7881413',
                                         'SRR7881414',
                                         'SRR7881415',
                                         'SRR7881416',
                                         'SRR7881417',
                                         'SRR7881418',
                                         'SRR7881419',
                                         'SRR7881420',
                                         'SRR7881421',
                                         'SRR7881422',
                                         'SRR7881423'), 
                     output.dir = './HemaScopeR/test/caseData/GSE120221/',
                     pythonPath = './anaconda3/envs/HemaScopeR/bin/python',
                     # quality control and preprocessing
                     gene.column = 2,
                     min.cells = 10,
                     min.feature = 200,
                     mt.pattern = '^MT-',
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
                     ViolinPlot.cellTypeOrders = as.character(1:22),
                     ViolinPlot.cellTypeColors = NULL,
                     Org = 'hsa',
                     loom.files.path = c( './HemaScopeR/test/caseData/GSE120221/1399/SRR7881399/SRR7881399/velocyto/SRR7881399.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1400/SRR7881400/SRR7881400/velocyto/SRR7881400.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1401/SRR7881401/SRR7881401/velocyto/SRR7881401.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1402/SRR7881402/SRR7881402/velocyto/SRR7881402.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1403/SRR7881403/SRR7881403/velocyto/SRR7881403.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1404/SRR7881404/SRR7881404/velocyto/SRR7881404.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1405/SRR7881405/SRR7881405/velocyto/SRR7881405.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1406/SRR7881406/SRR7881406/velocyto/SRR7881406.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1407/SRR7881407/SRR7881407/velocyto/SRR7881407.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1408/SRR7881408/SRR7881408/velocyto/SRR7881408.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1409/SRR7881409/SRR7881409/velocyto/SRR7881409.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1410/SRR7881410/SRR7881410/velocyto/SRR7881410.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1411/SRR7881411/SRR7881411/velocyto/SRR7881411.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1412/SRR7881412/SRR7881412/velocyto/SRR7881412.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1413/SRR7881413/SRR7881413/velocyto/SRR7881413.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1414/SRR7881414/SRR7881414/velocyto/SRR7881414.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1415/SRR7881415/SRR7881415/velocyto/SRR7881415.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1416/SRR7881416/SRR7881416/velocyto/SRR7881416.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1417/SRR7881417/SRR7881417/velocyto/SRR7881417.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1418/SRR7881418/SRR7881418/velocyto/SRR7881418.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1419/SRR7881419/SRR7881419/velocyto/SRR7881419.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1420/SRR7881420/SRR7881420/velocyto/SRR7881420.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1421/SRR7881421/SRR7881421/velocyto/SRR7881421.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1422/SRR7881422/SRR7881422/velocyto/SRR7881422.loom',
                                          './HemaScopeR/test/caseData/GSE120221/1423/SRR7881423/SRR7881423/velocyto/SRR7881423.loom'),
                     # cell cycle
                     cellcycleCutoff = NULL,
                     # cell chat
                     sorting = FALSE,
                     ncores = 10,
                     # Verbose = FALSE,
                     # activeEachStep
                     Whether_load_previous_results = FALSE,
                     Step1_Input_Data = TRUE,
                     Step1_Input_Data.type = 'cellranger-count',
                     Step2_Quality_Control = TRUE,
                     Step2_Quality_Control.RemoveBatches = FALSE,
                     Step2_Quality_Control.RemoveDoublets = TRUE,
                     Step3_Clustering = TRUE,
                     Step4_Identify_Cell_Types = TRUE,
                     Step4_Use_Which_Labels = 'clustering',
                     Step4_Cluster_Labels = NULL,
                     Step4_Changed_Labels = NULL,
                     Step4_run_sc_CNV = FALSE,
                     Step5_Visualization = TRUE,
                     Step6_Find_DEGs = TRUE,
                     Step7_Assign_Cell_Cycle = TRUE,
                     Step8_Calculate_Heterogeneity = TRUE,
                     Step9_Violin_Plot_for_Marker_Genes = TRUE,
                     Step10_Calculate_Lineage_Scores = TRUE,
                     Step11_GSVA = TRUE,
                     Step11_GSVA.identify.cellType.features=TRUE,
                     Step11_GSVA.identify.diff.features=FALSE,
                     Step11_GSVA.comparison.design=NULL,
                     Step12_Construct_Trajectories = TRUE,
                     Step12_Construct_Trajectories.clusters = c('11', '15', '9', '10', '6', '19', '3', '14'), 
                     Step12_Construct_Trajectories.monocle = TRUE,
                     Step12_Construct_Trajectories.slingshot = TRUE,
                     Step12_Construct_Trajectories.scVelo = TRUE,
                     Step13_TF_Analysis = TRUE,
                     Step14_Cell_Cell_Interaction = TRUE,
                     Step15_Generate_the_Report = TRUE
        )
```

### B. Run the scRNA-seq pipeline step by step

Follow the documents of each function in HemaScopeR.

## 2. st-seq pipeline

### A. Run the complete st-seq pipeline automatically

```R
st_10x_visium_pipeline(
    input.data.dir = 'Path/to/data',
    output.dir = 'Path/to/save',
    sampleName = 'Hema',
    Step2_QC = T,
    Step3_Clustering = T,
    Step4_Find_DEGs = T,
    Step5_SVFs = T,
    Step6_Interaction = T,
    Step7_CNV = T,
    Step8_Deconvolution = T,
    Step9_Cellcycle = T,
    Step10_Nich = T,
    
    # settings
    verbose = FALSE,
    species = 'mouse', # human or mosue
    genReport = TRUE
)
```

### B. Run the st-seq pipeline step by step

Follow the documents of each function in HemaScopeR.

## 3. Graphical user interface

### A. Start the shiny GUI

```R
library(HemaScopeR)
shinyApp(ui = ui, 
         server = server, 
         options = list(host = your_host, port = your_port))
```

### B. Follow the instructions on the GUI

## 4. Docker

### A. Pull the 'hemascoper' docker image

### B. Start a docker container

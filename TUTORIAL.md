# Usage instructions:

## 1. scRNA-seq pipeline

### A. Run the complete scRNA-seq pipeline automatically
```R
library(HemaScopeR)
scRNASeq_10x_pipeline(
                     # input and output
                     input.data.dirs = c('/input/path/dataset1',
                                         '/input/path/dataset2'),
                     project.names = c('project1',
                                       'project2'), 
                     output.dir = '/output/path/',
                     pythonPath = '/python/path/',
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
                     loom.files.path = c( '/input/velocyto/dataset1.loom',
                                          '/input/velocyto/dataset2.loom'),
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
                     Step12_Construct_Trajectories.clusters = c('cluster1', 'cluster2', 'cluster3'), 
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
library(HemaScopeR)
st_10x_visium_pipeline(
    input.data.dir = '/Path/to/data',
    output.dir = '/Path/to/save',
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

```shell
docker pull l1hj/hemascoper
```
### B. Start a docker container

```shell
docker run -it --security-opt seccomp=unconfined hemascoper /bin/bash
```

Jupyter notebook is also available, which can be started as follow (e.g. set the port as 8888):
```shell
docker run -i -t -p 8888:8888 continuumio/miniconda3 /bin/bash -c "\
    mkdir -p /opt/notebooks && \
    jupyter notebook \
    --notebook-dir=/opt/notebooks --ip='*' --port=8888 \
    --no-browser --allow-root"
```

## 5. Demo datasets
The demo datasets are available in our Cloud Drive via this link [https://cloud.tsinghua.edu.cn/d/3d363e32665249409571/](https://cloud.tsinghua.edu.cn/d/3d363e32665249409571/).

### A. scRNA-seq demo
We use the SRR7881414 dataset from demo datasets to demonstrate the scRNA-seq pipeline in HemaScopeR.
```R
library(HemaScopeR)
scRNASeq_10x_pipeline(
                     # input and output
                     input.data.dirs = c('./testData/SRR7881414_filtered_feature_bc_matrix'),
                     project.names = c('SRR7881414'), 
                     output.dir = './testData/SRR7881414_result',
                     pythonPath = '/python/path',
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
                     ViolinPlot.cellTypeOrders = NULL,
                     ViolinPlot.cellTypeColors = NULL,
                     Org = 'hsa',
                     loom.files.path = c('./testData/SRR7881414.loom'),
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
                     Step11_GSVA.comparison.design=NULL,
                     Step12_Construct_Trajectories = TRUE,
                     Step12_Construct_Trajectories.clusters = NULL, 
                     Step12_Construct_Trajectories.monocle = TRUE,
                     Step12_Construct_Trajectories.slingshot = TRUE,
                     Step12_Construct_Trajectories.scVelo = TRUE,
                     Step13_TF_Analysis = TRUE,
                     Step14_Cell_Cell_Interaction = TRUE,
                     Step15_Generate_the_Report = TRUE
        )
```

### B. st-seq demo
We use one sample from the GSE230207 dataset in the demo datasets to showcase the st-seq pipeline in HemaScopeR. The other st-seq datasets within the demo can be run using the same code. And you just need to adjust specific paths and parameters.
```R 
library(HemaScopeR)
st_10x_visium_pipeline(
    input.data.dir = './testData/GSE230207_RAW/hot',
    output.dir = './testData/GSE230207_hot_result',
    sampleName = 'GSE230207_hot',
    Step2_QC = TRUE,
    Step3_Clustering = TRUE,
    Step4_Find_DEGs = TRUE,
    Step5_SVFs = TRUE,
    Step6_Interaction = TRUE,
    Step7_CNV = TRUE,
    Step8_Deconvolution = TRUE,
    Step9_Cellcycle = TRUE,
    Step10_Nich = TRUE,
    
    # settings
    verbose = FALSE,
    species = 'human', # human or mosue
    genReport = TRUE
)
```

### C. Examples of HTML-formatted analysis reports
The HTML-formatted demo analysis reports can be accessed via these links [https://github.com/ZhenyiWangTHU/HemaScopeR/blob/main/sc_base.html](https://github.com/ZhenyiWangTHU/HemaScopeR/blob/main/sc_base.html), at the same time, we are continually updating and optimizing the format and content of these reports.
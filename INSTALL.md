# Installing HemaScopeR:

## Quick Installation

### 1. Install R, python and RStudio

**HemaScopeR** is a package written in R and Python, it also can be used in the RStudio interactive environment.

R can be obtained and installed from [https://cran.rstudio.com](https://cran.rstudio.com).

Python can be obtained and installed from [https://www.python.org/downloads/release/python-3114/](https://www.python.org/downloads/release/python-3114/).

Please install R upper 4.3.1 and python upper 3.11.4.

Following installation of R, Rstudio can be obtained and installed from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/). The free version of RStudio Desktop is sufficient.

### 2. Install required R packages and HemaScopeR

We wrote a script that will attempt to install all requirements for HemaScopeR and then HemaScopeR itself. Start RStudio and then run the installation script from the console:

```source("https://github.com/ZhenyiWangTHU/HemaScopeR/blob/master/HemaScopeR-install.R")```

### Manual installation

If something went wrong previously, you may have try installing some of HemaScopeR's dependencies manually:

###### A. Install and attach the *devtools* package

*devtools* is required to compile and install HemaScopeR, since it's distributed through GitHub.

```install.packages("devtools")```
     
###### B. Install required packages

Because these packages must be installed before installing HemaScopeR (otherwise its installation will fail.)

``` install.packages(setdiff(c('shiny','textshaping','shinyjs','Seurat','phateR','DoubletFinder','monocle','slingshot','URD','GSVA','limma','plyr','dplyr','org.Mm.eg.db','org.Hs.eg.db','CellChat','velocyto.R','SeuratWrappers','stringr','scran','ggpubr','viridis','pheatmap','parallel','reticulate','SCENIC','feather','AUCell','RcisTarget','Matrix','foreach','doParallel','clusterProfiler','RColorBrewer','Rfast2','SeuratDisk','abcCellmap','biomaRt','copykat','gelnet','ggplot2','parallelDist','patchwork','markdown'), .packages(all.available=T)))```
     
###### C. Install HemaScopeR

*HemaScopeR* can be installed directly from the GitHub repository

```library(devtools)```  
```install_github(repo="ZhenyiWangTHU/HemaScopeR",subdir = "/HemaScopeR")```

If you do not want to update the installed dependencies, you could run this command
```install_github(repo="ZhenyiWangTHU/HemaScopeR",subdir = "/HemaScopeR", dep = FALSE)```

Or download HemaScopeR_1.0.0.tar.gz and install in R

```install.packages('HemaScopeR_1.0.0.tar.gz')```

### 3. Pull Docker image 'HemaScopeR' from Docker Hub
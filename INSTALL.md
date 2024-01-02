# Installing HemaScopeR:

## Quick Installation

### 1. Install R, python and optionally RStudio

**HemaScopeR** is a package written in R and Python. It can also be utilized within the interactive environment of RStudio.

R can be obtained and installed from [https://cran.rstudio.com](https://cran.rstudio.com).

Python can be obtained and installed from [https://www.python.org](https://www.python.org).

Please install R version 4.2.3 and python version 3.11.4 or higher.

Following installation of R, Rstudio can be obtained and installed from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/). The free version of RStudio Desktop is sufficient.

### 2. Install required R packages and HemaScopeR

We wrote a script that will attempt to install all requirements for HemaScopeR and then HemaScopeR itself. Start RStudio and then run the installation script from the console:

```R
source("https://github.com/ZhenyiWangTHU/HemaScopeR/blob/master/HemaScopeR-install.R")
```

### Manual installation

If something went wrong previously, you may have try installing some of HemaScopeR's dependencies manually:

###### A. Install and attach the *devtools* package

*devtools* is required to compile and install HemaScopeR, since it's distributed through GitHub.

```R
install.packages("devtools")
```
     
###### B. Install required R packages

Installation of these packages is necessary prior to HemaScopeR installation; otherwise, the installation process will fail.

```R
c('shiny','textshaping','shinyjs','Seurat','phateR','DoubletFinder','monocle','slingshot','GSVA','limma','plyr','dplyr','org.Mm.eg.db','org.Hs.eg.db','CellChat','velocyto.R','SeuratWrappers','stringr','scran','ggpubr','viridis','pheatmap','parallel','reticulate','SCENIC','feather','AUCell','RcisTarget','Matrix','foreach','doParallel','clusterProfiler','OpenXGR','RColorBrewer','Rfast2','SeuratDisk','abcCellmap','biomaRt','copykat','gelnet','ggplot2','parallelDist','patchwork','markdown')
```

###### C. Install required Python packages

Because these packages must be installed for Python-based steps in HemaScopeR; otherwise, Python-based analysis steps cannot be performed.

```Python
numpy 1.23.5
pandas 1.3.5
scvelo 0.2.5
pickle 4.0
arboreto 0.1.6
ot 0.9.1
anndata 0.9.2
scanpy 1.9.4
scipy 1.11.2
seaborn 0.12.2
commot 0.0.3
matplotlib 3.7.2
cell2location 0.1.3
scvi 1.0.3
phate 1.0.11
```
     
###### D. Install HemaScopeR

*HemaScopeR* can be installed directly from the GitHub repository

```R
library(devtools)
```  

```R
install_github(repo="ZhenyiWangTHU/HemaScopeR",subdir = "/HemaScopeR")
```

If you do not want to update the installed dependencies, you could run this command

```R
install_github(repo="ZhenyiWangTHU/HemaScopeR",subdir = "/HemaScopeR", dep = FALSE)
```

Or download HemaScopeR_1.0.0.tar.gz and install in R

```R
install.packages('HemaScopeR_1.0.0.tar.gz')
```

### 3. Pull Docker image 'hemascoper' from Docker Hub

*HemaScopeR* can be accessed via Docker Hub at [https://hub.docker.com/r/l1hj/hemascoper](https://hub.docker.com/r/l1hj/hemascoper) or through Docker pull command 

```shell
docker pull l1hj/hemascoper
```

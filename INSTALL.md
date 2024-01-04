# Installing HemaScopeR:

## Quick Installation

### 1. Install R, python and optionally RStudio

**HemaScopeR** is a package written in R and Python. It can also be utilized within the interactive environment of RStudio.

R can be obtained and installed from [https://cran.rstudio.com](https://cran.rstudio.com).

Python can be obtained and installed from [https://www.python.org](https://www.python.org).

Please install R version 4.2.3 and python version 3.11.4 or higher.

Following installation of R, Rstudio can be obtained and installed from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/). The free version of RStudio Desktop is sufficient.

### 2. Install required R-packages, Python-packages and HemaScopeR

###### A. Install required R-packages and HemaScopeR automatically

We wrote a script that will attempt to install all required R-packages for HemaScopeR and then HemaScopeR itself. Start RStudio and then run the installation script from the console:

```R
source("https://github.com/ZhenyiWangTHU/HemaScopeR/blob/master/HemaScopeR-install.R")
```

If something went wrong previously, you may have try installing some of HemaScopeR's dependencies manually:

###### B. Install and attach the *devtools* package

*devtools* is required to compile and install HemaScopeR, since it's distributed through GitHub.

```R
install.packages("devtools")
```
     
###### C. Install required R-packages

Installation of these R-packages is necessary prior to HemaScopeR installation.

```R
shiny 1.7.5
textshaping 0.3.6
shinyjs 2.1.0
Seurat 4.3.0.1
phateR 1.0.7
DoubletFinder 2.0.3
monocle 2.28.0
slingshot 2.8.0
GSVA 1.48.3
limma 3.56.2
plyr 1.8.8
dplyr 1.1.2
org.Mm.eg.db 3.17.0
org.Hs.eg.db 3.17.0
CellChat 1.6.1
velocyto.R 0.6
SeuratWrappers 0.3.1
stringr 1.5.0
scran 1.28.2
ggpubr 0.6.0
viridis 0.6.4
pheatmap 1.0.12
parallel 4.3.1
reticulate 1.31
SCENIC 1.1.2.1
feather 0.3.5
AUCell 1.22.0
RcisTarget 1.20.0
Matrix 1.6.1
foreach 1.5.2
doParallel 1.0.17
clusterProfiler 4.8.3
OpenXGR 1.0.0
RColorBrewer 1.1.3
Rfast2 0.1.5.1
SeuratDisk 0.0.0.9020
abcCellmap 0.1.0
biomaRt 2.56.1
copykat 1.1.0
gelnet 1.2.1
ggplot2 3.4.4
parallelDist 0.2.6
patchwork 1.1.3
markdown 1.10
```

###### D. Install required Python-packages

These Python-packages must be installed for Python-based steps in HemaScopeR.

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
     
###### E. Install HemaScopeR

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

###### F. Download databases
In addition to the R-packages and Python-packages, you will also need to download the databases for HemaScopeR. The databases are available in our Cloud Drive via this link [https://cloud.tsinghua.edu.cn/d/759fd04333274d3f9946/](https://cloud.tsinghua.edu.cn/d/759fd04333274d3f9946/).

### 3. Pull Docker image 'hemascoper' from Docker Hub

*HemaScopeR* can be accessed via Docker Hub at [https://hub.docker.com/r/l1hj/hemascoper](https://hub.docker.com/r/l1hj/hemascoper) or through Docker pull command 

```shell
docker pull l1hj/hemascoper
```

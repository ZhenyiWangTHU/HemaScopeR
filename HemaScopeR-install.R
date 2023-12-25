# Install requirements: 'shiny','textshaping','shinyjs','Seurat','phateR','DoubletFinder','monocle','slingshot','GSVA','limma','plyr','dplyr','org.Mm.eg.db','org.Hs.eg.db','CellChat','velocyto.R','SeuratWrappers','stringr','scran','ggpubr','viridis','pheatmap','parallel','reticulate',
# 'SCENIC','feather','AUCell','GENIE3','RcisTarget','Matrix','foreach','doParallel','clusterProfiler','RColorBrewer','Rfast2','SeuratDisk','abcCellmap','biomaRt','copykat','gelnet','ggplot2','parallelDist','patchwork','markdown'
# Require R version >= 4.3.1

if (as.numeric(R.version$major) <= 4 && as.numeric(R.version$minor) < 3) {
  # R < 4.3
  message("Require R version >= 4.3.1")
} else {
  # R >= 4.3
  message("Installing required packages from cran.")
  install.packages(c('remotes', 'shiny', 'textshaping', 'shinyjs', 'Seurat', 'phateR', 'plyr', 'dplyr', 'stringr', 'ggpubr', 'viridis', 'pheatmap', 'reticulate', 'feather', 'Matrix', 'foreach', 'doParallel', 'RColorBrewer', 'Rfast2', 'gelnet', 'ggplot2', 'parallelDist', 'patchwork', 'markdown'))

  message("Installing required packages from Bioconductor.")
  if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
  BiocManager::install(c('monocle', 'slingshot', 'GSVA', 'limma', 'org.Mm.eg.db', 'org.Hs.eg.db', 'scran', 'AUCell', 'RcisTarget', 'GENIE3', 'clusterProfiler', 'biomaRt'))

  message("Installing required packages from GitHub.")
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::install_github("jinworks/CellChat")
    devtools::install_github("velocyto-team/velocyto.R")
    devtools::install_github("aertslab/SCENIC")
    devtools::install_github("pzhulab/abcCellmap")
    devtools::install_github("navinlabcode/copykat")
  } else {
    message("Installing devtools")
    install.packages("devtools")
    devtools::install_github("jinworks/CellChat")
    devtools::install_github("velocyto-team/velocyto.R")
    devtools::install_github("aertslab/SCENIC")
    devtools::install_github("pzhulab/abcCellmap")
    devtools::install_github("navinlabcode/copykat")
  }

  if (requireNamespace("remotes", quietly = TRUE)) {
    remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
    remotes::install_github('satijalab/seurat-wrappers')
    remotes::install_github("mojaveazure/seurat-disk")
  } else {
    message("Installing remotes")
    install.packages("remotes")
    remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
    remotes::install_github('satijalab/seurat-wrappers')
    remotes::install_github("mojaveazure/seurat-disk")
  }
}

# Check that packages installation went smoothly.
# if (!requireNamespace("parallel", quietly = TRUE)) {stop("Failed to install required package 'parallel'.")}

# Install HemaScopeR
  message("Installing HemaScopeR")
  devtools::install_github(repo="ZhenyiWangTHU/HemaScopeR",subdir = "/HemaScopeR")

# Check that HemaScopeR installed.
if (requireNamespace("HemaScopeR", quietly = TRUE)) {
  message("HemaScopeR installed successfully!")
  message('You can load it by typing: library("HemaScopeR")')
  message('Try "?HemaScopeR" for starting tips.')
} else {
  message("Something went wrong. It doesn't seem that HemaScopeR installed correctly.")
}
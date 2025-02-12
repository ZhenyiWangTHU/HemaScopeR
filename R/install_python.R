#### Global variables ####
.HemaScope_env <- new.env()
.HemaScope_env$pythonPath.sc <- file.path(reticulate::miniconda_path(),
                                          'envs/HemaScope_sc/bin/python3.9')
.HemaScope_env$pythonPath.ST <- file.path(reticulate::miniconda_path(),
                                          'envs/HemaScope_ST/bin/python3.11')
.HemaScope_env$pythonPath.stereo <- file.path(reticulate::miniconda_path(),
                                              'envs/HemaScope_stereo/bin/python3.8')

#' The path to the miniconda environment of scRNA-seq pipeline
#'
#' @export
python.path.sc <- function(){
    options(future.globals.maxSize=4000000000)
    if(!file.exists(.HemaScope_env$pythonPath.sc)){
        .HemaScope_env$pythonPath.sc <- file.path(reticulate::miniconda_path(),
                                                  'envs/HemaScope_sc/bin/python3')
    }
    if(!file.exists(.HemaScope_env$pythonPath.sc)){
        .HemaScope_env$pythonPath.sc <- file.path(reticulate::miniconda_path(),
                                                  'envs/HemaScope_sc/python.exe')
    }
    if(!file.exists(.HemaScope_env$pythonPath.sc)){
        stop('Please run init_miniconda function first.')
    }
    .HemaScope_env$pythonPath.sc
}

#' The path to the miniconda environment of ST pipeline
#'
#' @export
python.path.ST <- function(){
    options(future.globals.maxSize=4000000000)
    if(!file.exists(.HemaScope_env$pythonPath.ST)){
        .HemaScope_env$pythonPath.ST <- file.path(reticulate::miniconda_path(),
                                                  'envs/HemaScope_ST/bin/python3')
    }
    if(!file.exists(.HemaScope_env$pythonPath.ST)){
        .HemaScope_env$pythonPath.ST <- file.path(reticulate::miniconda_path(),
                                                  'envs/HemaScope_ST/python.exe')
    }
    if(!file.exists(.HemaScope_env$pythonPath.ST)){
        stop('Please run init_miniconda function first.')
    }
    .HemaScope_env$pythonPath.ST
}

#' The path to the miniconda environment of stereo-seq pipeline
#'
#' @export
python.path.stereo <- function(){
    options(future.globals.maxSize=4000000000)
    if(!file.exists(.HemaScope_env$pythonPath.stereo)){
        .HemaScope_env$pythonPath.stereo <- file.path(reticulate::miniconda_path(),
                                                      'envs/HemaScope_stereo/bin/python3')
    }
    if(!file.exists(.HemaScope_env$pythonPath.stereo)){
        .HemaScope_env$pythonPath.stereo <- file.path(reticulate::miniconda_path(),
                                                      'envs/HemaScope_stereo/python.exe')
    }
    if(!file.exists(.HemaScope_env$pythonPath.stereo)){
        stop('Please run init_miniconda_stereo function first.')
    }
    .HemaScope_env$pythonPath.stereo
}

#' Init miniconda for scRNA-seq pipeline and basic ST pipline
#'
#' @import reticulate
#'
#' @export
#'
init_miniconda <- function(){
    conda_path <- reticulate::miniconda_path()
    if(file.exists(conda_path)){
        cat('Miniconda is installed at: ', conda_path, '\n')
    }else{
        cat("Miniconda is not installed. We will first install miniconda.\n")
        reticulate::install_miniconda()
    }

    ## For scRNA-seq
    cat('Creating the environment for scRNA-seq pipeline...')
    reticulate::conda_create(envname='HemaScope_sc', python_version='3.9.12')
    reticulate::conda_install(envname='HemaScope_sc', packages='anndata', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='arboreto', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='karateclub', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='networkx', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='phate', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='pot', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='scanpy', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='scipy', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='seaborn', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='distributed', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='matplotlib==3.6.2', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='dask==2022.2.1', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='distributed==2022.2.1', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='scvelo==0.2.5', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='numpy==1.23.5', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_sc', packages='pandas==1.5.3', pip=TRUE)

    .HemaScope_env$pythonPath.sc <- file.path(reticulate::miniconda_path(),
                                              'envs/HemaScope_sc/bin/python3.9')

    ## For ST
    cat('Creating the environment for ST pipeline...')
    reticulate::conda_create(envname='HemaScope_ST', python_version='3.11')
    reticulate::conda_install(envname='HemaScope_ST', packages='torch', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='torchvision', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='torchaudio', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='cell2location==0.1.3', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='scvi-tools==1.0.4', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='jax==0.4.14', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='jaxlib==0.4.14', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='orbax-checkpoint==0.4.1', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='flax==0.7.4', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='scipy==1.11.3', pip=TRUE)
    reticulate::conda_install(envname='HemaScope_ST', packages='commot', pip=TRUE)

    .HemaScope_env$pythonPath.ST <- file.path(reticulate::miniconda_path(),
                                              'envs/HemaScope_ST/bin/python3.11')
}

#' Init miniconda for stereoseq data
#'
#' @import reticulate
#'
#' @export
#'
init_miniconda_stereo <- function(){
    conda_path <- reticulate::miniconda_path()
    if(file.exists(conda_path)){
        cat('Miniconda is installed at: ', conda_path, '\n')
    }else{
        cat("Miniconda is not installed. We will first install miniconda.\n")
        reticulate::install_miniconda()
    }

    cat('Creating the environment for stereo-seq data.')

    reticulate::conda_create(envname='HemaScope_stereo', python_version='3.8')
    reticulate::conda_install(envname='HemaScope_stereo', packages='stereopy', pip=TRUE)

    .HemaScope_env$pythonPath.stereo <- file.path(reticulate::miniconda_path(),
                                                  'envs/HemaScope_stereo/bin/python3.8')
}




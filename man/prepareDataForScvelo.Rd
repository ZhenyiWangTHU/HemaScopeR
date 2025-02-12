% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_Step12_prepareDataForScvelo.R
\name{prepareDataForScvelo}
\alias{prepareDataForScvelo}
\title{Prepare Data for scVelo Analysis}
\usage{
prepareDataForScvelo(
  sc_object = NULL,
  loom.files.path = NULL,
  scvelo.reduction = "pca",
  scvelo.column = NULL,
  output.dir = NULL
)
}
\arguments{
\item{sc_object}{A Seurat object containing the single-cell RNA-seq data.}

\item{loom.files.path}{A character vector specifying the path(s) to the loom files for scVelo analysis.}

\item{scvelo.reduction}{A character specifying the reduction method used for scVelo analysis (default is 'pca').}

\item{scvelo.column}{A character specifying the column in the Seurat object metadata containing cell types.}

\item{output.dir}{A character specifying the directory where the scVelo input files will be saved.}
}
\value{
None. This function saves the required input files for scVelo analysis in the output.dir.
}
\description{
Function \code{prepareDataForScvelo} prepares the data for scVelo analysis by formatting and saving the necessary input files.
}
\details{
This function takes a Seurat object, cell embeddings, and cell type information as input and performs the following steps:

1. Formats the Seurat object metadata, including cell types and sample names.
2. Extracts the spliced, unspliced, and ambiguous count matrices from the Seurat object.
3. Combines the metadata and cell embeddings.
4. Writes the necessary input files for scVelo analysis, including cell embeddings, count matrices, and metadata.
}
\author{
Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
}

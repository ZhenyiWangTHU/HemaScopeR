% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_Step11_run_GSVA.R
\name{run_GSVA}
\alias{run_GSVA}
\title{Run GSVA Analysis}
\usage{
run_GSVA(
  sc_object = NULL,
  GSVA.genelist = NULL,
  GSVA.cellTypes = NULL,
  GSVA.cellTypes.orders = NULL,
  GSVA.cellGroups = NULL,
  GSVA.identify.cellType.features = TRUE,
  GSVA.identify.diff.features = FALSE,
  GSVA.comparison.design = NULL,
  OrgDB = NULL,
  output.dir = NULL
)
}
\arguments{
\item{sc_object}{A Seurat object containing the single-cell RNA-seq data.}

\item{GSVA.genelist}{A list of gene sets for GSVA analysis.}

\item{GSVA.cellTypes}{A character vector specifying the cell types or labels for each cell.}

\item{GSVA.cellTypes.orders}{A character vector specifying the order of cell types for visualization.}

\item{GSVA.cellGroups}{A character vector specifying the cell groups or conditions for each cell.}

\item{GSVA.identify.cellType.features}{Logical. If TRUE, identify cell type-specific features.}

\item{GSVA.identify.diff.features}{Logical. If TRUE, identify differentially expressed features between cell groups.}

\item{GSVA.comparison.design}{A list specifying the experimental design for differential GSVA analysis.}

\item{OrgDB}{An organism-specific annotation database (OrgDb) for gene symbol conversion. e.g. org.Mm.eg.db or org.Hs.eg.db.}

\item{output.dir}{The directory where the GSVA results and output files will be saved.}
}
\value{
None. This function saves GSVA results and output files in the specified output directory.
}
\description{
Function \code{run_GSVA} performs GSVA (Gene Set Variation Analysis) analysis on a given single-cell RNA-seq dataset.
}
\details{
This function runs GSVA analysis, which calculates enrichment scores for gene sets in each cell using the provided gene list.
It also performs differential GSVA analysis between specified cell groups and generates heatmaps of the results. GSVA was referred to
https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf;http://revigo.irb.hr/;https://www.slideshare.net/dgrapov/proteomics-workshop-2014-lab-dmitry-grapov
}
\author{
Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
}

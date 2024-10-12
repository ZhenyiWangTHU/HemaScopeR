#' Loading stereo-seq data
#'
#' @param input.data.dirs A character of data path
#' @param output.dir A character of path to store the h5ad file
#' @param data_type 'gem', 'gef' or 'h5ad'
#' @param sep Separator string for 'gem' format
#' @param bin_type The type of bin, if file format is stereo-seq file. 'bins' or 'ell_bins'.
#' @param bin_size the size of bin to merge. The parameter only takes effect when the value of `bin_type` is 'bins'.
#' @param is_sparse The matrix is sparse matrix if is_sparse is True else np.ndarray
#' @param gene_list Restrict to this gene list
#' @param region Restrict to this region, [minX, maxX, minY, maxY]
#' @param pythonPath The path to the Python environment to use for the analysis
#'
#' @details
#' This function loads stereo-seq data in different formats
#'
#' @import Seurat
#' @import SeuratDisk
#'
#' @export
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
stereo_Loading_Data <- function(
        input.data.dir, 
        output.dir = '.',
        sampleName = 'Hema_ST',
        assay = 'Spatial',
        data_type = 'gem',
        sep = '\t',
        bin_type = 'bins',
        bin_size = 100,
        spot_diameter = 80,
        is_sparse = TRUE,
        gene_list = NULL, 
        region = NULL,
        pythonPath = NULL
){
    if(!dir.exists(output.dir)){
        dir.create(output.dir)
    }
    
    if(is.null(pythonPath)==FALSE){ reticulate::use_python(pythonPath) }else{stop('Please set the path of Python.')}
    
    # use_condaenv(condaenv)
    source_python(file.path(system.file(package = "HemaScopeR"),
                            "python/stereo.py"))
    
    read_stereo(
        data_path = input.data.dir, 
        save_path = output.dir,
        data_type = data_type,
        sep = sep,
        bin_type = bin_type,
        bin_size = bin_size,
        is_sparse = is_sparse,
        gene_list = gene_list, 
        region = region
    )
    
    st_obj <- anndata2rds(
        h5ad.path = file.path(output.dir, 'stereo_adata.h5ad'), 
        metadata.path = file.path(output.dir, 'stereo_metadata.csv'), 
        spot_diameter = spot_diameter,
        assay = assay
    )
    
    st_obj@project.name <- sampleName
    st_obj$orig.ident <- factor(sampleName)
    st_obj@active.ident <- st_obj$orig.ident
    
    saveRDS(st_obj, file.path(output.dir, 'stereo.rds'))
    
    Rds2H5(object = st_obj,
           file.dir = output.dir,
           assay = assay,
           slice = 'slice1')
    
    return(st_obj)
}


#' Convert anndata to rds
#' 
#' @import Seurat
#' @import SeuratDisk
#' @import rjson
#' @import dplyr
#' 
#' @details
#' https://github.com/STOmics/stereopy/blob/dev/docs/source/_static/annh5ad2rds.R
anndata2rds <- function(
        h5ad.path, 
        metadata.path, 
        spot_diameter = 80,
        assay = 'Spatial'
){
    suppressMessages(suppressWarnings(
        Convert(h5ad.path, dest = "h5seurat", assay = assay, overwrite = TRUE)
    ))
    h5file <- paste(paste(unlist(strsplit(h5ad.path, "h5ad", fixed = TRUE)), collapse='h5ad'), "h5seurat", sep="")
    
    suppressMessages(suppressWarnings(
        object <- LoadH5Seurat(h5file, meta.data = FALSE)
    ))
    metadata <- read.csv(
        metadata.path,
        row.names = 1
    )
    
    object <- object[, rownames(metadata)]
    object@meta.data <- metadata
    
    if(packageVersion('SeuratObject') >= '5.0.0'){
        counts <- GetAssayData(object, layer = 'data',
                               assay = assay)
    }else{
        counts <- GetAssayData(object, slot = 'data',
                               assay = assay)
    }
    object@meta.data$nCount_Spatial <- colSums(as.matrix(counts))
    object@meta.data$nFeature_Spatial <- colSums(as.matrix(counts > 0))
    
    if (!is.null(object@reductions$spatial)) {
        object@reductions$spatial <- NULL
    }
    
    cell_coords <- unique(metadata[, c('x', 'y')])
    cell_coords['cell'] <- row.names(cell_coords)
    cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
    cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1
    
    # object of images$slice1@image, all illustrated as 1 since no concrete pic
    tissue_lowres_image <- matrix(1, max(cell_coords$y) * 0.04 + 1, max(cell_coords$x) * 0.04 + 1)
    
    # object of images$slice1@coordinates, concrete coordinate of X and Y
    tissue_positions_list <- data.frame(row.names = cell_coords$cell,
                                        tissue = 1,
                                        row = cell_coords$y, col = cell_coords$x,
                                        imagerow = cell_coords$y, imagecol = cell_coords$x)
    # @images$slice1@scale.factors
    scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 177, tissue_hires_scalef = 0.04, 
                                     tissue_lowres_scalef = 0.04, spot_diameter_fullres = spot_diameter))
    
    # generate object @images$slice1
    generate_BGI_spatial <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
        if (filter.matrix) {
            tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
        }
        unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
        spot.radius <- unnormalized.radius / max(dim(x = image))
        return(new(Class = 'VisiumV1',
                   image = image,
                   scale.factors = scalefactors(spot = scale.factors$spot_diameter_fullres,
                                                fiducial = scale.factors$fiducial_diameter_fullres,
                                                hires = scale.factors$tissue_hires_scalef,
                                                lowres = scale.factors$tissue_lowres_scalef),
                   coordinates = tissue.positions,
                   spot.radius = spot.radius))
    }
    
    BGI_spatial <- generate_BGI_spatial(image = tissue_lowres_image,
                                        scale.factors = fromJSON(scalefactors_json),
                                        tissue.positions = tissue_positions_list)
    
    # can be thought of as a background of spatial
    # import image into seurat object
    object@images[['slice1']] <- BGI_spatial
    object@images$slice1@key <- "slice1_"
    object@images$slice1@assay <- assay
    
    return(object)
}


#' Loading MERFISH data
#'
#' @param input.data.dirs A character of data path
#' @param output.dir A character of path to store the h5ad file
#' @param tech 'Vizgen', 'Xenium', 'Nonastring', or 'Akoya'
#' @param fov A character of the name of the sample
#' @param rds.file Whether it is a rds file
#' @param assay Name of assay
#' @param Vizgen.z Z-index to load; must be between 0 and 6, inclusive
#' @param Akoya.type 'processor', 'inform', or 'qupath'
#'
#'
#' @details
#' This function loads MERFISH data
#'
#' @import Seurat
#'
#' @export
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
MERFISH_Loading_Data <- function(
        input.data.dir, 
        output.dir = '.',
        tech = 'Vizgen',
        fov = 'fov',
        
        rds.file = FALSE,
        assay = NULL,
        Vizgen.z = 3L, 
        Akoya.type = 'inform'
){
    if(!dir.exists(output.dir)){
        dir.create(output.dir)
    }
    if(is.null(assay)){
        assay = tech
    }
    
    if(rds.file){
        st_obj <- readRDS(input.data.dir)
    }else{
        if(tech == 'Vizgen'){
            st_obj <- LoadVizgen(
                data.dir = input.data.dir, 
                fov = fov, 
                assay = assay, 
                z = Vizgen.z
            )
        }else if(tech == 'Xenium'){
            st_obj <- LoadXenium(
                data.dir = input.data.dir, 
                fov = fov, 
                assay = assay
            )
        }else if(tech == 'Nonastring'){
            st_obj <- LoadNanostring(
                data.dir = input.data.dir, 
                fov = fov, 
                assay = assay
            )
        }else if(tech == 'Akoya'){
            st_obj <- LoadAkoya(
                filename = input.data.dir, 
                type = Akoya.type, 
                fov = fov, 
                assay = assay
            )
        }else{
            stop(paste0(tech, ' is not supported now.'))
        }
    }
    
    saveRDS(st_obj, 
            file.path(output.dir, 'merfish.rds'))
    
    return(st_obj)
}

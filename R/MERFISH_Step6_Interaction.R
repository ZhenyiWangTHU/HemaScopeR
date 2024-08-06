#' Run spatial interaction
#'
#' @param h5ad_path The path of MERFISH data in h5ad format
#' @param counts_path The path of counts data in csv format
#' @param counts_transpose Whether to transpose the counts matrix. If rows stand for genes, it will be `FALSE`.
#' @param coordinates_path The path of coordinates data in xlsx format
#' @param coordinates_index_col The column index of coordinates file
#' @param metadata_path The path of metadata.csv
#' @param label_key The name of the label to determine the clusters
#' @param save_path The path to save output files
#' @param species The species of the sample, 'human' or 'mouse'
#' @param signaling_type The parameter of `Commot` to determine the type of interaction.
#' Choose from 'Secreted Signaling', 'Cell-Cell Contact', and 'ECM-Receptor' for CellChatDB or
#' 'Secreted Signaling' and 'Cell-Cell Contact' for CellPhoneDB_v4.0. If None, all pairs in the database are returned.
#' @param database The parameter of `Commot` to determine the database used. 'CellChat' or 'CellPhoneDB_v4.0'
#' @param min_cell_pct The parameter of `Commot`. The minimum expression percentage required for LR pairs to be kept
#' @param dis_thr The parameter of `Commot`. The threshold of spatial distance of signaling
#' @param n_permutations The parameter of `Commot`. Number of label permutations for computing the p-value
#' @param pythonPath The path to the Python environment to use for the analysis.
#'
#'
#' @import reticulate
#'
#' @export
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#'
MERFISH_Interaction <- function(
        h5ad_path = NULL,
        counts_path = NULL,
        coordinates_path = NULL,
        coordinates_index_col = 0,
        counts_transpose = TRUE,
        metadata_path = NULL,
        library_id = 'Hema_ST',
        label_key = 'seurat_clusters',
        save_path = '.',
        species = 'mouse',
        signaling_type = 'Secreted Signaling',
        database = 'CellChat',
        min_cell_pct = 0.05,
        dis_thr = 500,
        n_permutations = 100,
        pythonPath = NULL
        # condaenv = 'r-reticulate'
){
    if(!dir.exists(save_path)){
        dir.create(save_path)
    }
    if(!dir.exists(file.path(save_path, 'Pathway_matrix'))){
        dir.create(file.path(save_path, 'Pathway_matrix'))
    }
    if(!dir.exists(file.path(save_path, 'Pathway_pvalue'))){
        dir.create(file.path(save_path, 'Pathway_pvalue'))
    }
    if(!dir.exists(file.path(save_path, 'LR_matrix'))){
        dir.create(file.path(save_path, 'LR_matrix'))
    }
    if(!dir.exists(file.path(save_path, 'LR_pvalue'))){
        dir.create(file.path(save_path, 'LR_pvalue'))
    }
    
    if(is.null(pythonPath)==FALSE){ reticulate::use_python(pythonPath) }else{stop('Please set the path of Python.')}
    
    # use_condaenv(condaenv)
    source_python(file.path(system.file(package = "HemaScopeR"),
                            "python/Commot_MERFISH.py"))
    
    if(!is.null(h5ad_path)){
        print(paste0('The h5ad file ', h5ad_path, ' will be loaded.'))
    }else{
        if(is.null(counts_path)){
            stop('Please provide the path to the counts file.')
        }
        if(is.null(coordinates_path)){
            stop('Please provide the path to the coordinate file.')
        }
    }
    
    run_Commot_MERFISH(
        h5ad_path = h5ad_path,
        counts_path = counts_path,
        coordinates_path = coordinates_path,
        index_col = coordinates_index_col,
        counts_transpose = counts_transpose,
        metadata_path = metadata_path,
        library_id = library_id,
        label_key = label_key,
        save_path = save_path,
        species = species,
        signaling_type = signaling_type,
        database = database,
        min_cell_pct = as.numeric(min_cell_pct),
        dis_thr = as.integer(dis_thr),
        n_permutations = as.integer(n_permutations)
    )
    
    #### Interaction within clusters ####
    ## Show pathways
    commot_matrix.pathways = as.data.frame(read.csv(file.path(save_path, "Pathway_within_cluster_matrix.csv"),
                                                    row.names = 1))
    commot_pvalue.pathways = as.data.frame(read.csv(file.path(save_path, "Pathway_within_cluster_pvalue.csv"),
                                                    row.names = 1))
    plotStrengthAndPvalue_within(commot_matrix = commot_matrix.pathways,
                                 commot_pvalue = commot_pvalue.pathways,
                                 save_path = save_path,
                                 figure_name = 'Pathway_interaction_within_Dotplot',
                                 height = 1.0 + ncol(commot_matrix.pathways) * 0.11,
                                 width = 2.5 + nrow(commot_matrix.pathways) * 0.2)
    
    ## Show LRs
    commot_matrix.LR = as.data.frame(read.csv(file.path(save_path, "LR_within_cluster_matrix.csv"),
                                              row.names = 1, check.names = F))
    commot_pvalue.LR = as.data.frame(read.csv(file.path(save_path, "LR_within_cluster_pvalue.csv"),
                                              row.names = 1, check.names = F))
    
    # filter pathways according to the strength and pvalue
    tmp <- colSums(commot_matrix.LR)
    pathways.strength <- names(tmp[tmp > 0])
    tmp <- colSums(commot_pvalue.LR <= 0.05)
    pathways.pvalue <- names(tmp[tmp > 0])
    pathways.filter <- intersect(pathways.strength, pathways.pvalue)
    commot_matrix.LR <- commot_matrix.LR[, pathways.filter]
    commot_pvalue.LR <- commot_pvalue.LR[, pathways.pvalue]
    
    plotStrengthAndPvalue_within(commot_matrix = commot_matrix.LR,
                                 commot_pvalue = commot_pvalue.LR,
                                 save_path = save_path,
                                 figure_name = 'LR_interaction_within_Dotplot',
                                 height = 1.0 + ncol(commot_matrix.LR) * 0.11,
                                 width = 4.0 + nrow(commot_matrix.LR) * 0.2)
    
    #### Interaction between clusters ####
    pathways.matrix <- gsub('.csv', '',
                            list.files(file.path(save_path, 'Pathway_matrix')))
    pathways.pvalue <- gsub('.csv', '',
                            list.files(file.path(save_path, 'Pathway_pvalue')))
    pathways <- intersect(pathways.matrix, pathways.pvalue)
    
    commot_matrix.pathway.list <- list()
    commot_pvalue.pathway.list <- list()
    
    for(pathway in pathways){
        commot_matrix <- read.csv(file.path(save_path, 'Pathway_matrix',
                                            paste0(pathway, '.csv')),
                                  row.names = 1, check.names = F)
        commot_pvalue <- read.csv(file.path(save_path, 'Pathway_pvalue',
                                            paste0(pathway, '.csv')),
                                  row.names = 1, check.names = F)
        
        commot_matrix.pathway.list[[pathway]] <- commot_matrix
        commot_pvalue.pathway.list[[pathway]] <- commot_pvalue
    }
    
    plotStrengthAndPvalue_total(commot_matrix.list = commot_matrix.pathway.list,
                                commot_pvalue.list = commot_pvalue.pathway.list,
                                save_path = save_path,
                                figure_name = 'Pathway_interaction_between_Dotplot',
                                between.cluster = T,
                                height = 1.0 + length(commot_matrix.pathway.list) * 0.11,
                                width = 2.5 + 0.2 * nrow(commot_matrix.pathway.list[[1]]) *
                                    (ncol(commot_matrix.pathway.list[[1]]) - 1))
    
    LRs.matrix <- gsub('.csv', '',
                       list.files(file.path(save_path, 'LR_matrix')))
    LRs.pvalue <- gsub('.csv', '',
                       list.files(file.path(save_path, 'LR_pvalue')))
    LRs <- intersect(LRs.matrix, LRs.pvalue)
    
    commot_matrix.LR.list <- list()
    commot_pvalue.LR.list <- list()
    
    for(LR in LRs){
        commot_matrix <- read.csv(file.path(save_path, 'LR_matrix',
                                            paste0(LR, '.csv')),
                                  row.names = 1, check.names = F)
        commot_pvalue <- read.csv(file.path(save_path, 'LR_pvalue',
                                            paste0(LR, '.csv')),
                                  row.names = 1, check.names = F)
        
        commot_matrix.LR.list[[LR]] <- commot_matrix
        commot_pvalue.LR.list[[LR]] <- commot_pvalue
    }
    
    plotStrengthAndPvalue_total(commot_matrix.list = commot_matrix.LR.list,
                                commot_pvalue.list = commot_pvalue.LR.list,
                                save_path = save_path,
                                figure_name = 'LR_interaction_between_Dotplot',
                                between.cluster = T,
                                height = 1.0 + length(commot_matrix.LR.list) * 0.11,
                                width = 4.0 + 0.2 * nrow(commot_matrix.LR.list[[1]]) *
                                    (ncol(commot_matrix.LR.list[[1]]) - 1))
    
}

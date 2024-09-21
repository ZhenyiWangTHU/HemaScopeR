#' Run spatial interaction
#'
#' @param st_data_path The path of original ST data
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
st_Interaction <- function(
        st_data_path,
        metadata_path,
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
                            "python/Commot.py"))

    run_Commot(
        st_data_path = st_data_path,
        metadata_path = metadata_path,
        # library_id = library_id,
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

#' Plot the strength and P-value
#'
#' @param commot_matrix The output file of `Commot`
#' @param commot_pvalue The output file of `Commot`
#' @param save_path The path to save outputs
#' @param figure_name The name of figures
#'
#' @import ggplot2
#' @import RColorBrewer
#'
plotStrengthAndPvalue_within <- function(
        commot_matrix,
        commot_pvalue,
        save_path = NULL,
        figure_name = 'Dotplot',
        height = 4,
        width = 4
){
    df <- data.frame(interaction_term = matrix(0, nrow = ncol(commot_matrix) * nrow(commot_matrix)),
                     pathway_name = matrix(0, nrow = ncol(commot_matrix) * nrow(commot_matrix)),
                     strength = matrix(0, nrow = ncol(commot_matrix) * nrow(commot_matrix)),
                     pvalue = matrix(0, nrow = ncol(commot_matrix) * nrow(commot_matrix)))

    ii = 1
    for(i in 1:nrow(commot_matrix)){
        for(j in 1:ncol(commot_matrix)){
            df$interaction_term[ii] <- paste0(rownames(commot_matrix)[i], '->',
                                              rownames(commot_matrix)[i])
            df$pathway_name[ii] <- paste0(colnames(commot_matrix[j]))
            df$strength[ii] <- commot_matrix[[i, j]]
            df$pvalue[ii] <- commot_pvalue[[i, j]]

            ii = ii + 1
        }
    }

    df$strength[df$strength > quantile(df$strength, 0.95, na.rm = T)] <- as.numeric(quantile(df$strength, 0.95, na.rm = T))
    df$strength[df$strength < quantile(df$strength, 0.05, na.rm = T)] <- as.numeric(quantile(df$strength, 0.05, na.rm = T))

    df$pvalue[df$pvalue > 0.05] = 1
    df$pvalue[df$pvalue > 0.01 & df$pvalue <= 0.05] = 2
    df$pvalue[df$pvalue <= 0.01] = 3
    values <- c(1, 2, 3)
    names(values) <- c("p > 0.05", "0.01 < p < 0.05",
                       "p < 0.01")

    suppressMessages(suppressWarnings(
        p <- ggplot(df, aes(x = interaction_term, y = pathway_name,
                            color = strength, size = pvalue)) + geom_point(pch = 16) +
            theme_linedraw() + theme(panel.grid.major = element_blank()) +
            scale_radius(range = c(min(df$pvalue), max(df$pvalue)),
                         breaks = sort(unique(df$pvalue)),
                         labels = names(values)[values %in% sort(unique(df$pvalue))], name = "p-value") +
            scale_colour_gradientn(colors = rev(brewer.pal(n = 99, name = "Spectral")),
                                   na.value = "white", limits = c(quantile(df$strength,
                                                                           0, na.rm = T),
                                                                  quantile(df$strength, 1, na.rm = T)),
                                   breaks = c(quantile(df$strength, 0, na.rm = T),
                                              quantile(df$strength, 1, na.rm = T)),
                                   labels = c("min", "max")) +
            guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Strength")) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                             vjust = 0.5), axis.title.x = element_blank(),
                  axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
    ))

    if(!is.null(save_path)){
        saveImage(
            output.dir = save_path,
            p = p,
            saveName = figure_name,
            height = height,
            width = width
        )
    }else{
        p
    }
}


#' Plot total strength and P-value
#'
#' @param commot_matrix The output file of `Commot`
#' @param commot_pvalue The output file of `Commot`
#' @param save_path The path to save outputs
#' @param figure_name The name of figures
#'
#' @import ggplot2
#' @import RColorBrewer
#'
plotStrengthAndPvalue_total <- function(
        commot_matrix.list,
        commot_pvalue.list,
        save_path = NULL,
        figure_name = 'Dotplot',
        between.cluster = TRUE,
        height = 4,
        width = 4
){
    if(between.cluster){
        df <- data.frame(interaction_term = matrix(0, nrow = length(commot_matrix.list) *
                                                       nrow(commot_matrix.list[[1]]) *
                                                       (ncol(commot_matrix.list[[1]]) - 1)),
                         pathway_name = matrix(0, nrow = length(commot_matrix.list) *
                                                   nrow(commot_matrix.list[[1]]) *
                                                   (ncol(commot_matrix.list[[1]]) - 1)),
                         strength = matrix(0, nrow = length(commot_matrix.list) *
                                               nrow(commot_matrix.list[[1]]) *
                                               (ncol(commot_matrix.list[[1]]) - 1)),
                         pvalue = matrix(0, nrow = length(commot_matrix.list) *
                                             nrow(commot_matrix.list[[1]]) *
                                             (ncol(commot_matrix.list[[1]]) - 1)))
    }else{
        df <- data.frame(interaction_term = matrix(0, nrow = length(commot_matrix.list) *
                                                       nrow(commot_matrix.list[[1]]) *
                                                       ncol(commot_matrix.list[[1]])),
                         pathway_name = matrix(0, nrow = length(commot_matrix.list) *
                                                   nrow(commot_matrix.list[[1]]) *
                                                   ncol(commot_matrix.list[[1]])),
                         strength = matrix(0, nrow = length(commot_matrix.list) *
                                               nrow(commot_matrix.list[[1]]) *
                                               ncol(commot_matrix.list[[1]])),
                         pvalue = matrix(0, nrow = length(commot_matrix.list) *
                                             nrow(commot_matrix.list[[1]]) *
                                             ncol(commot_matrix.list[[1]])))
    }

    ii = 1
    for(i in 1:length(commot_matrix.list)){
        commot_matrix <- commot_matrix.list[[i]]
        commot_pvalue <- commot_pvalue.list[[i]]

        for(j in 1:nrow(commot_matrix)){
            for(k in 1:ncol(commot_matrix)){

                if(between.cluster){
                    if(j == k){
                        next
                    }
                }

                df$interaction_term[ii] <- paste0(rownames(commot_matrix)[j], '->',
                                                  colnames(commot_matrix)[k])
                df$pathway_name[ii] <- names(commot_matrix.list)[i]
                df$strength[ii] <- commot_matrix[[j, k]]
                df$pvalue[ii] <- commot_pvalue[[j, k]]

                ii <- ii + 1
            }
        }
    }

    df$strength[df$strength > quantile(df$strength, 0.95, na.rm = T)] <- as.numeric(quantile(df$strength, 0.95, na.rm = T))
    df$strength[df$strength < quantile(df$strength, 0.05, na.rm = T)] <- as.numeric(quantile(df$strength, 0.05, na.rm = T))

    df$pvalue[df$pvalue > 0.05] = 1
    df$pvalue[df$pvalue > 0.01 & df$pvalue <= 0.05] = 2
    df$pvalue[df$pvalue <= 0.01] = 3
    values <- c(1, 2, 3)
    names(values) <- c("p > 0.05", "0.01 < p < 0.05",
                       "p < 0.01")

    suppressMessages(suppressWarnings(
        p <- ggplot(df, aes(x = interaction_term, y = pathway_name,
                            color = strength, size = pvalue)) + geom_point(pch = 16) +
            theme_linedraw() + theme(panel.grid.major = element_blank()) +
            scale_radius(range = c(min(df$pvalue), max(df$pvalue)),
                         breaks = sort(unique(df$pvalue)),
                         labels = names(values)[values %in% sort(unique(df$pvalue))], name = "p-value") +
            scale_colour_gradientn(colors = rev(brewer.pal(n = 99, name = "Spectral")),
                                   na.value = "white", limits = c(quantile(df$strength,
                                                                           0, na.rm = T),
                                                                  quantile(df$strength, 1, na.rm = T)),
                                   breaks = c(quantile(df$strength, 0, na.rm = T),
                                              quantile(df$strength, 1, na.rm = T)),
                                   labels = c("min", "max")) +
            guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Strength")) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                             vjust = 0.5), axis.title.x = element_blank(),
                  axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
    ))

    if(!is.null(save_path)){
        saveImage(
            output.dir = save_path,
            p = p,
            saveName = figure_name,
            height = height,
            width = width
        )
    }else{
        p
    }
}

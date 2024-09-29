#' The pipeline for analyzing stereo-seq data.
#'
#' Function \code{st_stereo_pipeline} encompasses the complete downstream
#' analysis pipeline for four kinds of MERFISH data within HemascopeR.
#' @param input.data.dirs A character of data path.
#' @param output.dir A character of path to store the results and figures
#' @param data_type 'gem', 'gef' or 'h5ad'
#' @param sep Separator string for 'gem' format
#' @param bin_type The type of bin, if file format is stereo-seq file. 'bins' or 'ell_bins'.
#' @param bin_size the size of bin to merge. The parameter only takes effect when the value of `bin_type` is 'bins'.
#' @param is_sparse The matrix is sparse matrix if is_sparse is True else np.ndarray
#' @param gene_list Restrict to this gene list
#' @param region Restrict to this region, [minX, maxX, minY, maxY]
#'
#' @param min.gene An integer representing the minimum number of genes detected in a spot
#' @param max.gene An integer representing the maximum number of genes detected in a spot
#' @param min.nUMI An integer representing the minimum number of nUMI detected in a spot
#' @param max.nUMI An integer representing the maximum number of nUMI detected in a spot
#' @param min.spot An integer representing the minimum number of spots expressing each gene
#' @param species A character representing the species of sample, `human` or `mouse`
#' @param bool.remove.mito A bool value indicating whether removing mitochondrial gene
#'
#' @param normalization.method 'SCTransform' or other choice of `normalization.method` in `NormalizaData`
#' @param n.dim.used An integer representing the number of dimension used for finding neighbors
#' @param resolution An float parameter of `FindClusters` of `Seurat`
#'
#' @param only.pos A bool value to indicate whether only return positive markers.
#' @param min.pct The parameter of `FindAllMarkers`.
#' @param logfc.threshold The parameter of `FindAllMarkers`.
#' @param test.use The test used in `FindAllMarkers`.
#'
#' @param selection.method Method for selecting spatially variable features,
#' `markvariogram` or `moransi`
#' @param n.top.show The number of top genes shown in the figure
#' @param n.col The number of columns in the output figure
#'
#' @param h5ad_path The path of MERFISH data in h5ad format
#' @param counts_path The path of counts data in csv format
#' @param counts_transpose Whether to transpose the counts matrix. If rows stand for genes, it will be `FALSE`.
#' @param coordinates_path The path of coordinates data in xlsx format
#' @param coordinates_index_col The column index of coordinates file
#' @param commot.signaling_type The parameter of `Commot` to determine the type of interaction.
#' Choose from 'Secreted Signaling', 'Cell-Cell Contact', and 'ECM-Receptor' for CellChatDB or
#' 'Secreted Signaling' and 'Cell-Cell Contact' for CellPhoneDB_v4.0. If None, all pairs in the database are returned.
#' @param commot.database The parameter of `Commot` to determine the database used. 'CellChat' or 'CellPhoneDB_v4.0'
#' @param commot.min_cell_pct The parameter of `Commot`. The minimum expression percentage required for LR pairs to be kept
#' @param commot.dis_thr The parameter of `Commot`. The threshold of spatial distance of signaling
#' @param commot.n_permutations The parameter of `Commot`. Number of label permutations for computing the p-value
#'
#' @param @param Step2_QC A bool value indicating whether performing QC.
#' @param Step3_Clustering A bool value indicating whether performing clustering.
#' @param Step4_Find_DEGs A bool value indicating whether finding DEGs.
#' @param Step5_SVFs A bool value indicating whether identifying spatially variable genes.
#' @param Step6_Interaction A bool value indicating whether performing spatial interaction analysis.
#'
#' @param verbose verbose as `Seurat`
#' @param pythonPath The path to the Python environment to use for the analysis.
#'
#' @details
#' The st_MERFISH_pipeline function encapsulates the complete downstream analysis pipeline
#' for MERFISH data within HemascopeR. It takes in various parameters to customize the analysis
#' workflow and provided analyzed results and publication-quality vector images for each step of the analysis workflow.
#'
#' The function performs the following steps:
#'
#' 1. Loading data.
#'
#' 2. QC (Quality Control) of genes and spots.
#'
#' 3. Normalization, PCA, and clustering.
#'
#' 4. Finding differentially expressed genes in each cluster.
#'
#' 5. Identifying spatially variable features.
#'
#' 6. Performing spatial interaction analysis using Commot.
#'
#'
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#' @import Seurat
#' @import stringr
#' @import ggplot2
#' @import patchwork
#' @import RColorBrewer
#' @import dplyr
#' @import knitr
#' @import kableExtra
#' @import SeuratDisk
#' @export
#'
st_stereo_pipeline <- function(
        input.data.dir,
        output.dir,
        sampleName = 'Hema_stereo',

        # For Step1 Loading
        data_type = 'gem',
        sep = '\t',
        bin_type = 'bins',
        bin_size = 100,
        spot_diameter = 80,
        is_sparse = TRUE,
        gene_list = NULL,
        region = NULL,
        assay = 'Spatial',

        # For Step2 QC
        Step2_QC = TRUE,
        min.gene = 20,
        min.nUMI = 50,
        max.gene = Inf,
        max.nUMI = Inf,
        min.spot = 0, # for genes
        bool.remove.mito = FALSE,
        species = 'mouse', # human or mosue

        # For Step3 Clustering
        Step3_Clustering = TRUE,
        normalization.method = 'SCTransform',
        npcs = 50,
        pcs.used = 1:10,
        resolution = 0.1,
        max.n.cluster = 30,

        # For Step4 Find DEGs
        Step4_Find_DEGs = TRUE,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25,
        test.use = 'wilcox',

        # For Step5 SVF
        Step5_SVFs = TRUE,
        selection.method = 'moransi',
        n.top.show = 10,
        n.col.show = 5,

        # For Step6 Interaction
        Step6_Interaction = TRUE,
        h5ad_path = NULL,
        counts_path = NULL,
        coordinates_path = NULL,
        coordinates_index_col = 0,
        counts_transpose = TRUE,
        commot.signaling_type = 'Secreted Signaling',
        commot.database = 'CellChat',
        commot.min_cell_pct = 0.05,
        commot.dis_thr = 500,
        commot.n_permutations = 100,

        # For Step7 Cellcycle
        Step7_Cellcycle = TRUE,
        s.features = NULL,
        g2m.features = NULL,

        verbose = FALSE,
        pythonPath = NULL
){
    ### Param ###
    SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
    FeatureColors.bi <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu')))
    FeatureColors.one <- colorRampPalette(colors = brewer.pal(n = 9, name = 'YlOrRd'))

    if(!dir.exists(file.path(output.dir, sampleName))){
        dir.create(file.path(output.dir, sampleName))
    }else{
        warning(paste0('The new results will overwrite the file under ',
                       file.path(output.dir, sampleName)))
    }
    output.dir <- file.path(output.dir, sampleName)

    print(paste0('The results will be saved in ', output.dir))

    #### Step1: Loading data ####
    print('Loading data...')
    st_obj <- stereo_Loading_Data(
        input.data.dir = input.data.dir,
        output.dir = file.path(output.dir, 'Step1_Loading_Data'),
        sampleName = sampleName,
        data_type = data_type,
        sep = sep,
        bin_type = bin_type,
        bin_size = bin_size,
        assay = assay,
        spot_diameter = spot_diameter,
        is_sparse = is_sparse,
        gene_list = gene_list,
        region = region,
        pythonPath = pythonPath
    )

    #### Step2: QC
    if(Step2_QC){
        print('Performing QC of genes and spots...')
        st_obj <- QC_Spatial(
            st_obj = st_obj,
            output.dir = file.path(output.dir, 'Step2_QC'),
            min.gene = min.gene,
            min.nUMI = min.nUMI,
            max.gene = max.gene,
            max.nUMI = max.nUMI,
            min.spot = min.spot,
            species = species,
            bool.remove.mito = bool.remove.mito,
            SpatialColors = FeatureColors.bi
        )
    }

    #### Step3: Normalization, PCA and Clustering ####
    if(Step3_Clustering){
        print('Performing normalization, PCA and clustering...')
        st_obj <- st_Clustering(
            st_obj = st_obj,
            output.dir = file.path(output.dir, 'Step3_Clustering'),
            assay = assay,
            normalization.method = normalization.method,
            npcs = npcs,
            pcs.used = pcs.used,
            resolution = resolution,
            max.n.cluster = max.n.cluster,
            verbose = verbose
        )
    }

    #### Step4: Differential expressed genes ####
    if(Step4_Find_DEGs){
        print('Finding differential expressed genes in each cluster...')
        st.markers <- st_Find_DEGs(
            st_obj = st_obj,
            output.dir = file.path(output.dir, 'Step4_Find_DEGs'),
            ident.label = 'seurat_clusters',
            only.pos = only.pos,
            min.pct = min.pct,
            logfc.threshold = logfc.threshold,
            test.use = test.use,
            verbose = verbose
        )
    }

    #### Step5: Spatially variable features ####
    if(Step5_SVFs){
        print('Identifying spatially variable features... This step may be memory-consuming.')
        st_obj <- st_SpatiallyVariableFeatures(
            st_obj = st_obj,
            output.dir = file.path(output.dir, 'Step5_SpatiallyVariableFeatures'),
            assay = st_obj@active.assay,
            selection.method = selection.method,
            n.top.show = n.top.show,
            n.col = n.col.show,
            verbose = verbose
        )
    }

    #### Step6: Spatial interaction ####
    if(Step6_Interaction){
        print('Performing spatial interaction analysis using Commot...')
        interaction_path = file.path(output.dir, 'Step6_Interaction')
        if(!dir.exists(interaction_path)){
            dir.create(interaction_path)
        }
        write.csv(st_obj@meta.data,
                  file.path(interaction_path, 'metadata.csv'),
                  row.names = TRUE)
        st_Interaction(
            st_data_path = input.data.dir,
            metadata_path = file.path(interaction_path, 'metadata.csv'),
            label_key = 'seurat_clusters',
            save_path = interaction_path,
            species = species,
            signaling_type = commot.signaling_type,
            database = commot.database,
            min_cell_pct = commot.min_cell_pct,
            dis_thr = commot.dis_thr,
            n_permutations = commot.n_permutations,
            pythonPath = pythonPath
            #condaenv = condaenv
        )
    }

    #### Step7: Cell cycle analysis ####
    if(Step7_Cellcycle){
        print('Performing cell cycle analysis...')
        st_obj <- MERFISH_Cell_cycle(
            st_obj = st_obj,
            save_path = file.path(output.dir, 'Step7_Cellcycle'),
            s.features = s.features,
            g2m.features = g2m.features,
            species = species,
            FeatureColors.bi = FeatureColors.bi
        )
    }

    #### Save data ####
    output.dir.final <- file.path(output.dir, 'Data_and_report')
    if(!dir.exists(output.dir.final)){
        dir.create(output.dir.final)
    }
    saveRDS(st_obj, file.path(output.dir.final, 'st_object.rds'))
}

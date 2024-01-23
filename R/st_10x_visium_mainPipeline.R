#' The pipeline for analyzing 10x Visium data.
#'
#' Function \code{st_10x_visium_pipeline} encompasses the complete downstream
#' analysis pipeline for 10x Visium data within HemascopeR.
#' @param input.data.dirs A character of data path, where there are filtered_feature_bc_matrix.h5
#' and spatial folder.
#' @param output.dir A character of path to store the results and figures.
#' @param sampleName A character of the name of the sample.
#' @param min.gene An integer representing the minimum number of genes detected in a spot.
#' @param max.gene An integer representing the maximum number of genes detected in a spot.
#' @param min.nUMI An integer representing the minimum number of nUMI detected in a spot.
#' @param max.nUMI An integer representing the maximum number of nUMI detected in a spot.
#' @param min.spot An integer representing the minimum number of spots expressing each gene.
#' @param normalization.method 'SCTransform' or other choice of `normalization.method` in `NormalizaData`.
#' @param resolution An float parameter of `FindClusters` of `Seurat`.
#' @param npcs The total number of PCs to compute.
#' @param pcs.used Dimensions of PCs to use as input.
#' @param only.pos A bool value to indicate whether only return positive markers.
#' @param min.pct The parameter of `FindAllMarkers`.
#' @param logfc.threshold The parameter of `FindAllMarkers`.
#' @param test.use The test used in `FindAllMarkers`.
#' @param selection.method Method for selecting spatially variable features, `markvariogram` or `moransi`.
#'
#' @param commot.signaling_type The parameter of `Commot` to determine the type of interaction.
#' Choose from 'Secreted Signaling', 'Cell-Cell Contact', and 'ECM-Receptor' for CellChatDB or.
#' 'Secreted Signaling' and 'Cell-Cell Contact' for CellPhoneDB_v4.0. If None, all pairs in the database are returned.
#' @param commot.database The parameter of `Commot` to determine the database used. 'CellChat' or 'CellPhoneDB_v4.0'.
#' @param commot.min_cell_pct The parameter of `Commot`. The minimum expression percentage required for LR pairs to be kept.
#' @param commot.dis_thr The parameter of `Commot`. The threshold of spatial distance of signaling.
#' @param commot.n_permutations The parameter of `Commot`. Number of label permutations for computing the p-value.
#'
#' @param copykat.genome 'hg20' or 'mm10'
#' @param copykat.LOW.DR The parameter of `copykat`. The minimal population fractions of genes for smoothing.
#' @param copykat.UP.DR The parameter of `copykat`. The minimal population fractions of genes for segmentation.
#' @param copykat.win.size The parameter of `copykat`. The minimal window sizes for segmentation.
#' @param copykat.distance The parameter of `copykat`. Distance methods include euclidean,
#' and correlation converted distance include pearson and spearman.
#' @param copykat.n.cores The parameter of `copykat`. The number of cores for parallel computing.
#'
#' @param cell2loc.sc.h5ad.dir A character of path of h5ad file of scRNA-seq reference data, default NULL and used default dataset.
#' @param cell2loc.sc.max.epoch A integer representing the maximum epochs of training scRNA-seq data.
#' @param cell2loc.st.max.epoch A integer representing the maximum epochs of training ST data.
#' @param cell2loc.use.gpu A bool value indicating whether to use GPU.
#' @param cell2loc.use.Dataset 'HematoMap' or 'LymphNode'.
#'
#' @param coexistence.method The method to analyze the coexistence of cell types. 'correlation' or 'Wasserstein'.
#' @param Niche.cluster.n The number of clusters to cluster.
#'
#' @param verbose verbose as `Seurat`.
#' @param species A character representing the species of sample, `human` or `mouse`.
#' @param bool.remove.mito A bool value indicating whether removing mitochondrial gene.
#' @param bool.regress.cycling A bool value indicating whether regressing cycling.
#' @param bool.regress.cycling.standard A bool value indicating the method using in
#' regressing cycling. See more at `https://satijalab.org/seurat/articles/cell_cycle_vignette.html`.
#' @param s.features A gene vector used to calculate the S phase score. If NULL,
#' the default gene set will be used.
#' @param g2m.features A gene vector used to calculate the G2/M phase score. If NULL,
#' the default gene set will be used.
#' @param genReport A bool value indicating whether generating the report.
#' @param Step2_QC A bool value indicating whether performing QC.
#' @param Step3_Clustering A bool value indicating whether performing clustering.
#' @param Step4_Find_DEGs A bool value indicating whether finding DEGs.
#' @param Step5_SVFs A bool value indicating whether identifying spatially variable genes.
#' @param Step6_Interaction A bool value indicating whether performing spatial interaction analysis.
#' @param Step7_CNV A bool value indicating whether performing CNV analysis.
#' @param Step8_Deconvolution A bool value indicating whether performing deconvolution analysis.
#' @param Step9_Cellcycle A bool value indicating whether analyzing cell cycle.
#' @param Step10_Niche A bool value indicating whether performing niche analysis.
#' @param pythonPath The path to the Python environment to use for the analysis.
#' @details
#' The st_10x_visium_pipeline function encapsulates the complete downstream analysis pipeline
#' for 10x Visium data within HemascopeR. It takes in various parameters to customize the analysis
#' workflow and generates a well-formatted HTML analysis report along with meticulously organized
#' and provided analyzed results and publication-quality vector images for each step of the analysis workflow.
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
#' 7. CNV (Copy Number Variation) analysis using copykat.
#'
#' 8. Deconvolution using cell2location.
#'
#' 9. Cell cycle analysis.
#'
#' 10. Niche analysis.
#'
#' Finally, it saves the analyzed data and generates the report.
#'
#' @return Return the user with a well-formatted HTML analysis report.
#' Additionally, meticulously organize and provide the analyzed results and
#' publication-quality vector images for each step of the analysis workflow.
#' @author Zhenyi Wang wangzy17@tsinghua.org.cn and Yuxin Miao miaoyx21@mails.tsinghua.edu.cn
#' @import Seurat
#' @import stringr
#' @import ggplot2
#' @import patchwork
#' @import RColorBrewer
#' @import SeuratDisk
#' @import dplyr
#' @export
#'

st_10x_visium_pipeline <- function(
        input.data.dir,
        output.dir,
        sampleName = 'Hema_ST',

        # For Step2 QC
        Step2_QC = TRUE,
        min.gene = 200,
        min.nUMI = 500,
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
        resolution = 0.8,

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
        commot.signaling_type = 'Secreted Signaling',
        commot.database = 'CellChat',
        commot.min_cell_pct = 0.05,
        commot.dis_thr = 500,
        commot.n_permutations = 100,

        # For Step7 CNV analysis
        Step7_CNV = TRUE,
        copykat.genome = NULL,
        copykat.LOW.DR = 0.05,
        copykat.UP.DR = 0.1,
        copykat.win.size = 25,
        copykat.distance = "euclidean",
        copykat.n.cores = 1,

        # For Step8 Deconvolution
        Step8_Deconvolution = TRUE,
        cell2loc.sc.h5ad.dir = NULL,
        cell2loc.sc.max.epoch = 1000,
        cell2loc.st.max.epoch = 10000,
        cell2loc.use.gpu = FALSE,
        cell2loc.use.Dataset = 'LymphNode',

        # For Step9 Cellcycle
        Step9_Cellcycle = TRUE,
        s.features = NULL,
        g2m.features = NULL,

        # For Step10 Nich
        Step10_Niche = TRUE,
        coexistence.method = 'correlation',
        Niche.cluster.n = 4,

        # settings
        #condaenv = 'r-reticulate',
        verbose = FALSE,
        genReport = TRUE,
        pythonPath = NULL
){
    ### Param ###
    SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
    FeatureColors.bi <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = 'RdYlBu')))
    FeatureColors.one <- colorRampPalette(colors = brewer.pal(n = 9, name = 'YlOrRd'))
    if(is.null(pythonPath)==FALSE){ reticulate::use_python(pythonPath) }else{print('Please set the path of Python.')}

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
    st_obj <- Load10X_Spatial(input.data.dir,
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "slice1",
                              filter.matrix = TRUE,
                              to.upper = FALSE,
                              image = NULL)
    st_obj@project.name <- sampleName
    st_obj$orig.ident <- factor(sampleName)
    st_obj@active.ident <- st_obj$orig.ident

    #### Step2: QC ####
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
            SpatialColors = SpatialColors
        )
    }


    #### Step3: Normalization, PCA and Clustering ####
    if(Step3_Clustering){
        print('Performing normalization, PCA and clustering...')
        st_obj <- st_Clustering(
            st_obj = st_obj,
            output.dir = file.path(output.dir, 'Step3_Clustering'),
            normalization.method = normalization.method,
            npcs = npcs,
            pcs.used = pcs.used,
            resolution = resolution,
            verbose = verbose
        )
    }

    #### Step4: Differentially expressed genes ####
    if(Step4_Find_DEGs){
        print('Finding differentially expressed genes in each cluster...')
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
        print('Identifying spatially variable features...')
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
        st_interaction(
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

    #### Step7: CNV analysis ####
    if(Step7_CNV){
        print('Performing CNV analysis using copykat...')
        st_obj <- st_CNV(
            st_obj = st_obj,
            save_path = file.path(output.dir, 'Step7_CNV_analysis'),
            assay = 'Spatial',
            LOW.DR = copykat.LOW.DR,
            UP.DR = copykat.UP.DR,
            win.size = copykat.win.size,
            distance = copykat.distance,
            genome = copykat.genome,
            n.cores = copykat.n.cores,
            species = species
        )
    }

    #### Step8: Deconvolution ####
    if(Step8_Deconvolution){
        print('Performing deconvolution using cell2location...')
        st_obj <- st_Deconvolution(
            st.data.dir = input.data.dir,
            sc.h5ad.dir = cell2loc.sc.h5ad.dir,
            st_obj = st_obj,
            save_path = file.path(output.dir, 'Step8_Deconvolution'),
            sc.labels.key = 'seurat_clusters',
            species = species,
            sc.max.epoch = cell2loc.sc.max.epoch,
            st.max.epoch = cell2loc.st.max.epoch,
            use.gpu = cell2loc.use.gpu,
            use.Dataset = cell2loc.use.Dataset,
            pythonPath = pythonPath
            # condaenv = condaenv
        )
    }

    #### Step9: Cell cycle analysis ####
    if(Step9_Cellcycle){
        print('Performing cell cycle analysis...')
        st_obj <- st_Cell_cycle(
            st_obj = st_obj,
            save_path = file.path(output.dir, 'Step9_Cellcycle'),
            s.features = s.features,
            g2m.features = g2m.features,
            species = species,
            FeatureColors.bi = FeatureColors.bi
        )
    }

    #### Step10: Niche analysis ####
    if(Step10_Niche){
      print('Performing niche analysis...')
      if(!file.exists(file.path(output.dir, 'Step8_Deconvolution', 'cell2loc_res.csv'))){
        stop('Please run step 8 deconvolution first.')
      }
      tmp <- read.csv(file.path(output.dir, 'Step8_Deconvolution', 'cell2loc_res.csv'),
                      row.names = 1)
      features <- colnames(tmp)

      if(!all(features %in% names(st_obj@meta.data))){
        common.barcodes <- intersect(colnames(st_obj), rownames(tmp))
        tmp <- tmp[common.barcodes, ]
        st_obj <- st_obj[, common.barcodes]
        st_obj <- AddMetaData(st_obj,
                              metadata = tmp)
      }

      rm(tmp)
      gc()

      st_obj <- st_NicheAnalysis(
        st_obj,
        features = features,
        save_path = file.path(output.dir, 'Step10_NicheAnalysis'),
        coexistence.method = coexistence.method,
        kmeans.n = Niche.cluster.n,
        st_data_path = file.path(output.dir, 'Step1_Loading_Data'),
        slice = slice,
        species = species,
        # condaenv = condaenv,
        pythonPath = pythonPath
      )
    }

    #### Save data ####
    output.dir.final <- file.path(output.dir, 'Data_and_report')
    if(!dir.exists(output.dir.final)){
        dir.create(output.dir.final)
    }
    saveRDS(st_obj, file.path(output.dir.final, 'st_object.rds'))
    suppressMessages(suppressWarnings(
        SaveH5Seurat(st_obj,
                     filename = file.path(output.dir.final, 'st_object.h5Seurat'),
                     overwrite = TRUE)
    ))
    suppressMessages(suppressWarnings(
        Convert(file.path(output.dir.final, 'st_object.h5Seurat'), dest = "h5ad",
                overwrite = TRUE)
    ))

    #### Generate the report ####
    print('Generating the report...')
    if(genReport){
        # knitr::knit(file.path(system.file(package = "HemaScopeR"), "rmd/st_base.Rmd"),
        #             file.path(output.dir.final, 'st_pipeline.md'))
        knitr::knit(file.path("/public/home/rp1008csj/rp1020wangzy/HemaScopeR/rmd/st_base.Rmd"),
                    file.path(output.dir.final, 'st_pipeline.md'))
        markdown::markdownToHTML(file.path(output.dir.final, 'st_pipeline.md'),
                                 file.path(output.dir.final, 'st_pipeline.html'))
    }
}




# We deposit assiatant functions in this file.
# And would not be exported to users.

# visualization theme-----------------------------------------------------------
mytheme <- theme(panel.grid.major =element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(size = 1,
                                          colour = "black"),
                 axis.title.x =element_text(size=20,
                                            family = "sans",
                                            color = "black",
                                            face = "bold"),
                 axis.text.x = element_text(size = 20,
                                            family = "sans",
                                            color = "black",
                                            face = "bold",
                                            vjust = 0,
                                            hjust = 0),
                 axis.text.y = element_text(size = 20,
                                            family = "sans",
                                            color = "black",
                                            face = "bold",
                                            vjust = 0,
                                            hjust = 1),
                 axis.title.y=element_text(size=20,
                                           family = "sans",
                                           color = "black",
                                           face = "bold"),
                 legend.text = element_text(size=15,
                                            family = "sans",
                                            color = "black",
                                            face = "bold"),
                 legend.title = element_text(size=15,
                                             family = "sans",
                                             color = "black",
                                             face = "bold"),
                 legend.background = element_blank(),
                 legend.key = element_blank()
)

# GO enrichment------------------------------------------------------------------
#' HemaScopeREnrichment - Function for Gene Ontology (GO) Enrichment Analysis
#' 
#' This function performs Gene Ontology (GO) enrichment analysis for a set of differentially expressed genes (DEGs) grouped by clusters. It generates GO enrichment results in terms of Biological Process (BP) category, including barplots and CSV files.
#'
#' @param DEGs A data frame containing differential gene expression information with columns "cluster" and "gene". The "cluster" column specifies the cluster to which each gene belongs, and the "gene" column contains the gene names.
#' @param OrgDb A database for gene annotation, used to map gene symbols to GO terms. e.g. It should be 'org.Mm.eg.db' for mmu and 'org.Hs.eg.db' for hsa.
#' @param output.dir The directory path where the output files (PDF and CSV) will be saved.
#'
#' @details The function takes a data frame of DEGs and performs GO enrichment analysis separately for each cluster defined in the "cluster" column. It utilizes the "enrichGO" function to calculate GO enrichment, and generates barplots of significant enriched GO terms for each cluster. The results are saved as both PDF files and CSV files, providing insights into the biological processes associated with the DEGs in different clusters.
#'
#' @return The function does not return any values directly, but it saves PDF files with barplots showing significant GO enrichment results for each cluster, and CSV files containing the enriched GO terms and associated statistics.
#'
#' @author Zhenyi Wang <wangzy17@tsinghua.org.cn> and Yuxin Miao <miaoyx21@mails.tsinghua.edu.cn>
#'
#' @export
HemaScopeREnrichment = function(DEGs = NULL,
                                OrgDb = NULL,
                                output.dir = NULL){
    DEGs$cluster <- make.names(DEGs$cluster)
    uniqueClusters <- unique(DEGs$cluster)
    for(i in uniqueClusters){
        DEGs_temp <- subset(DEGs, DEGs$cluster == i)
        DEGs_temp <- DEGs_temp$gene
        
        #enrich
        GO_BP <- enrichGO(gene = DEGs_temp,
                          keyType = "SYMBOL",
                          OrgDb = OrgDb, 
                          ont = "BP",
                          pAdjustMethod = "fdr",
                          pvalueCutoff = 0.2,
                          qvalueCutoff  = 0.2,
                          minGSSize = 3, 
                          maxGSSize = 500,
                          readable = FALSE)

        GO_BP.result <- as.data.frame(GO_BP@result)
        GO_BP.result <- GO_BP.result[which(GO_BP.result$pvalue <= 0.2),]
        GO_BP.result <- GO_BP.result[which(GO_BP.result$qvalue <= 0.2),]

        if(nrow(GO_BP.result) > 0){
            #plot
            pdf(file=paste(output.dir, as.character(i), '_BP_GOenrich.pdf',sep=''), width=6, height=15)
            print(barplot(GO_BP, 
                          x="Count", 
                          color="qvalue",
                          showCategory=30,
                          font.size=8,
                          title=paste(as.character(i),'_BP_GOenrich', sep=''))+scale_y_discrete(labels=function(x) str_wrap(x, width=60)))
            dev.off()
                                                                                                 
            #data frame
            GO_BP.result <- GO_BP.result[order(GO_BP.result[,9], decreasing = TRUE),]
            write.csv(GO_BP.result, file=paste(output.dir, as.character(i),'_BP_GOenrich.csv',sep=''))
        }                                                                                     
    }
}

# XGR enrichment-------------------------------------------------------------------------------------------------------------------
OpenXGR_SAG = function(sc_object.markers = NULL,
                       output.dir = NULL,
                       subnet.size = 10){
    for(i in unique(sc_object.markers$cluster)){
        SAGdata <- sc_object.markers[which(sc_object.markers$cluster == i), ]
        SAGdata <- SAGdata[,c('gene','p_val')]
        #parameter 
        placeholder <- "http://www.comptransmed.com/bigdata_openxgr"
        network <- "STRING_high"
        
        #subnetwork analysis 
        SAGig <- oDefineNet(network=network, STRING.only=c("experimental_score","database_score"), placeholder=placeholder)
        SAGig2 <- oNetInduce(SAGig, nodes_query=V(SAGig)$name, largest.comp=T) %>% as.undirected()
        SAGsubg <- oSubneterGenes(SAGdata, network.customised=SAGig2, subnet.size=subnet.size, placeholder=placeholder)
        if(!is.null(SAGsubg)){
            #crosstalk table 
            SAGdf_subg <- SAGsubg %>% oIG2TB("nodes") %>% 
                          transmute(Genes=name, Pvalue=as.numeric(significance),
                          Description=description) %>% arrange(Pvalue)
            write.csv(SAGdf_subg, paste0(output.dir,'/SAG-crosstalk','_cluster_', i ,'.csv'))
            if(nrow(SAGdf_subg) > 2){
                #network 
                SAGsubg <- SAGsubg %>% oLayout(c("layout_with_kk","graphlayouts.layout_with_stress")[2])
                SAGvec <- V(SAGsubg)$significance %>% as.numeric()
                SAGvec[SAGvec==0] <- min(SAGvec[SAGvec!=0])
                V(SAGsubg)$logP <- -log10(SAGvec)
                SAGvec <- -log10(SAGvec)
                if(max(SAGvec)<20){
                  zlim <- c(0, ceiling(max(SAGvec)))
                }else{
                  zlim <- c(0, floor(max(SAGvec)/10)*10)
                }
                SAGnetwork <- oGGnetwork(g=SAGsubg, node.label="name", node.label.size=3, node.label.color="black",
                                         node.label.alpha=0.95, node.label.padding=0.5, node.label.arrow=0, 
                                         node.label.force=0.4, node.shape=19, node.xcoord="xcoord", node.ycoord="ycoord",
                                         node.color="logP", node.color.title=expression(-log[10]("pvalue")), 
                                         colormap="spectral.top", zlim=zlim, node.size.range=5, title="", edge.color="steelblue4",
                                         edge.color.alpha=0.5, edge.size=0.3, edge.curve=0.05)
                ggsave(paste0(output.dir, '/SAG-network','_cluster_', i ,'.pdf'),SAGnetwork,width = 6, height = 6)
            }
          }
    }
}
                                                                                                
# # cbind matrix
# cbind_all = function(matrix_list) {
#   result <- matrix_list[[1]]  
#   for (i in 2:length(matrix_list)) {
#     result <- cbind(result, matrix_list[[i]])  
#   }
#   return(result)
# }

# # rbind matrix
# rbind_all = function(matrix_list) {
#   result <- matrix_list[[1]]  
#   for (i in 2:length(matrix_list)) {
#     result <- rbind(result, matrix_list[[i]])  
#   }
#   return(result)
# }

# Rename the genes in Seurat object
# This function is referenced from the forum https://www.jianshu.com/p/6495706bac53.
RenameGenesSeurat <- function(obj = NULL,
                              newnames = NULL,
                              gene.use = NULL, 
                              de.assay = "RNA",
                              lassays = "RNA") {
    # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration.
    # It only changes obj@assays$RNA@counts, @data and @scale.data.
    print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features,@reductions$pca@feature.loadings")
    #names(obj@assays)
    assay.use <- obj@reductions$pca@assay.used
    DefaultAssay(obj) <- de.assay
    if (is.null(gene.use)) {
     all_genenames <- rownames(obj)
    }else{
     all_genenames <- gene.use
     obj <- subset(obj,features=gene.use)
    }
    
    order_name <- function(v1,v2,ref){
     #v2 <- make.names(v2,unique=T)
     df1 <- data.frame(v1,v2)
     rownames(df1) <- df1$v1
     df1 <- df1[ref,]
     return(df1)
    }
    
    df1 <- order_name(v1=all_genenames,
                      v2=newnames,
                      ref=rownames(obj))
    all_genenames <- df1$v1
    newnames <- df1$v2
    
    if ('SCT' %in% lassays) {
      if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
        obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
        rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
      }
    }
    
    change_assay <- function(a1=de.assay,
                             obj,
                             newnames=NULL,
                             all_genenames=NULL){
    RNA <- obj@assays[a1][[1]]
      if (nrow(RNA) == length(newnames)) {
        if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
        if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
        if (length(RNA@var.features)) {
            df1 <- order_name(v1=all_genenames,v2=newnames,ref=RNA@var.features)
            all_genenames1 <- df1$v1
            newnames1 <- df1$v2
            RNA@var.features <- newnames1
        }
        if (length(RNA@scale.data)){
            df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(RNA@scale.data))
            all_genenames1 <- df1$v1
            newnames1 <- df1$v2
            rownames(RNA@scale.data) <- newnames1
        }
    
      } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
      obj@assays[a1][[1]] <- RNA
      return(obj)
    }
    
    for (a in lassays) {
      DefaultAssay(obj) <- a
      df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
      all_genenames1 <- df1$v1
      newnames1 <- df1$v2
      obj <- change_assay(obj=obj,a1=a,newnames=newnames1,all_genenames=all_genenames1)
    }
    hvg <- VariableFeatures(obj,assay=assay.use)
    if (length(obj@reductions$pca)){
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=hvg)
        df1 <- df1[rownames(obj@reductions$pca@feature.loadings),]
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        rownames(obj@reductions$pca@feature.loadings) <- newnames1
    }
    try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]])))
    return(obj)
}

# load the previous results----------------------------------------------------------------------------
Load_previous_results = function(previous_results_path=NULL){
       # Get a list of all .RDS files in the specified path
       rds_files_list <- list.files(previous_results_path, pattern = "\\.rds$", full.names = TRUE)
       # Loop through the list of RDS files and load them into variables
       for (rds_file in rds_files_list) {
         # Extract the variable name from the file name
         var_name <- tools::file_path_sans_ext(basename(rds_file))
         if(!(var_name %in% c(   # input and output
                                 'input.data.dirs',
                                 'project.names', 
                                 'output.dir',
                                 'pythonPath',
                                 # quality control and preprocessing
                                 'gene.column',
                                 'min.cells',
                                 'min.feature',
                                 'mt.pattern',
                                 'nFeature_RNA.limit',
                                 'percent.mt.limit',
                                 'scale.factor',
                                 'nfeatures',
                                 'ndims',
                                 'vars.to.regress',
                                 'PCs',
                                 'resolution',
                                 'n.neighbors',
                                 # remove doublets
                                 'doublet.percentage',
                                 'doublerFinderwraper.PCs',
                                 'doublerFinderwraper.pN',
                                 'doublerFinderwraper.pK',
                                 # phateR
                                 'phate.knn',
                                 'phate.npca',
                                 'phate.t',
                                 'phate.ndim',
                                 'min.pct',
                                 'logfc.threshold',
                                 # visualization
                                 'marker.genes',
                                 'ViolinPlot.cellTypeOrders',
                                 'ViolinPlot.cellTypeColors',
                                 'Org',
                                 'lineage.genelist',
                                 'lineage.names',
                                 'groups_colors',
                                 #slingshot
                                 'slingshot.start.clus',
                                 'slingshot.end.clus',
                                 'slingshot.colors',
                                 'loom.files.path',
                                 # cell cycle
                                 'cellcycleCutoff',
                                 # cell chat
                                 'sorting',
                                 'ncores',
                                 # Verbose = FALSE,
                                 # activeEachStep
                                 'Whether_load_previous_results',
                                 'Step1_Input_Data',
                                 'Step1_Input_Data.type',
                                 'Step2_Quality_Control',
                                 'Step2_Quality_Control.RemoveBatches',
                                 'Step2_Quality_Control.RemoveDoublets',
                                 'Step3_Clustering',
                                 'Step4_Identify_Cell_Types',
                                 'Step4_Use_Which_Labels',
                                 'Step4_Cluster_Labels',
                                 'Step4_Changed_Labels',
                                 'Step4_run_sc_CNV',
                                 'Step5_Visualization',
                                 'Step6_Find_DEGs',
                                 'Step7_Assign_Cell_Cycle',
                                 'Step8_Calculate_Heterogeneity',
                                 'Step9_Violin_Plot_for_Marker_Genes',
                                 'Step10_Calculate_Lineage_Scores',
                                 'Step11_GSVA',
                                 'Step11_GSVA.identify.cellType.features',
                                 'Step11_GSVA.identify.diff.features',
                                 'Step11_GSVA.comparison.design',
                                 'Step12_Construct_Trajectories',
                                 'Step12_Construct_Trajectories.monocle',
                                 'Step12_Construct_Trajectories.slingshot',
                                 'Step12_Construct_Trajectories.scVelo',
                                 'Step13_TF_Analysis',
                                 'Step13_TF_Analysis.cellTypes_colors',
                                 'Step13_TF_Analysis.groups_colors',
                                 'Step14_Cell_Cell_Interaction',
                                 'Step15_Generate_the_Report'
                                 ))){  
           # Load the RDS file and assign it to the variable with the same name 
           rds_file.temp <- base::readRDS(rds_file)
           assign(var_name, rds_file.temp, envir = .GlobalEnv)
          }
       }
    }
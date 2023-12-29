#' Function to define a gene network
#'
#' \code{oDefineNet} is supposed to define a gene network sourced from the STRING database or the Pathway Commons database. It returns an object of class "igraph". 
#'
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGGc_Metabolism' for pathways grouped into 'Metabolism', 'KEGGc_GeneticProcess' for 'Genetic Process' pathways, 'KEGGc_EnvironmentalProcess' for 'Environmental Process' pathways, 'KEGGc_CellularProcess' for 'Cellular Process' pathways, 'KEGGc_OrganismalSystem' for 'Organismal System' pathways, and 'KEGGc_HumanDisease' for 'Human Disease' pathways. For KEGG 'Environmental Process' pathways, it can be further into 'KEGGs_ImmuneSystem','KEGGs_EndocrineSystem','KEGGs_CirculatorySystem','KEGGs_DigestiveSystem','KEGGs_ExcretorySystem','KEGGs_NervousSystem','KEGGs_SensorySystem','KEGGs_DevelopmentAndRegeneration','KEGGs_Aging','KEGGs_EnvironmentalAdaptation'. 'REACTOME' for protein-protein interactions derived from Reactome pathways. 'HuRI' for human reference direct interactome (biophysical), while 'HuRI_union' for all HuRI identified so far, 'HuRI_litbm' for literature-derived direct interactions, and 'HuRI_all' for the union of 'HuRI_union' and 'HuRI_litbm'
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param weighted logical to indicate whether edge weights should be considered. By default, it sets to false. If true, it only works for the network from the STRING database 
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param placeholder the characters to tell the placeholder of RDS files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return
#' an object of class "igraph"
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{oRDS}}, \code{\link{oIG2TB}}, \code{\link{oTB2IG}}
#' @include oDefineNet.r
#' @examples
#' \dontrun{
#' # STRING (high quality)
#' g <- oDefineNet("STRING_high", placeholder=placeholder)
#' # STRING (high quality), with edges weighted 
#' g <- oDefineNet("STRING_high", weighted=TRUE, placeholder=placeholder)
#' # STRING (high quality), only edges sourced from experimental or curated data
#' g <- oDefineNet("STRING_high", STRING.only=c("experimental_score","database_score"), placeholder=placeholder)
#' 
#' # Pathway Commons 
#' g <- oDefineNet("PCommonsDN_medium", placeholder=placeholder)
#' 
#' # KEGG (all)
#' g <- oDefineNet("KEGG", placeholder=placeholder)
#' # KEGG category ('Organismal Systems')
#' g <- oDefineNet("KEGGc_OrganismalSystem", placeholder=placeholder)
#' # KEGG subcategory ('Immune Systems')
#' g <- oDefineNet("KEGGs_ImmuneSystem", placeholder=placeholder)
#' 
#' # REACTOME 
#' g <- oDefineNet("REACTOME", placeholder=placeholder)
#' }

oDefineNet <- function(network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG", "KEGGc_CellularProcess","KEGGc_EnvironmentalProcess","KEGGc_GeneticProcess","KEGGc_HumanDisease","KEGGc_Metabolism","KEGGc_OrganismalSystem", "KEGGs_Aging","KEGGs_AminoAcidMetabolism","KEGGs_BiosynthesisOfOtherSecondaryMetabolites","KEGGs_CancerOverview","KEGGs_CancerSpecificTypes","KEGGs_CarbohydrateMetabolism","KEGGs_CardiovascularDisease","KEGGs_CellGrowthAndDeath","KEGGs_CellMotility","KEGGs_CellularCommunityEukaryotes","KEGGs_CirculatorySystem","KEGGs_DevelopmentAndRegeneration","KEGGs_DigestiveSystem","KEGGs_DrugResistanceAntineoplastic","KEGGs_EndocrineAndMetabolicDisease","KEGGs_EndocrineSystem","KEGGs_EnergyMetabolism","KEGGs_EnvironmentalAdaptation","KEGGs_ExcretorySystem","KEGGs_FoldingSortingAndDegradation","KEGGs_GlycanBiosynthesisAndMetabolism","KEGGs_ImmuneDisease","KEGGs_ImmuneSystem","KEGGs_InfectiousDiseaseBacterial","KEGGs_InfectiousDiseaseParasitic","KEGGs_InfectiousDiseaseViral","KEGGs_LipidMetabolism","KEGGs_MetabolismOfCofactorsAndVitamins","KEGGs_MetabolismOfOtherAminoAcids","KEGGs_MetabolismOfTerpenoidsAndPolyketides","KEGGs_NervousSystem","KEGGs_NeurodegenerativeDisease","KEGGs_NucleotideMetabolism","KEGGs_ReplicationAndRepair","KEGGs_SensorySystem","KEGGs_SignalTransduction","KEGGs_SignalingMoleculesAndInteraction","KEGGs_SubstanceDependence","KEGGs_Translation","KEGGs_TransportAndCatabolism","KEGGs_XenobioticsBiodegradationAndMetabolism", "REACTOME", "HuRI", "HuRI_all","HuRI_union","HuRI_litbm"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], weighted=FALSE, verbose=TRUE, placeholder=NULL, guid=NULL)
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    
	if(verbose){
		message(sprintf("Load the network %s (%s) ...", network, as.character(Sys.time())), appendLF=TRUE)
	}
	
	if(network %>% stringr::str_detect('STRING')){

		## restrict to those edges with given confidence
		flag <- network %>% stringr::str_split("_", simplify=TRUE)
		combined_score <- NULL
		if(flag[2]=='highest'){
			g <- oRDS('org.Hs.string_highest', placeholder=placeholder, guid=guid, verbose=verbose)
		}else if(flag[2]=='high'){
			g <- oRDS('org.Hs.string_high', placeholder=placeholder, guid=guid, verbose=verbose)
		}else if(flag[2]=='medium'){
			g <- oRDS('org.Hs.string_medium', placeholder=placeholder, guid=guid, verbose=verbose)
		}else if(flag[2]=='low'){
			g <- oRDS('org.Hs.string_low', placeholder=placeholder, guid=guid, verbose=verbose)
		}
		
		## nodes
		g %>% oIG2TB("nodes") %>% dplyr::select_at(dplyr::vars("symbol","description")) %>% dplyr::distinct() -> nodes
		
		## extract edges (by symbol)
		### because of the way storing the network from the STRING database
		V(g)$name <- V(g)$symbol
		g %>% oIG2TB("edges") -> df_edges
		
		## further restricted by the evidence type
		default.STRING.only <- c("neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")
		ind <- match(default.STRING.only, STRING.only)
		STRING.only <- default.STRING.only[!is.na(ind)]
		if(length(STRING.only)>0){
			. <- NULL
			#df_edges %>% dplyr::filter_at(dplyr::vars(STRING.only),dplyr::any_vars(.>0)) -> df_edges
			df_edges <- df_edges %>% dplyr::filter_at(dplyr::all_of(STRING.only),dplyr::any_vars(.>0)) %>% suppressWarnings()
		}
		
		## weighted or not
		if(weighted){
			## keep the maximim weight
			from <- to <- combined_score <- NULL
			df_edges %>% dplyr::select(from,to,combined_score) %>% dplyr::group_by(from,to) %>% dplyr::summarise(weight=base::max(combined_score)) %>% dplyr::ungroup() -> edges
		}else{
			df_edges %>% dplyr::select_at(dplyr::vars("from","to")) %>% dplyr::distinct() %>% dplyr::mutate(weight=1) -> edges
		}
		
		## "intersected=TRUE": make sure vertices in nodes and edges are the same
		g <- oTB2IG(edges=edges, nodes=nodes, directed=FALSE, intersected=TRUE)
			
    }else if(network %>% stringr::str_detect('PCommonsUN')){
    
		g <- oRDS('org.Hs.PCommons_UN', placeholder=placeholder, guid=guid, verbose=verbose)
		
		## extract relations
		g %>% oIG2TB("edges") -> df_edges
		
		## restrict to those edges with given confidence
		flag <- network %>% stringr::str_split("_", simplify=TRUE)
		if(flag[2]=='high'){
			# restrict to those edges with physical interactions and with score>=102
			df_edges %>% dplyr::filter_at(dplyr::vars("in_complex_with","interacts_with"),dplyr::any_vars(.>102)) -> df_edges
		}else if(flag[2]=='medium'){
			# restrict to those edges with physical interactions and with score>=101
			df_edges %>% dplyr::filter_at(dplyr::vars("in_complex_with","interacts_with"),dplyr::any_vars(.>101)) -> df_edges
		}
		df_edges %>% dplyr::select_at(dplyr::vars("from","to")) %>% dplyr::mutate(weight=1) -> edges
		
		## extract nodes
		g %>% oIG2TB("nodes") -> df_nodes
		df_nodes %>% dplyr::select_at(dplyr::vars("symbol","description")) -> nodes
		
		## "intersected=TRUE": make sure vertices in nodes and edges are the same
		g <- oTB2IG(edges=edges, nodes=nodes, directed=FALSE, intersected=TRUE)
		
    }else if(length(grep('PCommonsDN',network,perl=TRUE)) > 0){
    	flag <- network %>% stringr::str_split("_", simplify=TRUE)
		if(flag[2]=='high'){
			g <- oRDS('org.Hs.PCommons_DN', placeholder=placeholder, guid=guid, verbose=verbose)
			## extract relations
			g %>% oIG2TB("edges") -> df_edges
			## restrict to those edges with high confidence score>=102
			df_edges %>% dplyr::filter_at(dplyr::vars("catalysis_precedes","controls_expression_of","controls_phosphorylation_of","controls_state_change_of","controls_transport_of"),dplyr::any_vars(.>102)) -> df_edges
		}else if(flag[2]=='medium'){
			g <- oRDS('org.Hs.PCommons_DN', placeholder=placeholder, guid=guid, verbose=verbose)
			## extract relations
			g %>% oIG2TB("edges") -> df_edges
			## restrict to those edges with median confidence score>=101
			df_edges %>% dplyr::filter_at(dplyr::vars("catalysis_precedes","controls_expression_of","controls_phosphorylation_of","controls_state_change_of","controls_transport_of"),dplyr::any_vars(.>101)) -> df_edges
		}else{
			ls_ig <- oRDS('org.Hs.PCommons_DN_source', placeholder=placeholder, guid=guid, verbose=verbose) %>% tibble::deframe()
			## specific source
			g <- ls_ig[[flag[2]]]
			g %>% oIG2TB("edges") -> df_edges
			## restrict to those edges with high confidence score>=101
			df_edges %>% dplyr::filter_at(dplyr::vars("catalysis_precedes","controls_expression_of","controls_phosphorylation_of","controls_state_change_of","controls_transport_of"),dplyr::any_vars(.>101)) -> df_edges
		}
		df_edges %>% dplyr::select_at(dplyr::vars("from","to")) %>% dplyr::mutate(weight=1) -> edges
		
		## extract nodes
		g %>% oIG2TB("nodes") -> df_nodes
		df_nodes %>% dplyr::select_at(dplyr::vars("symbol","description")) -> nodes
		
		## "intersected=TRUE": make sure vertices in nodes and edges are the same
		g <- oTB2IG(edges=edges, nodes=nodes, directed=TRUE, intersected=TRUE)
    
    }else if(network=='KEGG'){
    	g <- oRDS('ig.KEGG.merged', placeholder=placeholder, guid=guid, verbose=verbose)
    	g <- igraph::delete_vertex_attr(g, "hsa")
    	g <- igraph::delete_vertex_attr(g, "GeneID")
    	g <- igraph::delete_vertex_attr(g, "Symbol")
    	E(g)$weight <- 1
    	
    }else if(network %>% stringr::str_detect('KEGGc_')){
    	
    	ls_ig <- oRDS('ig.KEGG.category', placeholder=placeholder, guid=guid, verbose=verbose) %>% tibble::deframe()
    	
    	#oRDS('ig.KEGG.category', placeholder=placeholder, guid=guid, verbose=verbose) %>% arrange(category) %>% pull(category) %>% str_to_title() %>% str_replace_all(" ","") %>% str_c("KEGGc_",.) %>% str_c(collapse='","') %>% noquote()
    	#oRDS('ig.KEGG.category', placeholder=placeholder, guid=guid, verbose=verbose) %>% arrange(category) %>% pull(category) %>% str_to_title() %>% str_c(collapse='","') %>% noquote()
    	
    	default.KEGGc <- c("KEGGc_CellularProcess","KEGGc_EnvironmentalProcess","KEGGc_GeneticProcess","KEGGc_HumanDisease","KEGGc_Metabolism","KEGGc_OrganismalSystem")
    	names(default.KEGGc) <- c("Cellular Process","Environmental Process","Genetic Process","Human Disease","Metabolism","Organismal System")
    	ind <- match(network, default.KEGGc)
		g <- ls_ig[[names(default.KEGGc[ind])]]
		
		g <- igraph::delete_vertex_attr(g, "hsa")
		g <- igraph::delete_vertex_attr(g, "GeneID")
		g <- igraph::delete_vertex_attr(g, "Symbol")
		E(g)$weight <- 1
		
    }else if(network %>% stringr::str_detect('KEGGs_')){
    	
    	ls_ig <- oRDS('ig.KEGG.subcategory', placeholder=placeholder, guid=guid, verbose=verbose) %>% tibble::deframe()
    	
    	#oRDS('ig.KEGG.subcategory', placeholder=placeholder, guid=guid, verbose=verbose) %>% arrange(subcategory) %>% pull(subcategory) %>% str_to_title() %>% str_replace_all(" |:|-|,","") %>% str_c("KEGGs_",.) %>% str_c(collapse='","') %>% noquote()
    	#oRDS('ig.KEGG.subcategory', placeholder=placeholder, guid=guid, verbose=verbose) %>% arrange(subcategory) %>% pull(subcategory) %>% str_c(collapse='","') %>% noquote()
    	
    	default.KEGGs <- c("KEGGs_Aging","KEGGs_AminoAcidMetabolism","KEGGs_BiosynthesisOfOtherSecondaryMetabolites","KEGGs_CancerOverview","KEGGs_CancerSpecificTypes","KEGGs_CarbohydrateMetabolism","KEGGs_CardiovascularDisease","KEGGs_CellGrowthAndDeath","KEGGs_CellMotility","KEGGs_CellularCommunityEukaryotes","KEGGs_CirculatorySystem","KEGGs_DevelopmentAndRegeneration","KEGGs_DigestiveSystem","KEGGs_DrugResistanceAntineoplastic","KEGGs_EndocrineAndMetabolicDisease","KEGGs_EndocrineSystem","KEGGs_EnergyMetabolism","KEGGs_EnvironmentalAdaptation","KEGGs_ExcretorySystem","KEGGs_FoldingSortingAndDegradation","KEGGs_GlycanBiosynthesisAndMetabolism","KEGGs_ImmuneDisease","KEGGs_ImmuneSystem","KEGGs_InfectiousDiseaseBacterial","KEGGs_InfectiousDiseaseParasitic","KEGGs_InfectiousDiseaseViral","KEGGs_LipidMetabolism","KEGGs_MetabolismOfCofactorsAndVitamins","KEGGs_MetabolismOfOtherAminoAcids","KEGGs_MetabolismOfTerpenoidsAndPolyketides","KEGGs_NervousSystem","KEGGs_NeurodegenerativeDisease","KEGGs_NucleotideMetabolism","KEGGs_ReplicationAndRepair","KEGGs_SensorySystem","KEGGs_SignalTransduction","KEGGs_SignalingMoleculesAndInteraction","KEGGs_SubstanceDependence","KEGGs_Translation","KEGGs_TransportAndCatabolism","KEGGs_XenobioticsBiodegradationAndMetabolism")
    	names(default.KEGGs) <- c("Aging","Amino acid metabolism","Biosynthesis of other secondary metabolites","Cancer: overview","Cancer: specific types","Carbohydrate metabolism","Cardiovascular disease","Cell growth and death","Cell motility","Cellular community - eukaryotes","Circulatory system","Development and regeneration","Digestive system","Drug resistance: antineoplastic","Endocrine and metabolic disease","Endocrine system","Energy metabolism","Environmental adaptation","Excretory system","Folding, sorting and degradation","Glycan biosynthesis and metabolism","Immune disease","Immune system","Infectious disease: bacterial","Infectious disease: parasitic","Infectious disease: viral","Lipid metabolism","Metabolism of cofactors and vitamins","Metabolism of other amino acids","Metabolism of terpenoids and polyketides","Nervous system","Neurodegenerative disease","Nucleotide metabolism","Replication and repair","Sensory system","Signal transduction","Signaling molecules and interaction","Substance dependence","Translation","Transport and catabolism","Xenobiotics biodegradation and metabolism")
    	ind <- match(network, default.KEGGs)
		g <- ls_ig[[names(default.KEGGs[ind])]]
		
		g <- igraph::delete_vertex_attr(g, "hsa")
		g <- igraph::delete_vertex_attr(g, "GeneID")
		g <- igraph::delete_vertex_attr(g, "Symbol")
		E(g)$weight <- 1
		
    }else if(network=='HuRI'){
    	g <- oRDS('ig.HuRI', placeholder=placeholder, guid=guid, verbose=verbose)
    	if(weighted==FALSE){
    		E(g)$weight <- 1
    	}

    }else if(network %>% stringr::str_detect('HuRI_')){
		g <- oRDS('ig.HuRI_union', placeholder=placeholder, guid=guid, verbose=verbose)
		## extract relations
		g %>% oIG2TB("edges") -> df_edges
		## restrict to those edges with given source
		flag <- network %>% stringr::str_split("_", simplify=TRUE)
		if(flag[2]=='all'){
			df_edges %>% dplyr::filter_at(dplyr::vars("hiunion","litbm"),dplyr::any_vars(.==1)) -> df_edges
		}else if(flag[2]=='union'){
			df_edges %>% dplyr::filter_at(dplyr::vars("hiunion"),dplyr::any_vars(.==1)) -> df_edges
		}else if(flag[2]=='litbm'){
			df_edges %>% dplyr::filter_at(dplyr::vars("litbm"),dplyr::any_vars(.==1)) -> df_edges
		}
		df_edges %>% dplyr::select_at(dplyr::vars("from","to")) %>% dplyr::mutate(weight=1) -> edges
		
		## extract nodes
		g %>% oIG2TB("nodes") -> nodes

		## "intersected=TRUE": make sure vertices in nodes and edges are the same
		g <- oTB2IG(edges=edges, nodes=nodes, directed=FALSE, intersected=TRUE)		
    }
    
    
    invisible(g)
}

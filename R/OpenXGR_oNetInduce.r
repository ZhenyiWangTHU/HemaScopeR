#' Function to generate a subgraph induced by given vertices and their k nearest neighbors
#'
#' \code{oNetInduce} is supposed to produce a subgraph induced by given vertices and its k nearest neighbors. The input is a graph of "igraph" or "graphNET" object, a list of the vertices of the graph, and a k value for finding k nearest neighbors for these vertices. The output is a subgraph induced by given vertices plus their k neighbours. The resultant subgraph inherits the class from the input one. The induced subgraph contains exactly the vertices of interest, and all the edges between them. 
#'
#' @param ig an object of class "igraph"
#' @param nodes_query the vertices for which the calculation is performed
#' @param knn an integeter specifying how many k steps are used to find the nearest neighbours of the given vertices. By default, knn is set to zero; it means no neighbors will be considered. When knn is 1, the immediate neighbors of the given vertices will be also considered for inducing the subgraph. The same is true when knn is 2, etc
#' @param remove.loops logical to indicate whether the loop edges are to be removed. By default, it sets to false
#' @param largest.comp logical to indicate whether the largest component is only retained. By default, it sets to true for the largest component being left
#' @param min.comp.size an integer specifying the minimum size of component that will be retained. This parameter only works when setting the false to keep the largest component. By default, it sets to 1 meaning all nodes will be retained
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' \itemize{
#'  \item{\code{subg}: an induced subgraph, an object of class "igraph". Appended with a node attribute 'comp' if multiple components are kept}
#' }
#' @note The given vertices plus their k nearest neighbors will be used to induce the subgraph.
#' @export
#' @seealso \code{\link{oNetInduce}}
#' @include oNetInduce.r
#' @examples
#' \dontrun{
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#' V(g)$name <- as.character(V(g))
#'
#' # 2) select the first 10 vertices as the query nodes
#' nodes_query <- as.character(1:4)
#'
#' # 3) produce the induced subgraph only based on the nodes in query
#' subg <- oNetInduce(g, nodes_query, knn=0)
#'
#' # 4) produce the induced subgraph based on the nodes in query are their immediate neighbours
#' subg <- oNetInduce(g, nodes_query, knn=1)
#' }

oNetInduce <- function(ig, nodes_query, knn=0, remove.loops=FALSE, largest.comp=TRUE, min.comp.size=1, verbose=TRUE) 
{
    
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    if(!is(ig,"igraph")){
        warnings("The function must apply to the 'igraph' object.\n")
        return(NULL)
    }
    
    if(is(nodes_query,"igraph.vs")){
        nodes_query <- nodes_query$name
    }
    
    if(is.null(knn)){
        knn <- 0
    }
    knn <- as.integer(knn)
    if(knn < 0){
        knn <- 0
    }
    
    if(is.null(V(ig)$name)){
        V(ig)$name <- as.character(V(ig))
    }
    ## check nodes in query
    ind <- match(nodes_query, V(ig)$name)
    nodes_query <- nodes_query[!is.na(ind)]
    if(length(nodes_query)==0){
    	warnings("Nodes in query cannot be found in the input graph.\n")
        return(NULL)
    }
    
    nei <- unique(unlist(igraph::ego(ig, nodes=nodes_query, order=knn)))
    subg <- igraph::induced_subgraph(ig, vids=nei)
    
    if(remove.loops){
        subg <- igraph::simplify(subg, remove.loops=remove.loops)
    }
    
    if(largest.comp==TRUE){
        clust <- igraph::clusters(subg)
        cid <- which.max(clust$csize)
        subg <- igraph::induced.subgraph(subg, V(subg)[clust$membership==cid])
    }else{
    	min.comp.size <- as.integer(min.comp.size)
    	if(min.comp.size<1 & min.comp.size>igraph::vcount(subg)){
    		min.comp.size <- 1
    	}
    	clust <- igraph::clusters(subg)
		cid <- which(clust$csize >= min.comp.size)
		subg <- igraph::induced_subgraph(subg, V(subg)[clust$membership %in% cid])
		
		## append a node attribute 'comp'
		clust <- igraph::clusters(subg)
		if(0){
			V(subg)$comp <- clust$membership
		}else{
			### sort in a decreasing order
			cid <- clust$csize
			names(cid) <- 1:length(clust$csize)
			cid <- sort(cid,decreasing=T)
			ind <- match(clust$membership, names(cid))
			V(subg)$comp <- ind
		}
    }

####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)

    return(subg)
}
#' Function to convert an igraph into a tibble for nodes or edges
#'
#' \code{oIG2TB} is supposed to convert an igraph into a tibble for nodes or edges.
#'
#' @param ig an "igraph" object
#' @param what what to extract. It can be "edges" for edges and "nodes" for nodes
#' @return
#' a tibble object
#' @note none
#' @export
#' @seealso \code{\link{oIG2TB}}
#' @include oIG2TB.r
#' @examples
#' set.seed(825)
#' ig <- sample_pa(20)
#' V(ig)$name <- seq(1,vcount(ig))
#' ig %>% oIG2TB('edges')
#' ig %>% oIG2TB('nodes')

oIG2TB <- function(ig, what=c('edges','nodes'))
{
    what <- match.arg(what)
    
    name <- name_tmp <- NULL
    
   	if(is(ig,"igraph")){
   		
   		if(what=='edges'){
   			edges <- igraph::as_data_frame(ig, what="edges") %>% tibble::as_tibble()
   			return(edges)
   		}else if(what=='nodes'){
   			if(V(ig)$name %>% duplicated() %>% any()){
   				V(ig)$name_tmp <- V(ig)$name
   				V(ig)$name <- seq_len(igraph::vcount(ig))
   				nodes <- igraph::as_data_frame(ig, what="vertices") %>% tibble::as_tibble() %>% dplyr::select(-name) %>% dplyr::rename(name=name_tmp) %>% dplyr::select(name, dplyr::everything())
   			}else{
   				nodes <- igraph::as_data_frame(ig, what="vertices") %>% tibble::as_tibble()
   			}
   			return(nodes)
   		}
	}else{
		return(NULL)
	}
}

#' Function to visualise an igraph object using ggnetwork
#'
#' \code{oGGnetwork} is supposed to visualise an igraph object using ggnetwork.
#'
#' @param g an object of class "igraph". For an advanced use, it can be a list of igraph objects; in this case, multiple panels will be shown (particularly useful when visualising the same network but color-coded differently)
#' @param node.label either a vector labelling nodes or a character specifying which node attribute used for the labelling. If NULL (by default), no node labelling
#' @param label.wrap.width a positive integer specifying wrap width of node labelling
#' @param label.wrap.lineheight line height spacing for text in ggplot. By default it is 0.8
#' @param node.label.size a character specifying which node attribute used for node label size
#' @param node.label.fontface a character specifying which node attribute used for node label fontface ('plain', 'bold', 'italic', 'bold.italic')
#' @param node.label.color a character specifying which node attribute used for the node label color
#' @param node.label.alpha the 0-1 value specifying transparency of node labelling
#' @param node.label.padding the padding around the labeled node
#' @param node.label.arrow the arrow pointing to the labeled node
#' @param node.label.force the repelling force between overlapping labels
#' @param node.shape an integer specifying node shape or a character specifying which node attribute used for the node shape (no matter whether it is numeric or character)
#' @param node.shape.title a character specifying the title for node shaping
#' @param node.xcoord a vector specifying x coordinates. If NULL, it will be created using igraph::layout_as_tree
#' @param node.ycoord a vector specifying y coordinates. If NULL, it will be created using igraph::layout_as_tree
#' @param node.color a character specifying which node attribute used for node coloring
#' @param node.color.title a character specifying the title for node coloring
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{https://html-color-codes.info/color-names/index.html}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum values for which colors should be plotted
#' @param na.color the color for NAs. By default, it is 'grey80'
#' @param node.color.alpha the 0-1 value specifying transparency of node colors
#' @param node.size either a vector specifying node size or a character specifying which node attribute used for the node size
#' @param node.size.title a character specifying the title for node sizing
#' @param node.size.range the range of actual node size. Can be two values (range) or a value (fixed size)
#' @param slim the minimum and maximum values for which sizes should be plotted
#' @param title a character specifying the title for the plot
#' @param edge.size a numeric value specifying the edge size. By default, it is 0.5. It can be a character specifying which edge attribute defining the edge size (though without the legend)
#' @param edge.color a character specifying the edge color. By default, it is "black". It can be a character specifying which edge attribute defining the edge color (though 
#' @param edge.color.alpha the 0-1 value specifying transparency of edge color. By default, it is 0.5. It can be a character specifying which edge attribute defining the transparency of edge color (though without the legend)
#' @param edge.curve a numeric value specifying the edge curve. 0 for the straight line
#' @param edge.arrow a numeric value specifying the edge arrow. By default, it is 2
#' @param edge.arrow.gap a gap between the arrow and the node
#' @param ncolumns an integer specifying the number of columns for facet_wrap. By defaul, it is NULL (decided on according to the number of groups that will be visualised)
#' @return
#' a ggplot object, appended with 'data_nodes' and 'data_edges'
#' @note none
#' @export
#' @seealso \code{\link{oGGnetwork}}
#' @include oGGnetwork.r
#' @examples
#' \dontrun{
#' ###########################
#' # load REACTOME
#' # restricted to Immune System ('R-HSA-168256') or Signal Transduction ('R-HSA-162582')
#' g <- oRDS('ig.REACTOME')
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="R-HSA-168256", mode="out")
#' nodeInduced <- V(g)[unique(unlist(neighs.out))]$name
#' ig <- igraph::induced.subgraph(g, vids=nodeInduced)
#'
#' # visualise the graph with vertices being color-coded
#' V(ig)$degree <- igraph::degree(ig)
#' gp <- oGGnetwork(g=ig, node.label='term_id', label.wrap.width=30, node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0, node.label.force=1, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='degree', node.color.title='Degree', colormap='grey-orange-darkred', ncolors=64, zlim=c(0,10), node.size.range=3, edge.color="black",edge.color.alpha=0.3,edge.curve=0.05,edge.arrow.gap=0.02, title='')
#' # advanced use: visualise the list of graphs
#' ls_ig <- list(ig, ig)
#' gp <- oGGnetwork(g=ls_ig, node.label='term_id', label.wrap.width=30, node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0, node.label.force=1, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='degree', node.color.title='Degree', colormap='grey-orange-darkred', ncolors=64, zlim=c(0,10), node.size.range=3, edge.color="black",edge.color.alpha=0.3,edge.curve=0.05,edge.arrow.gap=0.02, title='')
#' 
#' ###########################
#' # load PhasedTargets
#' # restricted to disease ('EFO:0000408') or immune system disease ('EFO:0000540')
#' g <- oRDS(RData.customised='ig.PhasedTargets', placeholder=placeholder)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="EFO:0000408", mode="out")
#' nodeInduced <- V(g)[unique(unlist(neighs.out))]$name
#' ig <- igraph::induced.subgraph(g, vids=nodeInduced)
#'
#' # append with the number of approved and phased targets
#' dag <- ig
#' V(dag)$num_approved <- sapply(V(ig)$max_phase,function(x) sum(x$max_phase>=4))
#' V(dag)$num_phased <- sapply(V(ig)$max_phase,function(x) sum(x$max_phase>=0))
#' # keep nodes with num_approved >=20
#' dag_ig <- igraph::induced.subgraph(dag, vids=which(V(dag)$num_approved>=20))
#' # (optional) further restricted to the direct children of the root
#' root <- dnet::dDAGroot(dag_ig)
#' neighs.out <- igraph::neighborhood(dag_ig, order=1, nodes=root, mode="out")
#' nodeInduced <- V(dag_ig)[unique(unlist(neighs.out))]$name
#' dag_ig <- igraph::induced.subgraph(dag_ig, vids=nodeInduced)
#' # nodes colored by num_approved
#' V(dag_ig)$node_color <- log2(V(dag_ig)$num_approved)
#' glayout <- igraph::layout_with_kk(dag_ig)
#' V(dag_ig)$xcoord <- glayout[,1]
#' V(dag_ig)$ycoord <- glayout[,2]
#' gp <- oGGnetwork(g=dag_ig, node.label='term_name', label.wrap.width=30, node.label.size=2, node.label.color='black', node.label.alpha=0.9, node.label.padding=0, node.label.arrow=0, node.label.force=0.5, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='node_color', node.color.title='Approved\n(log2-scale)', colormap='ggplot2.top', ncolors=64, node.size.range=3, edge.color="orange",edge.color.alpha=0.5,edge.curve=0.05,edge.arrow.gap=0.02, title='')
#' 
#' ###########################
#' # visualise gene network
#' glayout <- igraph::layout_with_kk(g)
#' V(g)$xcoord <- glayout[,1]
#' V(g)$ycoord <- glayout[,2]
#' V(g)$degree <- igraph::degree(g)
#' gp <- oGGnetwork(g=g, node.label='name', node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0, node.label.force=0.01, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='priority', node.color.title='5-star\nrating', colormap='yellow-red', ncolors=64, zlim=c(0,5), node.size='degree', node.size.title='Degree', slim=c(0,5), edge.color="orange",edge.color.alpha=0.5,edge.curve=0,edge.arrow.gap=0.025, title='')
#' gp_rating <- oGGnetwork(g=g, node.label='name', node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0.1, node.label.arrow=0, node.label.force=0.01, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='priority', node.color.title='5-star\nrating', colormap='white-yellow-red', ncolors=64, zlim=c(0,5), node.size.range=5, edge.color="orange",edge.color.alpha=0.3,edge.curve=0,edge.arrow.gap=0.02, title='')
#'
#' ###########################
#' # use edge weight to color/size edges (without legends)
#' # edge color
#' #E(g)$color <- oColormap(colormap='RdYlBu', data=E(g)$weight)
#' #E(g)$size <- (E(g)$weight - min(E(g)$weight)) / (max(E(g)$weight) - min(E(g)$weight))
#' e.color <- subset(gp$data, !is.na(na.y))$e.color
#' gp + ggnetwork::geom_edges(color=e.color, show.legend=FALSE)
#' # edge size/thickness
#' e.size <- subset(gp$data, !is.na(na.y))$e.size
#' gp + ggnetwork::geom_edges(size=e.size, show.legend=FALSE)
#' }

oGGnetwork <- function(g, node.label=NULL, label.wrap.width=NULL, label.wrap.lineheight=0.8, node.label.size=NULL, node.label.fontface='plain', node.label.color='black', node.label.alpha=0.90, node.label.padding=1, node.label.arrow=0.01, node.label.force=1, node.shape=19, node.shape.title=NULL, node.xcoord=NULL, node.ycoord=NULL, node.color=NULL, node.color.title=NULL, colormap='grey-orange-darkred', ncolors=64, zlim=NULL, na.color='grey80', node.color.alpha=0.9, node.size=NULL, node.size.title=NULL, node.size.range=c(2,5), slim=NULL, title='', edge.size=0.5, edge.color="steelblue4", edge.color.alpha=0.5, edge.curve=0.1, edge.arrow=2, edge.arrow.gap=0.02, ncolumns=NULL)
{
    
   	if(is(g,"igraph")){
		ls_ig <- list(g)
	}else if(is(g,"list")){
		## Remove null elements in a list
		ls_ig <- base::Filter(base::Negate(is.null), g)
		if(length(ls_ig)==0){
			return(NULL)
		}
	}else{
		stop("The function must apply to 'list' of 'igraph' objects or a 'igraph' object.\n")
	}
	
	## check list_names
	ls_names <- names(ls_ig)
	if(is.null(ls_names)){
		ls_names <- paste('IG', 1:length(ls_ig), sep='_')
		names(ls_ig) <- ls_names
	}
    
	ls_df <- lapply(1:length(ls_ig), function(i){
    	ig <- ls_ig[[i]]
    	    	
		if(igraph::vcount(ig)==0){
			return(NULL)
		}
    	
    	###############################
    	if(is.null(V(ig)$name)){
    		V(ig)$name <- 1:vcount(ig)
    	}
    	###############################
    	
    	########## remove any attributes with the list data type
    	node_attrs <- igraph::vertex_attr_names(ig)
    	for(k in 1:length(node_attrs)){
    		if(is(igraph::vertex_attr(ig, node_attrs[k]),'list')){
    			ig <- ig %>% igraph::delete_vertex_attr(node_attrs[k])
    		}
    	}
    	##########

		nnode <- igraph::vcount(ig)
		## node.xcoord (by default, NULL)
		if(length(node.xcoord)!=nnode | length(node.ycoord)!=nnode){
			if(!is.null(node.xcoord)){
				node.xcoord <- igraph::vertex_attr(ig, node.xcoord)
			}
			if(!is.null(node.ycoord)){
				node.ycoord <- igraph::vertex_attr(ig, node.ycoord)
			}
			
			if(is.null(node.xcoord) | is.null(node.ycoord)){
				## layout
				#glayout <- igraph::layout_with_kk(ig)
				glayout <- igraph::layout_as_tree(ig,root=dnet::dDAGroot(ig),circular=TRUE,flip.y=TRUE)
				if(all(is.na(glayout))){
					glayout <- igraph::layout_with_kk(ig)
				}
				glayout <- glayout[,c(2:1)]
				node.xcoord <- glayout[,1]
				node.ycoord <- glayout[,2]
			}
		}
		## scale into [-1,1]
		if(max(node.xcoord) != min(node.xcoord)){
			node.xcoord <- (node.xcoord - min(node.xcoord)) / (max(node.xcoord) - min(node.xcoord)) * 2 - 1
		}
		if(max(node.ycoord) != min(node.ycoord)){
			node.ycoord <- (node.ycoord - min(node.ycoord)) / (max(node.ycoord) - min(node.ycoord)) * 2 - 1
		}
	
		## node.label (by default, NULL)
		if(length(node.label)!=nnode){
			if(!is.null(node.label)){
				node.label <- igraph::vertex_attr(ig, node.label)
			}
			if(is.null(node.label)){
				node.label <- rep('', nnode)
			}
		}
		## text wrap
		if(!is.null(label.wrap.width)){
			width <- as.integer(label.wrap.width)
			res_list <- lapply(node.label, function(x){
				if(!is.na(x)){
					x <- gsub('_', ' ', x)
					y <- strwrap(x, width=width)
					if(length(y)==2){
						paste(y, collapse='\n')
					}else if(length(y)>2){
						#paste0(y[1], '...')
						paste0(paste(y[1:2],collapse='\n'),'...' )
					}else{
						y
					}
				}else{
					x
				}
			})
			node.label <- unlist(res_list)
		}
		V(ig)$n.label <- node.label
	
		## node.label.size (by default, 0)
		if(length(node.label.size)!=nnode){
			if(!is.null(node.label.size)){
				tmp.node.label.size <- igraph::vertex_attr(ig, node.label.size)
			}else{
				tmp.node.label.size <- rep(0, nnode)
			}
			if(is.null(tmp.node.label.size)){
				node.label.size <- rep(node.label.size, nnode)
			}else{
				node.label.size <- tmp.node.label.size
			}
		}
		V(ig)$n.label.size <- node.label.size
	
		## node.label.fontface (by default, 0)
		if(length(node.label.fontface)!=nnode){
			if(!is.null(node.label.fontface)){
				tmp.node.label.fontface <- igraph::vertex_attr(ig, node.label.fontface)
			}else{
				tmp.node.label.fontface <- rep(0, nnode)
			}
			if(is.null(tmp.node.label.fontface)){
				node.label.fontface <- rep(node.label.fontface, nnode)
			}else{
				node.label.fontface <- tmp.node.label.fontface
			}
		}
		V(ig)$n.label.fontface <- node.label.fontface
	
		## node.label.color (by default, 0)
		if(length(node.label.color)!=nnode){
			if(!is.null(node.label.color)){
				tmp.node.label.color <- igraph::vertex_attr(ig, node.label.color)
			}else{
				tmp.node.label.color <- rep(0, nnode)
			}
			if(is.null(tmp.node.label.color)){
				node.label.color <- rep(node.label.color, nnode)
			}else{
				node.label.color <- tmp.node.label.color
			}
		}
		V(ig)$n.label.color <- node.label.color
	
		## node.color (by default, 0)
		if(length(node.color)!=nnode){
			if(!is.null(node.color)){
				node.color <- igraph::vertex_attr(ig, node.color)
			}
			if(is.null(node.color)){
				node.color <- rep(0, nnode)
			}
		}
		## zlim
		if(is.null(zlim)){
			zlim <- c(min(node.color), max(node.color))
		}
		node.color[node.color<=zlim[1]] <- zlim[1]
		node.color[node.color>=zlim[2]] <- zlim[2]
		V(ig)$n.color <- node.color
	
		## node.size (by default, 0)
		if(length(node.size)!=nnode){
			if(!is.null(node.size)){
				tmp.node.size <- igraph::vertex_attr(ig, node.size)
			}else{
				tmp.node.size <- rep(0, nnode)
			}
			if(is.null(tmp.node.size)){
				node.size <- rep(node.size, nnode)
			}else{
				node.size <- tmp.node.size
			}
		}
		## slim
		if(is.null(slim)){
			slim <- c(min(node.size), max(node.size))
		}
		node.size[node.size<=slim[1]] <- slim[1]
		node.size[node.size>=slim[2]] <- slim[2]
		V(ig)$n.size <- node.size
		
		## node.shape (by default, 19)
		if(length(node.shape)!=nnode){
			if(!is.null(node.shape)){
				tmp.node.shape <- igraph::vertex_attr(ig, node.shape)
			}else{
				tmp.node.shape <- rep(19, nnode)
			}
			if(is.null(tmp.node.shape)){
				node.shape <- rep(node.shape, nnode)
			}else{
				node.shape <- tmp.node.shape
			}
		}
		V(ig)$n.shape <- node.shape
		
		###########################
		nedge <- igraph::ecount(ig)
		## edge.color (by default, 'black')
		if(length(edge.color)!=nedge){
			if(!is.null(edge.color)){
				tmp.edge.color <- igraph::edge_attr(ig, edge.color)
			}else{
				tmp.edge.color <- rep('black', nedge)
			}
			if(is.null(tmp.edge.color)){
				edge.color <- rep(edge.color, nedge)
			}else{
				edge.color <- tmp.edge.color
			}
		}
		E(ig)$e.color <- edge.color
		## edge.color.alpha (by default, 0.5)
		if(length(edge.color.alpha)!=nedge){
			if(!is.null(edge.color.alpha)){
				tmp.edge.color.alpha <- igraph::edge_attr(ig, edge.color.alpha)
			}else{
				tmp.edge.color.alpha <- rep(0.5, nedge)
			}
			if(is.null(tmp.edge.color.alpha)){
				edge.color.alpha <- rep(edge.color.alpha, nedge)
			}else{
				edge.color.alpha <- tmp.edge.color.alpha
			}
		}
		E(ig)$e.color.alpha <- edge.color.alpha
		## edge.size (by default, 0.5)
		if(length(edge.size)!=nedge){
			if(!is.null(edge.size)){
				tmp.edge.size <- igraph::edge_attr(ig, edge.size)
			}else{
				tmp.edge.size <- rep(0.5, nedge)
			}
			if(is.null(tmp.edge.size)){
				edge.size <- rep(edge.size, nedge)
			}else{
				edge.size <- tmp.edge.size
			}
		}
		E(ig)$e.size <- edge.size
		###########################
		
		#gnet <- ggnetwork::ggnetwork(intergraph::asNetwork(ig), layout=cbind(node.xcoord,node.ycoord), arrow.gap=edge.arrow.gap, cell.jitter=0.75)
		gnet <- ggnetwork::ggnetwork(ig, layout=cbind(node.xcoord,node.ycoord), arrow.gap=edge.arrow.gap, cell.jitter=0.75)
		data.frame(gnet, group=rep(names(ls_ig)[i],nrow(gnet)), stringsAsFactors=FALSE)
	})
    df <- do.call(rbind, ls_df)
    
    ## To replace only factors:
    i <- sapply(df, is.factor)
	df[i] <- lapply(df[i], as.character)

    ## ordered according to the input
    df$group <- factor(df$group, levels=names(ls_ig))

    ## make sure numeric
    df$n.color <- as.numeric(df$n.color)
    df$n.size <- as.numeric(df$n.size)

    ## make sure factor    
    df$n.shape <- factor(df$n.shape, levels=sort(unique(df$n.shape)))
    
    if(0){
		n.label.size <- n.label.fontface <- n.label.color <- NULL
		df$n.label.size <- as.numeric(df$n.label.size)
		df$n.label.fontface <- factor(df$n.label.fontface, levels=sort(unique(df$n.label.fontface)))
		df$n.label.color <- factor(df$n.label.color, levels=sort(unique(df$n.label.color)))
    }
    
    #############################################################
    n.color <- n.size <- n.shape <- n.label <- n.label.size <- n.label.color <- NULL
    x <- y <- xend <- yend <- NULL
    
	## ggplot
	gp <- ggplot(df, aes(x=x,y=y,xend=xend,yend=yend)) 
	
	########
	## edges
	########
	whether_node <- NULL
	df$whether_node <- FALSE
	df$whether_node[is.na(df$e.color)] <- TRUE
	e.color <- subset(df, whether_node==FALSE)$e.color
	e.color.alpha <- subset(df, whether_node==FALSE)$e.color.alpha
	e.size <- subset(df, whether_node==FALSE)$e.size
	if(igraph::is_directed(ls_ig[[1]])){
		gp <- gp + ggnetwork::geom_edges(color=e.color, size=e.size,  alpha=e.color.alpha, curvature=edge.curve, arrow=arrow(length=unit(edge.arrow,"pt"),type="closed"), show.legend=FALSE)
	}else{
		gp <- gp + ggnetwork::geom_edges(color=e.color, size=e.size, alpha=e.color.alpha, curvature=edge.curve, show.legend=FALSE)
	}

	########
	## nodes
	########
	if(length(unique(df$n.shape))==1){
		#####################################
		if(!is.numeric(node.shape)){
			node.shape <- 19
		}
		#####################################		
		gp <- gp + ggnetwork::geom_nodes(aes(color=n.color,size=n.size), shape=node.shape, alpha=node.color.alpha)
	}else{
		gp <- gp + ggnetwork::geom_nodes(aes(color=n.color,size=n.size, shape=n.shape), alpha=node.color.alpha)
		#######
		# node shape
		gp <- gp + scale_shape(guide=guide_legend(node.shape.title,title.position="top",ncol=1))
		#######
	}
	
	if(is.null(node.color.title)){
		node.color.title <- 'Node color'
	}
	if(is.null(zlim)){
		zlim <- range(df$n.color[!is.na(df$n.color)])
	}
	if(zlim[1] != zlim[2]){
		#gp <- gp + scale_colour_gradientn(colors=oColormap(colormap)(ncolors), limits=zlim, breaks=c(zlim[1],sum(zlim)/2,zlim[2]), guide=guide_colorbar(title=node.color.title,title.position="top",barwidth=0.5,nbin=64,draw.ulim=FALSE,draw.llim=FALSE), na.value=na.color)
		gp <- gp + scale_colour_gradientn(colors=oColormap(colormap)(ncolors), limits=zlim, guide=guide_colorbar(title=node.color.title,title.position="top",barwidth=0.5,nbin=64,draw.ulim=FALSE,draw.llim=FALSE), na.value=na.color)
	}else{
		gp <- gp + scale_colour_gradientn(colors=oColormap(colormap)(ncolors), guide=guide_colorbar(title=node.color.title,title.position="top",barwidth=0.5))
	}
	
	if(is.null(node.size.title)){
		node.size.title <- 'Node size'
	}
	if(is.null(slim)){
		slim <- range(df$n.size)
	}
	if(slim[1] != slim[2]){
		gp <- gp + scale_size_continuous(limits=slim, range=node.size.range, guide=guide_legend(node.size.title,title.position="top",ncol=1))
	}else{
		gp <- gp + scale_size_continuous(limits=slim, range=node.size.range, guide=guide_legend(node.size.title,title.position="top",ncol=1))
	}
	
	gp <- gp + ggnetwork::theme_blank()
	gp <- gp + theme(text=element_text(family="sans")) + labs(title=title) + theme(plot.title=element_text(hjust=0.5,size=10,face="bold"), plot.margin=unit(rep(0,4),rep("lines",4)))
	
	if((zlim[1]!=zlim[2]) & (slim[1]!=slim[2])){
		gp <- gp + theme(legend.position="right")
    }else{
		if(slim[1]==slim[2]){
			gp <- gp + guides(size="none")
		}
		if(zlim[1]==zlim[2]){
			gp <- gp + guides(color="none")
		}
		
		if(length(unique(df$n.shape))==1){
			gp <- gp + guides(shape="none")
		}
		
    }
    gp <- gp + theme(legend.title=element_text(size=8,face="bold"),legend.text=element_text(size=6))
    
    ## facet by 'group' artificially added 
    if(length(ls_ig)>1){
		if(is.null(ncolumns)){
			ncolumns <- ceiling(sqrt(length(ls_ig)))
		}
		group <- NULL
		gp <- gp + facet_wrap(~group, ncol=ncolumns)
		gp <- gp + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=12,face="bold"), strip.placement="inside", panel.spacing=unit(0,"lines"))
    }
    
	if(sum(df$n.label=='')!=nrow(df)){
	
		###########
		StatNodes <- ggplot2::ggproto("StatNodes", ggplot2::Stat,
					compute_layer = function(data, scales, params) {
						if(all(c("xend", "yend") %in% names(data))) {
							unique(subset(data, select = c(-xend, -yend)))
						} else {
							unique(data)
						}
					}
				)
		## my_geom_nodetext_repel
		my_geom_nodetext_repel <- function (mapping = NULL, data = NULL, parse = FALSE, ..., box.padding = unit(0.25, "lines"), point.padding = unit(1e-06, "lines"), segment.size = 0.5, arrow = NULL, force = 1, max.iter = 2000, nudge_x = 0, nudge_y = 0, na.rm = FALSE, show.legend = FALSE, inherit.aes = TRUE){
			ggplot2::layer(data = data, mapping = mapping, stat = StatNodes, 
				geom = ggrepel::GeomTextRepel, position = "identity", 
				show.legend = show.legend, inherit.aes = inherit.aes, 
				params = list(parse = parse, na.rm = na.rm, box.padding = box.padding, 
					point.padding = point.padding, 
					segment.size = segment.size, arrow = arrow, force = force, 
					max.iter = max.iter, nudge_x = nudge_x, nudge_y = nudge_y, max.overlaps=Inf,
					...)
			)
		}
		###########	
		n.label.size <- n.label.fontface <- n.label.color <- NULL
		###
		# not working properly
		###
		n.label.size <- subset(df, whether_node==TRUE)$n.label.size
		n.label.fontface <- subset(df, whether_node==TRUE)$n.label.fontface
		n.label.color <- subset(df, whether_node==TRUE)$n.label.color
		###
		
		gp <- gp + my_geom_nodetext_repel(data=df, aes(label=n.label), lineheight=label.wrap.lineheight, size=n.label.size, color=n.label.color, fontface=n.label.fontface, alpha=node.label.alpha, box.padding=unit(0.5,"lines"), point.padding=unit(node.label.padding,"lines"), segment.alpha=0.2, segment.size=0.2, arrow=arrow(length=unit(node.label.arrow,'npc')), force=node.label.force)
	}
    
    ####################
    
    if(1){
		df <- gp$data
		
		## append data_nodes
		df_sub <- subset(df, whether_node==TRUE)
		ind <- match(colnames(df_sub), c('whether_node','e.color','edge.color','e.size'))
		df_sub <- df_sub[,is.na(ind)]
		df_sub <- df_sub[!duplicated(df_sub),]
		gp$data_nodes <- df_sub
		
		## append data_edges
		df_sub <- subset(df, whether_node==FALSE)
		ind <- match(colnames(df_sub), c('x','y','xend','yend','e.color','edge.color','e.size'))
		df_sub <- df_sub[,!is.na(ind)]
		df_sub <- df_sub[!duplicated(df_sub),]
		gp$data_edges <- df_sub
    }
    
    invisible(gp)
}

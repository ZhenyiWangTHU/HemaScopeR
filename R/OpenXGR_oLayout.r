#' Function to define graph node coordinates according to igraph- or sna-style layout
#'
#' \code{oLayout} is supposed to define graph node coordinates according to igraph- or sna-style layout.
#'
#' @param g an object of class "igraph" (or "graphNEL") for a graph
#' @param layout a character specifying graph layout function. This character can be used to indicate igraph-style layout ("layout_nicely","layout_randomly","layout_in_circle","layout_on_sphere","layout_with_fr","layout_with_kk","layout_as_tree","layout_with_lgl","layout_with_graphopt","layout_with_sugiyama","layout_with_dh","layout_with_drl","layout_with_gem","layout_with_mds","layout_as_bipartite"), or sna-style layout ("gplot.layout.adj","gplot.layout.circle","gplot.layout.circrand","gplot.layout.eigen","gplot.layout.fruchtermanreingold","gplot.layout.geodist","gplot.layout.hall","gplot.layout.kamadakawai","gplot.layout.mds","gplot.layout.princoord","gplot.layout.random","gplot.layout.rmds","gplot.layout.segeo","gplot.layout.seham","gplot.layout.spring","gplot.layout.springrepulse","gplot.layout.target"), or graphlayouts-style layout ("graphlayouts.layout_with_stress","graphlayouts.layout_as_backbone"), or ForeceAtlas2 layout used in Dephi ("gephi.forceatlas2")
#' @param seed an integer specifying the seed
#' @param flip logical to indicate whether x- and y-coordiates flip. By default, it sets to false
#' @return It returns an igraph object, appended by node attributes including "xcoord" for x-coordinates, "ycoord" for y-coordiates.
#' @export
#' @seealso \code{\link{oGGnetwork}}
#' @include oLayout.r
#' @examples
#' \dontrun{
#' # load REACTOME
#' # restricted to Immune System ('R-HSA-168256') or Signal Transduction ('R-HSA-162582')
#' g <- oRDS('ig.REACTOME', placeholder=placeholder)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="R-HSA-168256", mode="out")
#' nodeInduced <- V(g)[unique(unlist(neighs.out))]$name
#' ig <- igraph::induced.subgraph(g, vids=nodeInduced)
#' 
#' # compare Fruchterman and Reingold force-directed placement algorithm
#' ## based on igraph layout
#' ig1 <- oLayout(ig, layout="layout_with_fr")
#' gp1 <- oGGnetwork(ig1, node.xcoord="xcoord", node.ycoord="ycoord")
#' ## based on sna layout
#' ig2 <- oLayout(ig, layout="gplot.layout.fruchtermanreingold")
#' gp2 <- oGGnetwork(ig2, node.xcoord="xcoord", node.ycoord="ycoord")
#' 
#' # compare Kamada-Kawai force-directed placement algorithm
#' ## based on igraph layout
#' ig1 <- oLayout(ig, layout="layout_with_kk")
#' gp1 <- oGGnetwork(ig1, node.xcoord="xcoord", node.ycoord="ycoord")
#' ## based on sna layout
#' ig2 <- oLayout(ig, layout="gplot.layout.kamadakawai")
#' gp2 <- oGGnetwork(ig2, node.xcoord="xcoord", node.ycoord="ycoord")
#' ## do together
#' layouts <- c("layout_with_fr","gplot.layout.fruchtermanreingold","layout_with_kk","gplot.layout.kamadakawai", "gephi.forceatlas2")
#' ls_ig <- lapply(layouts, function(x) oLayout(ig, layout=x))
#' names(ls_ig) <- layouts
#' gp <- oGGnetwork(ls_ig, node.xcoord='xcoord', node.ycoord='ycoord', ncolumns=5)
#' }

oLayout <- function(g, layout=c("layout_nicely","layout_randomly","layout_in_circle","layout_on_sphere","layout_with_fr","layout_with_kk","layout_as_tree","layout_with_lgl","layout_with_graphopt","layout_with_sugiyama","layout_with_dh","layout_with_drl","layout_with_gem","layout_with_mds","layout_as_bipartite", "gplot.layout.adj","gplot.layout.circle","gplot.layout.circrand","gplot.layout.eigen","gplot.layout.fruchtermanreingold","gplot.layout.geodist","gplot.layout.hall","gplot.layout.kamadakawai","gplot.layout.mds","gplot.layout.princoord","gplot.layout.random","gplot.layout.rmds","gplot.layout.segeo","gplot.layout.seham","gplot.layout.spring","gplot.layout.springrepulse","gplot.layout.target", "graphlayouts.layout_with_stress","graphlayouts.layout_as_backbone", "gephi.forceatlas2"), seed=825, flip=F)
{
    
    layout <- layout[1]
    
    if(is(g,"graphNEL")){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if(!is(ig,"igraph")){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
	
	glayout <- NULL
	
	if(grepl("gplot",layout)){
		## based on sna
		m <- as.matrix(oConverter(ig, from='igraph', to='dgCMatrix', verbose=F))
		set.seed(seed)
		eval(parse(text=paste0('glayout <- sna::',layout,'(m, NULL)')))
		
	}else if(grepl("graphlayouts",layout)){
		## based on graphlayouts
		layout <- gsub('graphlayouts.','',layout)
		set.seed(seed)
		eval(parse(text=paste0('glayout <- graphlayouts::',layout,'(ig)')))
		
		if(layout=='layout_as_backbone'){
			glayout <- glayout$xy
		}
		
	}else if(layout=='gephi.forceatlas2'){
		#################################################
		## ForceAtlas2, a Continuous Graph Layout Algorithm for Handy Network Visualization Designed for the Gephi Software
		iterations = 100 #Number of iterations to be performed
		linlog = FALSE #If TRUE the algorithm uses logarithmic attraction force `F <- log (1+F)`
		pos = NULL #If NULL the initial positions are random (that is, the initial locations of points, number of nodes X dimension)
		nohubs = FALSE #If TRUE nodes with high indegree have more central position than nodes with outdegree (for directed graphs)
        k = 400 #the repel constant: the greater the constant k the stronger the repulsion force between points
        gravity = 1 #the gravity constant: indicates how strongly the nodes should be attracted to the center of gravity
        ks = 0.1 #the speed constant: the greater the value of `ks` the more movement the nodes make under the acting forces
        ksmax = 10 #limits the speed
        delta = 1 #modify attraction force, the weights are raised to the power of `delta`
        center = NULL #the center of gravity
        tolerance = 0.1 #the tolerance to swinging constant
        dim = 2 #dimension of the positions

	  	if(length(names(igraph::edge.attributes(ig)))==0){
			attr_g <- NULL
	  	}else{
			attr_g <- names(igraph::edge.attributes(ig))
	  	}
	  	A <- igraph::get.adjacency(ig, type="both", attr=attr_g, edges=F, names=T, sparse=F)
		
		## center of gravity is by default set to the origin
	  	if(is.null(center)) center <- rep(0,dim)

	  	nnodes <- nrow(A)
	  	## Binary will be a matrix of simple incidence (0-not connected, 1-connected)
	  	Binary <- A
	  	Binary[Binary!=0] <- 1
	  	## Deg will be a vector of the degrees of vertices
	  	Deg <- rowSums(Binary)
	  	#### Forces1 will be a table containing all the sums of forces acting on points at a step
	  	Forces1 <- matrix(0, nrow=dim, ncol=nnodes)

	  	## If there are no initial coordinates of points,
	  	## they are chosen at random from 1000^dim square
	  	if(is.null(pos)){
			difference <- 2000/(nnodes*dim)
			set.seed(seed)
			position <- matrix(sample(seq(-1000,1000,difference),nnodes*dim),nnodes,dim)
	  	}else{
			position <- pos
	  	}

	  	## None of the nodes should be exactly at the center of gravity
	  	temp <- which(position[,1] == center[1])
	  	for(index in 2:ncol(position)){
			temp <- intersect(temp,which(position[,index] == center[index]))
	  	}
	  	position[temp,] <- center + 0.01
	  	#rm(index,temp)

	  	## displacement will be a matrix of points' movement at the current iteration
	  	displacement <- matrix(0, nrow=dim, ncol=nnodes)
	  	m <- nrow(position)
	  	
	  	for(iteration in 1:iterations){
			displacement <- displacement * 0
			#### Forces2 is the table of the forces from previous step
			#### Forces1 is the table of the forces from current step
			Forces2 <- Forces1
			Forces1 <- matrix(, nrow = dim, ncol = 0)

			#### Calculate the Forces for each node
			### Distance matrix between all nodes
			distances <- as.matrix(stats::dist(position))
			distances[which(distances < 0.01)] <- 0.01 #We impose a minimum distance
			### Each element of the list contains a matrix with the j = 1,2,..., dim dimension of the unitary vector 1
			mylist <- vector("list",dim)
			for (j in 1:dim){
				mylist[[j]] <- (tcrossprod(position[,j],rep(1,m))-tcrossprod(rep(1,m),position[,j]))/distances
			}
			### Calculate the repulsion Force
			Fr <- k*((tcrossprod(rep(1,m),Deg)+1)*(tcrossprod(Deg,rep(1,m))+1))/distances

			#The classical attraction force is just based on distance
			Fa <- distances
			#The linlog mode calculates the attraction force as log(1+d(n1,n2))
			if(linlog){
				Fa <- log(1+Fa)
			}
			#Edge weights. The edges are weighted based on parameter delta. delta=0 implies no weight
			Fa <- (A^delta)*Fa

			#Dissuade Hubs. This mode is meant to grant authorities (nodes with high indegree)
			#a more central position than hubs (nodes with high outdegree)
			if(nohubs){
			  	Fa <- Fa/(tcrossprod(Deg,rep(1,m))+1)
			}

			### Function to calculate the Attraction and Repulsion forces
			Farfunction <- function(x) rowSums(x*(Fr-Fa),na.rm=T)
			### And we aggregate it over all dimensions
			Far <- do.call(rbind,lapply(mylist,Farfunction))
			### Unitary Vector 2, the directions between each point and the center
			uv2 <- apply(matrix(rep(center,m),nrow=m,byrow=T)-position,1,function(x) x/sqrt(sum(x^2)))
			### The gravity force
			#Fg <- uv2*matrix(rep(gravity*(rowSums(A)+1),dim),nrow=dim,byrow=T)
			Fg <- uv2*matrix(rep(gravity*(Deg+1),dim),nrow=dim,byrow=T)
			### Forces 1 is the sum between all forces: Far (Fa + Fr) and Fg
			Forces1 <- Far+Fg
			Forces1 <- round(Forces1,2) #Use the first two decimals for the Forces.

			#### Swing is the vector of the swingings of all points
			swing <- abs(colSums((Forces1-Forces2)^2)^(1/2))
			Global_swing <- sum((Deg + 1)*swing)

			#If the swing of all nodes is zero, then convergence is reached and we break.
			if(all(swing==0)){
				message(sprintf("Convergence reached at step %d (%s)", iteration, as.character(Sys.time())), appendLF=T)
			  	break
			}

			#### tra is the vector of the traction of all points
			tra <- abs(colSums((Forces1+Forces2)^2)^(1/2))/2
			Global_tra <- sum((Deg+1)*tra)

			#### Global speed calculation
			Global_speed <- tolerance * Global_tra/Global_swing
			#### speed is the vector of individual speeds of points
			speed <- ks * Global_speed /  (1 + Global_speed * (swing)^(1/2))

			#### Imposing constrains on speed
			speed_constrain <- ksmax/abs(colSums((Forces1^2))^(1/2))
			speed <- ifelse(speed>=speed_constrain,speed_constrain,speed)

			#### calculating displacement and final position of points after iteration
			displacement <- Forces1 * t(matrix(rep(speed,dim),nnodes,dim))
			position <- position + t(displacement)
		}
	  	glayout <- position
	  	#################################################
	  
	}else{
		## based on igraph
		set.seed(seed)
		eval(parse(text=paste0('glayout <- ',layout,'(ig)')))
	}
	
	if(!is.null(glayout)){
		## scale into [-1,1]
		
		### whether flip coordinates (vertically or horizontally)
		if(flip){
			node.xcoord <- glayout[,2]
			node.ycoord <- glayout[,1]		
		}else{
			node.xcoord <- glayout[,1]
			node.ycoord <- glayout[,2]
		}
		
		if(max(node.xcoord) != min(node.xcoord)){
			node.xcoord <- (node.xcoord - min(node.xcoord)) / (max(node.xcoord) - min(node.xcoord)) * 2 - 1
		}
		if(max(node.ycoord) != min(node.ycoord)){
			node.ycoord <- (node.ycoord - min(node.ycoord)) / (max(node.ycoord) - min(node.ycoord)) * 2 - 1
		}
		glayout <- cbind(node.xcoord, node.ycoord)
		
		## append 'xcoord' and 'ycoord' into ig
		V(ig)$xcoord <- glayout[,1]
		V(ig)$ycoord <- glayout[,2]
	}
    
    invisible(ig)
}



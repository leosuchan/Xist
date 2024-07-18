###########################################################################################################################
# R code supplementary to the paper "Distributional limits of graph cuts on discretized grids" by Suchan, Li, Munk (2024) #
###########################################################################################################################


install.packages(c("igraph", "ggplot2", "latex2exp", "mvtnorm", "RANN"))
# Required packages:
# igraph: Required for essentially everything; igraph supplies the graph structure and many useful functions.
# ggplot2: Required for plotting only.
# latex2exp: Used for better ggplot labels only.
# mvtnorm: Used for efficient generation of samples from a multivariate normal distribution; only used in the function rnorm_multinomial_mvtnorm
#           which can be substituted by rnorm_multinomial_cholesky, ..._iterative or ..._matricized (see below).
# RANN: Used for its k-nearest neighbour function RANN::nn2; only required if generate_sample_graph is used with the nearest_neighbors option.


### PART I: FUNCTION DEFINITIONS ###


library(igraph)
# Implements the graph structure and provides an algorithm to compute st-MinCuts.
library(ggplot2)
# Used for most of the plotting.
library(latex2exp)
# Enables greek letters and other fancy plot labels.

# NOTE:
# Descriptions of the following self-implemented functions are given immediately *after* the function.
# The inner workings of the functions are sparsely commented to explain the way these functions operate.
# It is recommended to use an IDE that allows folding of {...} environments to be able to fold all function
# definitions for clarity and easier future referencing.

generate_grid_graph <- function(sample, t=1, sigm=NA, weights_by_sample=TRUE){
  
  if(t < 1) simpleError("t must be greater or equal than 1")
  if(length(dim(sample)) != 2) simpleError("sample must be 2-dimensional")
  n <- length(sample)
  if(is.vector(sample)) sample <- matrix(sample, round(sqrt(n)))
  sample <- t(sample)
  ndim <- dim(sample)
  # Determining (or guessing) the appropriate grid size.
  
  edgevec <- weightvec <- x1coordvec <- x2coordvec <- labelvec <- c()
  degreevec <- rep(0, n)
  # Initializing the attribute vectors that will be assigned to the graph.
  
  for(x1 in 1:ndim[1]){
    
    for(x2 in 1:ndim[2]){
      
      x1coordvec <- c(x1coordvec, x1)
      x2coordvec <- c(x2coordvec, x2)
      # Keeping track of x1 and x2 coordinates for later.
      
      labelvec <- c(labelvec, toString(x1 + (x2-1)*ndim[2]))
      # Keeping track of labels for later.
      
      if(sample[x1,x2] > 0){
      # There will be no connections if the sample at (x1,x2) is 0.
      
        y1vec <- x1 + (-floor(t)):floor(t)
        # First coordinate tries to connect to all left and right vertices within distance t.
        
        for(y1 in y1vec[y1vec >= 1 & y1vec <= ndim[1]]){
          
          y2vec <- c(x2 + 0:floor(sqrt(t^2 - (x1-y1)^2)))
          # To remove double connections in the second coordinate require that y2vec >= x2.
          
          for(y2 in y2vec[y2vec <= ndim[2]]){
          # yvec1 and yvec2 are all vertices that (x1,x2) will be connected to in a radius of t (within sample).
            
            if((x2 != y2 || x1 < y1) && sample[y1,y2] > 0){
            # If x2 = y2 require that x1 < y1 to remove double connections in the first coordinate.
            # Also, there will be no edge if the sample at (y1,y2) is 0.
              
              edgevec <- c(edgevec, x2 + (x1-1)*ndim[2], y2 + (y1-1)*ndim[2])
              # Adding the edge between x and y. Note that matrices in R are collapsed/indexed by columns.
              
              if(is.na(sigm) && weights_by_sample) current_weight <- sample[x1,x2] * sample[y1,y2]
              else if(is.na(sigm)) current_weight <- 1
              else if(weights_by_sample) current_weight <- sample[x1,x2] * sample[y1,y2] * exp(-sqrt((x1-y1)^2 + (x2-y2)^2)/sigm)
              else current_weight <- exp(-sqrt((x1-y1)^2 + (x2-y2)^2)/sigm)
              # Assigns weights according to the combination of sigm and weights_by_sample
              
              weightvec <- c(weightvec, current_weight)
              degreevec[x2 + (x1-1)*ndim[2]] <- degreevec[x2 + (x1-1)*ndim[2]] + current_weight
              degreevec[y2 + (y1-1)*ndim[2]] <- degreevec[y2 + (y1-1)*ndim[2]] + current_weight
              # Updates the degree of the graph (to avoid having to compute it later).
              
            }
          }
        }
      }
    }
  }
  
  g <- graph(edgevec, n = n, directed = FALSE)
  # Constructs an undirected graph with n vertices (none removed!) out of edgevec.
  
  E(g)$weight <- weightvec
  # Assigns capacities (=weights) of the edges from edge_list. This is needed for igraph::min_cut.
  
  V(g)$x <- x1coordvec
  V(g)$y <- x2coordvec
  V(g)$label <- labelvec
  V(g)$sample <- as.vector(t(sample))
  V(g)$degree <- degreevec
  # Assigns coordinates, labels (likely unnecessary), underlying distribution and degree to the graph.
  
  return(simplify(delete_vertices(g, which(V(g)$degree == 0))))
  
}
# Generates a graph out of a sample formatted as a matrix (e.g. grey values of an image) using t-neighborhood,
# with weights given by the product of the incident vertices (if weights_by_sample) and an exponentially decaying distance factor sigm.
# Input: Matrix sample, neighborhood radius t, exponential distance scaling parameter sigm, and
#        weights_by_sample indicates whether edge weights should be scaled by the sample values at the incident vertices.
# Output: igraph graph g.

generate_sample_graph <- function(sample, t, sigm=NA, nearest_neighbors=NA){
  
  d <- dim(sample)[2]
  n <- dim(sample)[1]
  # Determines sample size n and number of dimensions d.
  
  weightvec <- edge_mat <- c()
  degreevec <- rep(0,n)
  # Initializing. edge_mat will be a matrix containing all edges and weightvec will contain edge weights
  
  if(!is.na(nearest_neighbors) && nearest_neighbors >= 1) knn_sample <- RANN::nn2(sample, k=floor(nearest_neighbors)+1)
  # Computes the k nearest neighbors of each point and their respective distances.
  # The nearest point is always the point itself which will be ignored (therefore k+1, not k).
  
  if(is.na(t)){
    if(is.na(nearest_neighbors) || nearest_neighbors < 1) t <- median(RANN::nn2(sample, k=2)$nn.dists[,2])
    else t <- median(knn_sample$nn.dists[,2])
    print(paste("Choosing t =",t))
  }
  # Automatically chooses t if not specified. May lead to a disconnected graph; not recommended.
  
  for(x in 1:n){
    
    if(!is.na(nearest_neighbors) && nearest_neighbors >= 1){
    
      for(i in 2:dim(knn_sample$nn.idx)[2]){
      # Ignore the first entry of knn_sample since this is always the point i itself
        
        i_idx <- knn_sample$nn.idx[x,i]
        # Selects the (index of the) i-1-th nearest neighbor of x.
        
        if(sum((sample[x,] - sample[i_idx,])^2) > t){
        # Only add an edge here if it would not have been added through t-neighborhood later.
        
          edge_mat <- rbind(edge_mat, c(min(x, i_idx), max(x, i_idx)))
          # Adds an edge between x and i_idx. min/max to avoid counting edges double.
          
          current_weight <- ifelse(is.na(sigm), 1, exp(-sqrt(sum((sample[x,] - sample[i_idx,])^2))/sigm))
          weightvec <- c(weightvec, current_weight)
          degreevec[x] <- degreevec[x] + current_weight
          degreevec[i_idx] <- degreevec[i_idx] + current_weight
          # Should the graph be built using gaussian weights? If so, sigm is the parameter.
          # If not, the weights are given by w_{ij} = 1.
        
        }
        
      }
      
    }
    
    if(x < n){
    
      for(y in (x+1):n){
      # Tries not to add edges twice by requiring that x < y.
        
        if(sqrt(sum((sample[x,] - sample[y,])^2) <= t)){
          # Are nodes x and y connected?
          
          edge_mat <- rbind(edge_mat, c(min(x, y), max(x, y)))
          # Adds an edge between x and y. min/max are just for safety (since always x<y).
          
          current_weight <- ifelse(is.na(sigm), 1, exp(-sqrt(sum((sample[x,] - sample[y,])^2))/sigm))
          weightvec <- c(weightvec, current_weight)
          degreevec[x] <- degreevec[x] + current_weight
          degreevec[y] <- degreevec[y] + current_weight
          # Should the graph be built using gaussian weights? If so, sigm is the parameter.
          # If not, the weights are given by w_{xy} = 1.
          
        }
      }
    }
  }
  
  g <- graph(edges=c(t(edge_mat)), n=n, directed = FALSE)
  # Creates a graph with n vertices and only the unique edges from edge_mat.
  
  E(g)$weight <- weightvec
  E(g)$capacity <- weightvec
  # Assigns capacities (=weights) of the edges from edge_list. This is needed for igraph::min_cut.
  
  V(g)$x <- sample[,1]
  if(d > 1) V(g)$y <- sample[,2]
  if(d > 2) V(g)$z <- sample[,3]
  # Assigns the first three components (if applicable) as coordinates for plotting
  
  V(g)$label <- sapply(1:n, toString)
  V(g)$sample <- rep(1, n)
  # Assigns labels to all vertices.
  
  V(g)$degree <- degreevec
  # Computes the degree (i.e. the sum of the weights of all neighbouring vertices) for each vertex.
  
  return(simplify(g))
  # Returns the graph g, removing multiedges that appear for mutual k-nearest neighbors.
}
# Generates a graph out of a sample (matrix rows are sample point coordinates) using t-neighborhood and nearest neighbors.
# Input: n x d matrix sample containing n d-dimensional points, neighborhood radius t, exponential distance scaling parameter sigm,
#        and nearest_neighbours (if integer) indicates whether edges should (also) be defined between nearest_neighbors-nearest neighbors.
# Output: igraph graph g.

xcut_via_stmincut <- function(g, progress_bar=TRUE){
  
  mcut_min <- rcut_min <- ncut_min <- ccut_min <- Inf
  mcut <- rcut <- ncut <- ccut <- c()
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = length(V(g))*(length(V(g))+1)/2, style = 3)
  # Initializes a progress bar of appropriate length.
  
  for(i in 1:(length(V(g))-1)){
    
    for(j in (i+1):length(V(g))){
      
      S_ij <- min_cut(g, source = i, target = j, capacity = E(g)$weight, value.only = FALSE)$partition1
      # Computes a partitions attaining the ij-MinCut.
      
      ij_xcut <- xcut_partition_S(g, S_ij)
      # Computes the XCut of the partition S_ij.
      
      if(!is.na(ij_xcut$mcut_min) & ij_xcut$mcut_min < mcut_min){
        mcut_min <- ij_xcut$mcut_min
        mcut <- S_ij
      }
      if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < rcut_min){
        rcut_min <- ij_xcut$rcut_min
        rcut <- S_ij
      }
      if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ncut_min){
        ncut_min <- ij_xcut$ncut_min
        ncut <- S_ij
      }
      if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ccut_min){
        ccut_min <- ij_xcut$ccut_min
        ccut <- S_ij
      }
      # Updates the best XCut and the partition that attains it
      
      if(progress_bar) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
      # Updates progress bar.
    }
  }
  
  if(progress_bar) close(progress)
  
  return(list(mcut=mcut, rcut=rcut, ncut=ncut, ccut=ccut,
              mcut_min=mcut_min, rcut_min=rcut_min, ncut_min=ncut_min, ccut_min=ccut_min))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}
# Realizes the Xvst (XCut via st-MinCuts) algorithm.
# Input: igraph graph g, boolean progress_bar whether to display a progress bar.
# Output: List containing the XCut values and corresponding partitions for MinCut, Ratio Cut, NCut and Cheeger Cut balancing terms.

xcut_exact <- function(g, progress_bar=TRUE, normalization=FALSE, degree_up_to_date=FALSE){
  
  mcut_min <- rcut_min <- ncut_min <- ccut_min <- Inf
  mcut <- rcut <- ncut <- ccut <- list()
  # Keeps track of the minimum respective XCut. Again, apologies for ugliness. Should be efficient, though.
  
  vvec <- V(g)[.inc(E(g))]
  n_eff <- length(vvec)
  # Computes the set of vertices with at least one edge attached to it.
  
  if(n_eff < length(V(g)) && !all(sapply(vvec$sample, function(a)is.null(a) || is.na(a) || a==0))){
    simpleWarning("Isolated vertices detected. All graph cuts are zero.")
    return(list(mcut=vvec, rcut=vvec, ncut=vvec, ccut=vvec,
                mcut_min=0, rcut_min=0, ncut_min=0, ccut_min=0))
  }
  # Terminates if vvec != V(g) and at least one isolated vertex has weight (conveyed through sample), i.e. cannot be ignored.
  
  if(n_eff > 30) simpleWarning("Warning: More than 30 nodes. Computation may take very long (exponential complexity).")
  # Hard-coded check whether the number of nodes isn't too large in order to keep computations manageable.
  
  if(!degree_up_to_date || any(is.null(V(g)$degree))) V(g)$degree <- sapply(V(g), function(i)sum(E(g)[.inc(i)]$weight))
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = 2^(n_eff-1)-1, style = 3)
  # Initializes a progress bar of appropriate length.
  
  for(k in 1:(2^(n_eff-1)-1)){#c(7,27,56,232,3392)){#
    # Iterates over all partitions; k, converted to binary, is the representation of the partition.
    
    S <- as.logical(intToBits(k)[1:n_eff])
    # Converts a number k into a partition S by binary encoding.
    
    S_C <- !S
    # S_C is the complement of S in g.
    
    S_vol <- sum(V(g)[S]$degree)
    S_C_vol <- sum(V(g)[S_C]$degree)
    # Computes the volumes of S and S_C using vertex degrees.
    
    mcut_min_S <- sum(E(g)[S %--% S_C]$weight) / ifelse(normalization, sum(E(g)$weight), 1)
    rcut_min_S <- mcut_min_S / (sum(S) * sum(S_C)) * ifelse(normalization, length(V(g))^2, 1)
    ncut_min_S <- mcut_min_S / (S_vol * S_C_vol) * ifelse(normalization, (S_vol + S_C_vol)^2, 1)
    ccut_min_S <- mcut_min_S / min(S_vol, S_C_vol) * ifelse(normalization, (S_vol + S_C_vol), 1)
    # Computes the various XCuts of S.
    
    if(mcut_min_S == mcut_min) mcut <- c(mcut, list(S))
    else if(mcut_min_S < mcut_min){
      mcut_min <- mcut_min_S
      mcut <- list(S)
    }
    if(rcut_min_S == rcut_min) rcut <- c(rcut, list(S))
    else if(rcut_min_S < rcut_min){
      rcut_min <- rcut_min_S
      rcut <- list(S)
    }
    if(ncut_min_S == ncut_min) ncut <- c(ncut, list(S))
    else if(ncut_min_S < ncut_min){
      ncut_min <- ncut_min_S
      ncut <- list(S)
    }
    if(ccut_min_S == ccut_min) ccut <- c(ccut, list(S))
    else if(ccut_min_S < ccut_min){
      ccut_min <- ccut_min_S
      ccut <- list(S)
    }
    # Updates the optimal XCuts and the partition attaining them.
    
    if(progress_bar) setTxtProgressBar(progress, k)
    
  }
  return(list(mcut=mcut, rcut=rcut, ncut=ncut, ccut=ccut,
              mcut_min=mcut_min, rcut_min=rcut_min, ncut_min=ncut_min, ccut_min=ccut_min))
}
# Computes the exact XCut value of a (small) graph by brute force (i.e. exponential complexity). May take a long time.
# Input: igraph graph g, boolean progress_bar whether to display a progress bar, boolean normalization of whether to normalize
#         the cut values to allow comparison across different graphs, boolean degree_up_to_date of whether V(g)$degree contains
#         the correct degrees of the vertices.
# Output: List containing the XCut values and corresponding partitions for MinCut, Ratio Cut, NCut and Cheeger Cut balancing terms.

xcut_partition_S <- function(g, S, normalize=FALSE, degree_up_to_date=FALSE){
  
  if(is.logical(S)) S <- which(S)
  S_C <- V(g)[-S]
  # S_C is the complement of S in g.
  
  if(degree_up_to_date && !any(is.null(V(g)$degree))){
    
    S_vol <- sum(V(g)[S]$degree)
    S_C_vol <- sum(V(g)[S_C]$degree)
    # Computes the volumes of S and S_C using vertex degrees if these are present and up to date.
    
  } else {
    
    S_vol <- sum(E(g)[.inc(S)]$weight) + sum(E(g)[S %--% S]$weight)
    S_C_vol <- sum(E(g)[.inc(S_C)]$weight) + sum(E(g)[S_C %--% S_C]$weight)
    # Computes the volumes of S and S_C if no vertex degrees are present or they are not up to date.
    
  }
  
  mcut_min <- sum(E(g)[S %--% S_C]$weight) / ifelse(normalize, sum(E(g)$weight), 1)
  rcut_min <- mcut_min / (length(S) * length(S_C)) * ifelse(normalize, length(V(g))^2, 1)
  ncut_min <- mcut_min / (S_vol * S_C_vol) * ifelse(normalize, (S_vol + S_C_vol)^2, 1)
  ccut_min <- mcut_min / min(S_vol, S_C_vol) * ifelse(normalize, (S_vol + S_C_vol), 1)
  # Computes the various XCuts of S.
  
  return(list(mcut_min=mcut_min, rcut_min=rcut_min, ncut_min=ncut_min, ccut_min=ccut_min))
  # Returns a list of the computed XCut values of S.
}
# Computes the XCut values (MinCut, Ratio Cut, NCut and Cheeger Cut) for a graph g and partition S.
# normalize indicates whether normalization is desired to allow for XCut value comparison across different graphs.
# If enabled, degree_up_to_date uses the vertex attribute degree instead of computing it manually, thus saving time.
# Input: igraph graph g, subset of vertices S, boolean normalize, boolean degree_up_to_date.
# Output: List of XCut values (MinCut, Ratio Cut, NCut and Cheeger Cut) of partition S on g.

plot_graph_partition <- function(g, partition=NULL, vertex_scale=400/sqrt(length(V(g))), axes=FALSE, edge_scale=20,
                                 vertex_offset=8/log(length(V(g))), edge_offset=1, main=NULL,
                                 xlim=range(V(g)$x), ylim=range(V(g)$y)){
  
  if(!is.null(partition)){
  # If no partition is given, the graph will be plotted without colors (or with pre-assigned colors).
    
    if(length(partition) < length(V(g))){
    # Under this condition it is assumed that partition is a list/vector consituting a subset of vertices. 
      
      V(g)$color <- 'green3'
      V(g)[unlist(partition)]$color <- 'blue'
      # Assigns colors to the cut.
      
    } else if(length(partition) == length(V(g))){
      
      if(is.numeric(partition)){
      # If partition only contains numbers, greyscale colors are automatically assigned.
        
        allcolors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
        V(g)$color <- c("green3", "blue", sample(allcolors, length(unique(partition))-2))[partition]
        # Assigns colors to partition according to a greyscale vector (and green and blue for no reason).
        
      } else V(g)$color <- partition
      # Otherwise, it is assumed that partition is a vector of color assignments.
      
    }
  }
  
  if(!is.null(E(g)$weight)) edge_weights <- E(g)$weight
  else edge_weights <- 1
  if(!is.null(V(g)$sample)) vertex_weights <- V(g)$sample
  else vertex_weights <- 1
  # Determines vertex and edge weights for plotting.
  
  par(cex=0.7, mai=c(0,0.1,0,0.1))
  # Sets the plot boundaries, adjusting the spacing.
  
  if(axes){
  # Plots with x and y axes if desired. Then, the exact positions of the vertices are preserved.
    
    plot(g, rescale=FALSE, axes=TRUE, xlim=xlim, ylim=ylim, asp = 0, vertex.label=NA,
         edge.width = edge_scale*edge_weights/max(edge_weights)+edge_offset,
         vertex.size = vertex_scale*vertex_weights/max(vertex_weights)+vertex_offset, main=main)
    
  } else {
    
    plot(g, vertex.label=NA, edge.width = edge_scale*edge_weights/max(edge_weights)+edge_offset,
         vertex.size = vertex_scale*vertex_weights/max(vertex_weights)+vertex_offset, main=main)
    
  }
  # Plots the graph with appropriately scaled vertices and edges, where size emphasizes weight.
  
}
# Plots a graph, with partition colored appropriately and vertices/edges scaled proportionally to their weight.
# Input: igraph graph g, partition being either a subset of vertices or a V(g)-length vector of either numbers or colors,
#        boolean axes whether to draw axes, and plotting parameters vertex_offset, vertex_scale, edge_offset, edge_scale, xlim, ylim, main.
# Output: none. Plots are generated.

xist <- function(g, locmax_by_degree=FALSE, normalize=FALSE, progress_bar=TRUE){
  
  if(is.null(E(g)$weight)) E(g)$weight <- 1
  
  n <- length(V(g))
  
  mcut_min <- rcut_min <- ncut_min <- ccut_min <- Inf
  mcut <- rcut <- ncut <- ccut <- c()
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  if(locmax_by_degree){
  # Under locmax_by_degree, Vloc is constituted by vertices whose degree is higher than that of all of its neighbors.
    
    if(is.null(V(g)$degree)) V(g)$degree <- sapply(1:n, function(i)sum(E(g)[.inc(i)]$weight))
    # Computes the vertex degree if it has not been computed before.
    
    Vloc <- which(sapply(1:n, function(i)all(sapply(V(g)[.nei(i)], function(j)V(g)$degree[i] >= V(g)$degree[j]))))
    # Computes Vloc, the vertex set of local maxima
  
  } else {
  # Otherwise, Vloc is constituted by vertices whose sample (i.e. vertex weight) is higher than that of all of its neighbors.
    
    if(is.null(V(g)$sample)) V(g)$sample <- 1
    # However, if no sample values are provided, then all vertices are assigned sample value 1, s.t. Vloc = V.
    # Note that in this case, Xist yields the same output as Xvst if the XCut minimizing partition attains its st-MinCut uniquely (see paper).
    
    Vloc <- which(sapply(1:n, function(i)all(sapply(V(g)[.nei(i)], function(j)V(g)$sample[i] >= V(g)$sample[j]))))
    # Computes Vloc as local maxima through vertex weights.
    
  }
  # Computes the set of local maxima of g.
  
  N <- length(Vloc)
  
  eqvec <- rep(1, N)
  # eqvec is tau from the paper, initialized with all 1's.
  
  if(N == 1){
    
    simpleWarning("Only one local maximum computed; no clusters were computed. Returning infinte cut values.")
    return(list(mcut = c(), mcut_min = Inf, rcut = c(), rcut_min = Inf, ncut = c(), ncut_min = Inf, ccut = c(), ccut_min = Inf))
    # If there is just one local maximum, Xist terminates without further computation - g should not be divided in this case.
    
  }
  
  if(progress_bar) progress <- txtProgressBar(min = 1, max = N-1, style = 3)
  # Updates the progress bar.
  
  for(s in 2:N){
    
    t <- eqvec[s]
    # eqvec determines the next Vloc vertex index t to compute the st-MinCut with.
    
    ij_cut <- min_cut(g, source = Vloc[s], target = Vloc[t], capacity = E(g)$weight, value.only = FALSE)
    
    S_ij <- ij_cut$partition1
    
    if(!(Vloc[s] %in% S_ij)) S_ij  <- V(g)[-S_ij]
    # Adjusts S_ij such that it contains Vloc[s].
    
    for(i in s:N) if(Vloc[i] %in% S_ij && eqvec[i] == t) eqvec[i] <- s
    # Associate s as well as any i that has not yet been part of a MinCut (but is on the same side as s) to s.
    # This serves to prevent the computation of an ij-MinCut where i and j have been on different sides of a previous st-MinCut.
    
    ij_xcut <- xcut_partition_S(g, S_ij, normalize=normalize)
    # Computes the XCut of the partition S_ij.
    
    if(!is.na(ij_xcut$mcut_min) & ij_xcut$mcut_min < mcut_min){
      mcut_min <- ij_xcut$mcut_min
      mcut <- S_ij
    }
    if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < rcut_min){
      rcut_min <- ij_xcut$rcut_min
      rcut <- S_ij
    }
    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ncut_min){
      ncut_min <- ij_xcut$ncut_min
      ncut <- S_ij
    }
    if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ccut_min){
      ccut_min <- ij_xcut$ccut_min
      ccut <- S_ij
    }
    # Updates the respective current best XCut value and corresponding partition.
    
    if(progress_bar) setTxtProgressBar(progress, s)
    # Updates progress bar.
    
  }
  
  if(progress_bar) close(progress)
  
  return(list(mcut=mcut, rcut=rcut, ncut=ncut, ccut=ccut,
              mcut_min=mcut_min, rcut_min=rcut_min, ncut_min=ncut_min, ccut_min=ccut_min))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
  
}
# Realizes the Xist (XCut imitation through st-MinCuts) algorithm. locmax_by_degree indicates whether Vloc should be determined
# by degree as opposed to sample size (i.e. vertex weight), normalize indicates whether normalization is desired
# to allow for XCut value comparison across different graphs, progress_bar enables a progress bar.
# Input: igraph graph g, booleans locmax_by_degree, normalize, progress_bar.
# Output: List of XCut values (MinCut, Ratio Cut, NCut and Cheeger Cut) of partition S on g.

xcut_multipartition_S <- function(g, partitions, normalize=FALSE, degree_up_to_date=FALSE){
  
  partitions <- as.integer(factor(partitions))
  # Assuming that partitions is a V(g)-length vector containing partition assignments, standardize these from 1 to k.
  
  k <- max(partitions)
  # Determines the number k of desired partitions.
  
  if(k==1) return(list(mcut_min=Inf, rcut_min=Inf, ncut_min=Inf, ccut_min=Inf))
  # If every vertex belongs to the same partition, there is nothing to cut.
  
  mcut_min <- rcut_min <- ncut_min <- ccut_min <- 0
  # Initializes the cut values.
  
  for(i in 1:k){
    
    S <- which(partitions == i)
    S_C <- V(g)[-S]
    # For each cluster label i, determines the corresponding partition S and its complement S_C
    
    if(degree_up_to_date && !any(is.null(V(g)$degree))){
      
      S_vol <- sum(V(g)[S]$degree)
      S_C_vol <- sum(V(g)[S_C]$degree)
      # Computes the volumes of S and S_C using vertex degrees if these are present and up to date.
      
    } else {
      
      S_vol <- sum(E(g)[.inc(S)]$weight) + sum(E(g)[S %--% S]$weight)
      S_C_vol <- sum(E(g)[.inc(S_C)]$weight) + sum(E(g)[S_C %--% S_C]$weight)
      # Computes the volumes of S and S_C if no vertex degrees are present or they are not up to date.
      
    }
    # Computing the volumes of S and S_C (also computable using sum(E(g)[.inc(S)]$weight) if vertex degree attribute is missing)
    
    if(normalize) edgesum <- sum(E(g)$weight)
    # The total sum of all weights is only necessary if normalization is desired.
    
    mc_S <- sum(E(g)[S %--% S_C]$weight) / ifelse(normalize, edgesum, 1)
    # Computes the MinCut value of S, accounting for normalization through normalize.
    
    mcut_min <- mcut_min + mc_S
    rcut_min <- rcut_min + mc_S / (length(S) * length(S_C)) * ifelse(normalize, length(V(g))^2, 1)
    ncut_min <- ncut_min + mc_S / (S_vol * S_C_vol) * ifelse(normalize, (S_vol + S_C_vol)^2, 1)
    ccut_min <- ccut_min + mc_S / min(S_vol, S_C_vol) * ifelse(normalize, (S_vol + S_C_vol), 1)
    # Computes the various XCuts of S.
    
  }
  
  return(list(mcut_min=mcut_min/2, rcut_min=rcut_min/2, ncut_min=ncut_min/2, ccut_min=ccut_min/2))
  # Returns a list of the computed XCut values of S.
}
# Computes the XCut values (MinCut, Ratio Cut, NCut and Cheeger Cut) for a graph g and multipartition partitions.
# normalize indicates whether normalization is desired to allow for XCut value comparison across different graphs.
# If enabled, degree_up_to_date uses the vertex attribute degree instead of computing it manually, thus saving time.
# Input: igraph graph g, V(g)-length vector of partition assignments partitions, boolean normalize, boolean degree_up_to_date.
# Output: List of XCut values (MinCut, Ratio Cut, NCut and Cheeger Cut) of multipartition partitions on g.

multixist <- function(g, k, cut="ncut", locmax_by_degree=FALSE, normalize=TRUE, colorvec=NULL){
  
  if(cut %in% c("mcut", "rcut", "ncut", "ccut")) simpleError("'cut' must be either one of 'mcut', 'rcut', 'ncut' or 'ccut'")
  # If cut doesn't match any of the allowed outputs, terminate.
  
  labels_list <- list(V(g))
  partitions_list <- list(rep(1,length(V(g))))
  index_best <- 1
  # Initializes the list of partitions labels_list that g is currently divided into,
  # partitions_list containing partitions from all iterations of the algorithm, and index_best is the partition about to be cut.
  
  cut_list <- list(xist(g, locmax_by_degree=locmax_by_degree, normalize=normalize, progress_bar=FALSE))
  # Initializes the list containing all cuts, starts with a simple Xist call.
  
  glist <- list(g)
  # Initializes the list containing all (sub-)graphs.
  
  for(i in 2:k){
    
    g_current <- glist[[index_best]]
    partition_current <- cut_list[[index_best]][[cut]]
    # Selects the current (cut-wise best) graph g_current and corresponding cut partition_current.
    
    labels_list[[i]] <- labels_list[[index_best]][-partition_current]
    labels_list[[index_best]] <- labels_list[[index_best]][partition_current]
    # Stores the original labels corresponding to both sides of the cut on g_current in labels_list, i.e. updates it.
    
    glist[[i]] <- delete_vertices(g_current, V(g_current)[partition_current])
    glist[[index_best]] <- delete_vertices(g_current, V(g_current)[-partition_current])
    # Updates the respective graph by removing all vertices that have been moved.
    
    cut_list[[i]] <- xist(glist[[i]], locmax_by_degree=locmax_by_degree, normalize=normalize, progress_bar=FALSE)
    cut_list[[index_best]] <- xist(glist[[index_best]], locmax_by_degree=locmax_by_degree, normalize=normalize, progress_bar=FALSE)
    # Applies to Xist to both of the updated graphs.
    
    partitions_list[[i]] <- partitions_list[[i-1]]
    partitions_list[[i]][labels_list[[i]]] <- i
    # Stores the current global partition assignment, updating it according to the current labels.
    
    index_best <- which.min(sapply(1:i, function(j)cut_list[[j]][[paste0(cut,"_min")]]))
    # Selects the next best index index_best, i.e. the one which minimizes the XCut.
    
  }
  
  return(partitions_list)
  # Returns the list of clustering assignments.
  
}
# Realizes the Multi-Xist algorithm, i.e. clusters g into k partitions using Xist iteratively.
# Input: igraph graph g, number of clusters k, cut specifying which balanced graph cut to use, locmax_by_degree
#        indicating whether to construct Vloc using vertex degree or vertex weight (i.e. through the underlying sample),
#        normalize indicating whether to normalize the cut values for comparison across different graphs (passed to xcut_partition_S).
# Output: List of partitions containing clustering assignments for each of the k iterations of Multi-Xist.

rnorm_xcut_limit <- function(size, g, partitions=NULL, force_unequal_volumes=FALSE){
  
  pvec <- V(g)$sample
  m <- length(V(g))
  
  if(is.null(partitions)){
    xcut_g <- xcut_exact(g, progress_bar=FALSE)
    partitions <- unique(c(xcut_g$mcut, xcut_g$rcut, xcut_g$ncut, xcut_g$ccut))
    partitions_which <- lapply(list(xcut_g$mcut, xcut_g$rcut, xcut_g$ncut, xcut_g$ccut),
                               function(ls)sapply(ls, function(S)which(sapply(partitions, function(T)all(T==S)))))
    partitions_null <- TRUE
    partitions <- lapply(partitions, which)
#    print(partitions)
    
  } else partitions_null <- FALSE
  # If partitions==NULL, determine which partitions attain the respective XCuts (avoiding double computations),
  # storing which partitions belong to which cuts in partitions_which
  
  M <- length(partitions)
  
  Z <- rnorm_multinomial_matricized(size, pvec)
  # Generates a Gaussian sample with multinomial covariance matrix according to the underlying probability vector pvec.
  
  xcut_limit_samples <- list(mcut_sample=matrix(nrow=size,ncol=M), rcut_sample=matrix(nrow=size,ncol=M),
                             ncut_sample=matrix(nrow=size,ncol=M), ccut_sample=matrix(nrow=size,ncol=M),
                             mult_sample=matrix(Z, nrow=size))
  # Initializes the list of arrays where the XCut distributions (and that of Z) are stored.
  
  for(k in 1:M){
  # Iterates through all partitions in partitions.
    
    if(length(partitions[[k]]) < m){
      S <- partitions[[k]]
      S_boolean <- rep(FALSE, m)
      S_boolean[S] <- TRUE
    } else {
      S_boolean <- partitions[[k]]
      S <- which(S_boolean)
    }
    S_C <- V(g)[-S]
    # Defines the partition S and its complement S_C as well as a corresponding boolean assignment S_boolean.
    
#    xc_S <- xcut_partition_S(g, S)
    q_S <- sapply(1:m, function(i)ifelse(S_boolean[i], sum(E(g)[i %--% S_C]$weight), sum(E(g)[i %--% S]$weight))/pvec[i])
    vol_S <- sum(E(g)[.inc(S)]$weight) + sum(E(g)[S %--% S]$weight)
    vol_S_C <- sum(E(g)[.inc(S_C)]$weight) + sum(E(g)[S_C %--% S_C]$weight)
    mcut_S <- sum(q_S[S] * pvec[S])
    ncut_S <- mcut_S / (vol_S * vol_S_C)
    ccut_S <- mcut_S / min(vol_S, vol_S_C)
    # Computes the quantities necessary for computation of the limit distribution (see paper for notation).
    
    coeff_ncut <- (q_S - (q_S * ifelse(S_boolean,1,-1) * (vol_S - vol_S_C) + 2 * ifelse(S_boolean,vol_S_C,vol_S)
                          * sapply(1:m, function(i)sum(pvec[V(g)[.nei(i)]]))) * ncut_S) / (vol_S * vol_S_C)
    # Computes the coefficients for the NCut limit distribution.
    
    vol_comparator <- round(vol_S - vol_S_C, max(0,round(2*log10(min(pvec)^2)))+5)
    # Determines whether vol_S == vol_S_C or not. This construction hopefully avoids volume misclassification as unequal due to rounding errors.
    
    if(force_unequal_volumes && vol_comparator==0) vol_comparator <- 1
    # If volumes are equal, but this is behavior is not wanted, the wrong limit distribution will be computed by e.g. setting vol_comparator to 1.
    # This may be desired for illustration purposes, i.e. to emphasize the importance of considering the case vol_S == vol_S_C separately.
    
    if(vol_comparator > 0){
      S_pZ <- !S_boolean
      vol_S_pZ <- vol_S_C
    }
    else if(vol_comparator < 0){
      S_pZ <- S_boolean
      vol_S_pZ <- vol_S
    }
    # If vol(S) != vol(S^C), defines S_pZ to be the partition with minimum volume of the two.
    
    for(j in 1:size){
      
      xcut_limit_samples$mcut_sample[j,k] <- sum(q_S * Z[j,])
      xcut_limit_samples$rcut_sample[j,k] <- sum(q_S / (length(S) * (m - length(S))) * Z[j,])
      xcut_limit_samples$ncut_sample[j,k] <- sum(coeff_ncut * Z[j,])
      # Computes a sample of the MCut, RCut and NCut limit distributions, respectively.
      
      if(vol_comparator == 0){
        
        vol_decider <- function(i) sum(Z[j,i] * pvec[V(g)[.nei(i)]]) + sum(Z[j,V(g)[.nei(i)]] * pvec[i])
        
        if(sum(sapply(S, vol_decider)) < sum(sapply(S_C, vol_decider))){
          S_pZ <- S_boolean
          vol_S_pZ <- vol_S
        } else {
          S_pZ <- !S_boolean
          vol_S_pZ <- vol_S_C
        }
        # In case of equal volumes, the sum of vol_decider over S and S_C determines whether S_pZ will be S or S_C.
        
      }
      coeff_ccut <- q_S - ccut_S * sapply(1:m, function(i)ifelse(S_pZ[i], sum(pvec[V(g)[.nei(i)]]), 0) + sum(pvec[V(g)[.nei(i)]] * S_pZ[V(g)[.nei(i)]]))
      
      xcut_limit_samples$ccut_sample[j,k] <- sum(coeff_ccut / vol_S_pZ * Z[j,])
      # With S_pZ defined, this computes the sample of the CCut limit distribution.
      
    }
    
  }
  
  if(partitions_null){
    
    cut_names <- c("mcut_sample", "rcut_sample", "ncut_sample", "ccut_sample")
#    print(partitions_which)
#    print(lapply(xcut_limit_samples[-5], head))
    xcut_limit_samples_exact <- lapply(1:length(cut_names), function(cn)sapply(1:size, function(j)min(xcut_limit_samples[[cut_names[cn]]][j, partitions_which[[cn]]])))
    # Computes the minimum of the partitions that attain the true XCut (partitions_which).
    
    names(xcut_limit_samples_exact) <- cut_names
    return(c(xcut_limit_samples_exact, list(mult_sample=matrix(Z, nrow=size))))
    
  } else return(xcut_limit_samples)
  
}
# Generates a sample from the XCut limit distribution from the underlying graph g.
# Imput: Sample size size, graph g, list of partitions partitions dictating what partitions the limit distribution should consider simultaneously,
#        (if partitions=NULL, the exact XCut is taken instead), boolean force_unequal_volumes determining whether to compute the wrong limit distribution of XC_S(g) if vol(S)!=vol(S^C) for a partition S.
# Output: List of (1.-4.) samples of the limit distributions of the MCut, RCut, NCut and CCut statistics,
#         and of (5.) the corresponding N(0,Sigma) sample for the multinomial covariance matrix Sigma.

rnorm_xcut_empirical <- function(size, g, n, t=1, partitions=NULL, progress_bar=FALSE, g_truth=g, degree_up_to_date=FALSE){
  
  pvec <- V(g)$sample
  m <- length(V(g))
  M <- max(1, length(partitions))
  nrows <- length(unique(V(g)$x))
  cutnames <- c("mcut","rcut","ncut","ccut")
  # Initializes some necessary quantities.
  
  if(is.null(partitions) && m > 30) simpleWarning("Exact XCut computation not advisable for graphs with >30 vertices.")
  # Throws a warning if computation of the exact XCut statistics are desired and the graph is too large for this to be feasible.
  
  xcut_limit_samples <- list(mcut_sample=matrix(nrow=size,ncol=M), rcut_sample=matrix(nrow=size,ncol=M),
                             ncut_sample=matrix(nrow=size,ncol=M), ccut_sample=matrix(nrow=size,ncol=M),
                             mult_sample=matrix(nrow=size,ncol=m))
  # Initializes the list of arrays where the XCut distributions (and that of Z) are stored.
  
  if(!is.null(partitions)){
    
    xcuts_all_truth <- tapply(unlist(lapply(partitions, function(S)xcut_partition_S(g_truth, S)),FALSE,FALSE),
                              rep(cutnames,M), function(x)unlist(x,FALSE,FALSE))
    # Computes the (joint) XCut values of the partitions from partitions for the underlying distribution as given by g_truth.
    
  } else xcuts_all_truth <- xcut_exact(g_truth, progress_bar=FALSE, degree_up_to_date=degree_up_to_date)
  #If no partitions are given, the exact XCut values of the underlying graph g_truth are computed.
  
  if(progress_bar) progress <- txtProgressBar(min = 1, max = size, style = 3)
  
  for(i in 1:size){
    
    multinom_sample <- rmultinom(1, n, pvec) / n
#    print(multinom_sample)
    g_sample <- generate_grid_graph(matrix(multinom_sample, nrows),t=t)
#    print(xcut_exact(g_sample, progress_bar=FALSE))
    # Generates a graph from a sample of the multinomial distribution generated from pvec.
    
    if(is.null(partitions)){
      
      xcuts_all <- mapply(function(x,y)sqrt(n)*(x-y), xcut_exact(g_sample, progress_bar=FALSE, degree_up_to_date=degree_up_to_date)[paste0(cutnames,"_min")],
                                                      xcuts_all_truth[paste0(cutnames,"_min")], SIMPLIFY=FALSE)
#      print(xcuts_all)
      for(cn in cutnames) xcut_limit_samples[[paste0(cn,"_sample")]][i,1] <- xcuts_all[[paste0(cn,"_min")]]
      # If no partition is given, the empirical distribution of the exact XCut statistic is computed.
      
    } else for(j in 1:M){
      
      S <- partitions[[j]]
      
      xcuts_S <- mapply(function(x,y)sqrt(n)*(x-y), xcut_partition_S(g_sample, S), lapply(xcuts_all_truth, function(x)x[j]), SIMPLIFY=FALSE)
      for(cn in cutnames) xcut_limit_samples[[paste0(cn,"_sample")]][i,j] <- xcuts_S[[paste0(cn,"_min")]]
      # Saves the XCut statistics of partition S.
      
    }
    
    xcut_limit_samples$mult_sample[i,] <- sqrt(n) * (multinom_sample - V(g_truth)$sample)
    # Saves the multinomial sample statistic used in the XCut computations.
    
    if(progress_bar) setTxtProgressBar(progress, i)
    
  }
  
  if(progress_bar) close(progress)
  
  return(xcut_limit_samples)
  
}
# Generates a sample from the empirical distribution of the XCut statistic sqrt(n)*(XC(g_sample)- XC(g_truth)) for the sample graph g_sample.
# Input: Sample size size, graph g, number of observations n for each sample, list of partitions partitions dictating what partitions the empirical distribution should consider simultaneously
#        - if partitions=NULL, the exact XCut is taken instead -, and boolean progress_bar determining whether to draw a progress bar.
# Output: List of (1.-4.) samples of the empirical distributions of the MCut, RCut, NCut and CCut statistics,
#         and of (5.) the corresponding N(0,Sigma) sample for the multinomial covariance matrix Sigma.

rnorm_multinomial_mvtnorm <- function(n, pvec) mvtnorm::rmvnorm(n, sigma = diag(pvec) - outer(pvec,pvec))
rnorm_multinomial_cholesky <- function(n, pvec){
  d <- length(pvec)
  qvec <- c(1, 1 - cumsum(pvec))
  cholesky_sigma <- outer(1:d, 1:d, function(i,j)ifelse(i>j, -pvec[i]*sqrt(pvec[j]/qvec[j]/qvec[j+1]),
                                                        ifelse(i<j, 0, sqrt(pvec[i]*qvec[i+1]/qvec[i]))))
  return(t(replicate(n, CSmat_manual %*% rnorm(d), simplify=T)))
}
rnorm_multinomial_iterative <- function(n, pvec){
  d <- length(qvec)
  qvec <- c(1, 1 - cumsum(pvec))
  hvec <- pvec / qvec[-(d+1)]
  dvec <- hvec * qvec[-1]
  genfs <- function(){
    fs <- c()
    for(i in 1:d) fs[i] <- - sum(fs[1:(i-1)]) * hvec[i] + sqrt(dvec[i]) * rnorm(1)
    return(fs)
  }
  return(t(replicate(n, genfs())))
}
rnorm_multinomial_matricized <- function(n, pvec){
  d <- length(pvec)
  qvec <- c(1, 1 - cumsum(pvec))
  hvec <- pvec / qvec[-(d+1)]
  dvec <- hvec * qvec[-1]
  fs <- matrix(sqrt(dvec[1]) * rnorm(n), n)
  sumvec <- 0
  for(i in 2:d){
    sumvec <- sumvec + fs[,i-1]
    fs <- cbind(fs, sqrt(dvec[i]) * rnorm(n) - sumvec * hvec[i])
  }
  return(fs)
}
# Four methods of generating a sample from the limit distribution of a multinomial vector.
# For small length(pvec), _mvtnorm is best; otherwise choose _matricized.
# Input: Sample size n, vector of probabilities pvec.
# Output: A sample of N(0,Sigma) of size n, where Sigma is the covariance matrix of the multinomial distribution Mult(1, pvec).

generate_qq_plots_xcut_limit_distribution <- function(BigM, g, nvec, colorvec=NULL, progress_bar=TRUE, legend_pos=c(0.75,0.25)){
  
  cutnames <- c("mcut","rcut","ncut","ccut")
  pvec <- V(g)$sample
  
  if(progress_bar) cat("Phase 1/3: Computing the exact XCut\n")
  
  g_xcut_exact <- xcut_exact(g, progress_bar=progress_bar)
  
  if(progress_bar) cat("\nPhase 2/3: Computing the asymptotic distributions\n")
  
  asympsample_xcuts <- sapply(cutnames, function(cn)rnorm_xcut_limit(BigM, g, g_xcut_exact[[cn]])[[paste0(cn,"_sample")]], simplify=FALSE, USE.NAMES=TRUE)
  asympsample_xcuts_min <- sapply(cutnames, function(cn)sort(apply(asympsample_xcuts[[cn]], 1, min)), simplify=FALSE, USE.NAMES=TRUE)
  # Computes the XCut limit distribution according to the explicit expressions derived in the paper.
  
  asympsample_xcuts_unique <- sapply(cutnames, function(cn)rnorm_xcut_limit(BigM, g, g_xcut_exact[[cn]][1])[[paste0(cn,"_sample")]], simplify=FALSE, USE.NAMES=TRUE)
  # Computes the XCut limit distribution under the erroneous assumption that only one (the first) partition in g2_xcut_exact attains the global XCut minimum.
  
  asympsample_xcuts_uneqvol <- sapply(cutnames, function(cn)rnorm_xcut_limit(BigM, g, g_xcut_exact[[cn]], force_unequal_volumes=TRUE)[[paste0(cn,"_sample")]], simplify=FALSE, USE.NAMES=TRUE)
  asympsample_xcuts_min_uneqvol <- sapply(cutnames, function(cn)apply(asympsample_xcuts_uneqvol[[cn]], 1, min), simplify=FALSE, USE.NAMES=TRUE)
  # Computes the XCut limit distribution under the erroneous assumption that vol(S) != vol(S^C) for all XCut attaining partitions S.
  
  if(progress_bar) cat(paste("---> Done!\nPhase 3/3: Computing the empirical XCut distributions - has", length(nvec), "subphases\n"))
  
  empsample_raw <- lapply(nvec, function(n)rnorm_xcut_empirical(BigM, g, n, progress_bar=progress_bar))
  empsample_ordered <- sapply(cutnames, function(cn)sapply(1:length(nvec), function(i)sort(empsample_raw[[i]][[paste0(cn,"_sample")]], decreasing=FALSE)), simplify=FALSE, USE.NAMES=TRUE)
  # Generates a sample from the empirical XCut distribution for each n, orders them, and sorts them by cut.
  
  empsample_df <- data.frame(theoretical=c(sapply(cutnames,function(cn)rep(asympsample_xcuts_min[[cn]],length(nvec)))),
                             sample=c(sapply(empsample_ordered,c)),
                             n=sapply(rep(rep(nvec,each=BigM),length(cutnames)),toString), cut=rep(cutnames,each=BigM*length(nvec)))
  # Compiles the ordered samples into a dataframe in an appropriate form for plotting.
  
  qq_line_params <- function(sample, theoretical){
    y <- quantile(sample[!is.na(sample)], c(0.25, 0.75))
    x <- quantile(theoretical[!is.na(theoretical)], c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1] - slope * x[1]
    return(as.numeric(c(slope, int)))
  }
  # Manually calculates the parameters of a quantile-quantile line going through the 0.25 and the 0.75 quantile of
  # 'sample' against (a list of quantiles or a sample of) the theoretical distribution 'theoretical'.
  
  qq_line_unique <- sapply(cutnames, function(cn)qq_line_params(asympsample_xcuts_unique[[cn]], asympsample_xcuts_min[[cn]]), simplify=FALSE, USE.NAMES=TRUE)
  qq_line_uneqvol <- sapply(cutnames, function(cn)qq_line_params(asympsample_xcuts_min_uneqvol[[cn]], asympsample_xcuts_min[[cn]]), simplify=FALSE, USE.NAMES=TRUE)
  # Quantile-quantile line connecting the respective 0.25 and 0.75 quantiles of the sample ncnormals and (a sample from)
  # the true limiting distribution ncut_true_sample.
  
  if(is.null(colorvec)) colorvec <- c("green4", "lawngreen", "yellow", "darkgoldenrod", "chocolate1", "orangered1")
  # Defines a vector of colors if not already given.
  
  qq_plots <- list()
  # Initializes the list where all quantile-quantile plots will be stored.
  
  for(i in 1:3){
    
    cn <- c("mcut", "rcut", "ncut")[i]
    
    qq_plots[[paste0(cn,"_plot")]] <- ggplot(empsample_df[empsample_df$cut==cn,]) + geom_point(aes(x=theoretical, y=sample, color=n), key_glyph="point") +
      geom_abline(mapping=aes(slope=1,intercept=0, colour="black")) +
      geom_abline(mapping=aes(slope=qq_line_unique[[cn]][1],intercept=qq_line_unique[[cn]][2], colour="red"), linetype="dashed") +
      labs(x="Theoretical", y="Sample") + theme_minimal() +
      scale_color_manual(name=c("Minimum Cut", "Ratio Cut", "Normalized Cut")[i], values = c(colorvec[1:length(nvec)], "black", "red"),
                         labels=c(paste("n =", nvec), "Asymptotic distribution", "Assuming uniqueness"),
                         guide = guide_legend(override.aes = list(linetype = c(rep("blank", length(nvec)), "solid", "dashed"), shape = c(rep(16, length(nvec)), NA, NA)))) +
      theme(legend.position.inside = legend_pos, legend.box.background = element_rect(colour = "black"))
    
  }
  # For MCut, RCut and NCut, shows a quantile-quantile plot of the empirical distribution against the theoretical limit of the XCut statistic for different sample sizes from nvec.
  # The correct theoretical distribution is represented by the black line, and the distribution that assumes uniqueness of the XCut attaining partition is drawn in red.
  
  qq_plots$ccut_plot <- ggplot(empsample_df[empsample_df$cut=="ccut",]) + geom_point(aes(x=theoretical, y=sample, color=n), key_glyph="point") +
    geom_abline(mapping=aes(slope=1,intercept=0, colour="black")) +
    geom_abline(mapping=aes(slope=qq_line_unique$ccut[1],intercept=qq_line_unique$ccut[2], colour="red"), linetype="dashed") +
    geom_abline(mapping=aes(slope=qq_line_uneqvol$ccut[1],intercept=qq_line_uneqvol$ccut[2], colour="purple"), linetype="dashed") +
    labs(x="Theoretical", y="Sample") + theme_minimal() +
    scale_color_manual(name="Cheeger Cut", values = c(colorvec[1:length(nvec)], "black", "red", "purple"),
                       labels=c(paste("n =", nvec), "Asymptotic distribution", "Assuming unequal vol.", "Assuming uniqueness"),
                       guide = guide_legend(override.aes = list(linetype = c(rep("blank", length(nvec)), "solid", "dashed", "dashed"), shape = c(rep(16, length(nvec)), NA, NA, NA)))) +
    theme(legend.position.inside = legend_pos, legend.box.background = element_rect(colour = "black"))
  # Does the same for Cheeger Cut, with an additional line representing the assumption that any of the XCut attaining partitions is perfectly balanced volume-wise.
  
  return(c(empsample_ordered, qq_plots))
  # Returns the ordered sample in case that it will be used beyond the plotting.
  
}
# Generates quantile-quantile plots for the empirical against the theoretical limit distribution of the XCut statistics.
# Input: Number of sample points drawn for each distribution BigM, underlying graph g, vector of sample sizes n, colorvec containing possible point color values.
# Output: The sample used for the plotting. Additionally, this function yields the aforementioned QQ-Plots.

xcut_bootstrap_sample <- function(n, t, pmat, M_n, BigM, size){
  ks_arr <- unlist(lapply(1:size, function(i){
    g_smpl <- generate_grid_graph(matrix(rmultinom(1, n, pmat)/n, dim(pmat)[1]),t = t)
    return(rnorm_xcut_empirical(BigM, g_smpl, M_n, t)[-5])
  }), recursive = FALSE)
  return(split(unlist(ks_arr, use.names=FALSE), rep(names(ks_arr), lengths(ks_arr))))
}
# Generates an XCut bootstrap sample, i.e. from the distribution sqrt(M_n)*(XC(G_{M_n,n}) - XC(G_n)) for the graph G_{M_n,n}
# generated from a multinomial sample from the sample Y and the graph G_n generated from Y.
# Input: Sample size n, grid graph distance t, underlying probability vector pmat in matrix form (i.e. forming a grid),
#         bootstrap sample size M_n, number of bootstrap samples BigM, number of Y samples size
# Output: A list, indexed in the respective XCuts, containing BigM*size samples from the bootstrap distribution for each XCut.

rej_KS <- function(n, m, signif) sqrt(-log(signif/2)/2 * (n+m)/(n*m))
# Realizes the rejection threshold for a of a two-sided Kolmogorov-Smirnov test.
# Input: Sample sizes n and m of both samples, significance level alpha.
# Output: Rejection threshold for the Kolmogorov-Smirnov statistic D, i.e. reject if D > output.

ks_statistic_over_n_t <- function(nvec, tvec, pmat, M_nvec, BigM, size, size_theory, signif=0.05, progress_bar=TRUE){
  
  if(progress_bar) progress <- txtProgressBar(min = 1, max = length(nvec)*length(tvec), style = 3)
  
  xcut_ks_dfs <- rep(list(data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("value", "n", "t"))))), 4)
  # Initializes the dataframe in which the KS statistics will be stored for all cuts
  
  for(t_id in 1:length(tvec)){
    
    ks_theoretical_sample <- rnorm_xcut_limit(size_theory, generate_grid_graph(pmat, t=tvec[t_id]))
    # Generates the theoretical limit distribution.
    
    for(n_id in 1:length(nvec)){
      
      if(!is.nan(M_nvec)) ks_empirical_sample <- xcut_bootstrap_sample(nvec[n_id], tvec[t_id], pmat, M_nvec[n_id], BigM, size)
      else ks_empirical_sample <- rnorm_xcut_empirical(BigM, generate_grid_graph(pmat, t=tvec[t_id]), nvec[n_id], t=tvec[t_id])
      # Generates BigM bootstrap samples of size size each; these are simply concatenated in the resulting dataframe.
      
      for(cut in 1:4){
        xcut_ks_sample <- sapply(1:size, function(i)unname(ks.test(ks_empirical_sample[[cut]][(i-1)*BigM+(1:BigM)], ks_theoretical_sample[[cut]])$statistic))
        xcut_ks_dfs[[cut]][nrow(xcut_ks_dfs[[cut]])+(1:size),] <- cbind(xcut_ks_sample, rep(nvec[n_id],size), rep(tvec[t_id],size))
        # Computes size KS test statistics and stores them appropriately together with the respective n and t.
      }
      
      if(progress_bar) setTxtProgressBar(progress, n_id + (t_id - 1)*length(nvec))
      
    }
  }
  if(progress_bar) close(progress)
  
  names(xcut_ks_dfs) <- c("mcut_sample", "rcut_sample", "ncut_sample", "ccut_sample")
  return(xcut_ks_dfs)
}
# Generates instances of the Kolmogorov-Smirnov test of a XCut bootstrap (or empirical) sample against its limit distribution over n and t.
# Input: Vector of sample sizes nvec, vector of grid graph distances tvec, underlying probability matrix pmat, vector of bootstrap sample sizes M_nvec
#         (or, if the integer BigM denoting the size of the 
#         sample used for the KS test, integer size of number of KS tests to perform for each combination of n and t, integer size_theory
#         determining the size of the limit distribution sample used for the KS test which is performed at significance level signif.
# Output: A list, indexed in the respective XCuts, containing BigM*size samples from the bootstrap distribution sqrt(M_n)*(XC(G_{M_n,n}) - XC(G_n)),
#         where G_n is the graph generated from a multinomial sample of size nvec, and G_{M_n,n} is the graph generated from a bootstrap sample
#         of size M_nvec, and M_n is an entry of M_nvec

discretize_sample <- function(sample, xrange, yrange, xdiv, ydiv){
  
  smpl_disc <- matrix(rep(0,xdiv*ydiv),nrow=xdiv,ncol=ydiv)
  
  smpl_resc <- ceiling(cbind((sample[,1]-xrange[1])/xrange[2]*xdiv, (sample[,2]-yrange[1])/yrange[2]*ydiv))
  
  for(i in 1:dim(smpl_resc)[1]) smpl_disc[smpl_resc[i,1], smpl_resc[i,2]] <- smpl_disc[smpl_resc[i,1], smpl_resc[i,2]] + 1
  
  return(t(smpl_disc))
  
}
# Discretizes a (two-dimensional) sample onto a grid.
# Input: Two-column sample matrix sample, vector xrange of length 2 giving the range that the first component of sample will be discretized in,
#         analogous vector yrange for the second component of sample, numbers xdiv and ydiv of brackets to discretize sample into.
# Output: Matrix of size xdiv x ydiv containing the discretized sample, i.e. each matrix entry contains the number of sample points falling into that bin.



### PART II: APPLICATION EXAMPLES ###


# A toy example to showcase our discretization ----


pts_exmp <- matrix(c(0.6,2.2,0.7,5.2,1.5,3.3,2.7,2.6,2.2,3.4,3.2,4.3,2.7,1.2,3.6,3.1,1.1,4.5,1.9,4.8,0.5,1.0,
                     1.3,2.4,0.5,3.5,1.6,1.7,7.1,1.4,8.7,1.5,5.0,2.7,5.4,3.7,6.4,2.3,7.6,2.5,6.7,3.3,8.5,3.7,
                     6.6,4.4,7.5,4.7,4.3,4.4,8.4,2.5)/2, ncol=2, byrow=TRUE)
# Realizes the toy example from the paper.

g_exmp <- generate_sample_graph(pts_exmp, 0.55)
xc_exmp <- xcut_exact(g_exmp)
plot_graph_partition(g_exmp, xc_exmp$ncut, vertex_scale=10, edge_scale=5)
# Computes (using brute force) and then plots the NCut partition of the toy example above.

pts_disc <- discretize_sample(pts_exmp, c(0,5), c(0,3), 5, 3)
g_exmp_disc <- generate_grid_graph(pts_disc, t=1)
xc_exmp_disc <- xcut_exact(g_exmp_disc)
plot_graph_partition(g_exmp_disc, xc_exmp_disc$ncut, vertex_scale = 40, edge_scale = 25)
# Computes (again using brute force) and then plots the NCut partition of the discretized graph built from the toy example.


# Kolmogorov-Smirnov distance between the empirical XCut distribution and its theoretical limit ----


nvec <- seq(1000, 20000, 1000)
eps <- 0.4
prob_binom <- rep(c(sqrt(1+eps^2+1.5*eps), eps, 1), each=3) %>% (function(x)x/sum(x))
g_bimodal <- generate_grid_graph(matrix(prob_binom, 3), t=1)
g_coords <- cbind(V(g_bimodal)$x, V(g_bimodal)$y)
tvec <- sort(unique(c(outer(1:dim(g_coords)[1], 1:dim(g_coords)[1], Vectorize(function(i,j)sqrt(sum((g_coords[i,]-g_coords[j,])^2)))))))[-1]
# Sets up the KS test for the bimodal distribution over n (in nvec) and t (in tvec).

ks_ccut_df <- ks_statistic_over_n_t(nvec, tvec, matrix(prob_binom,3), M_nvec=NaN, BigM=100, size=1, size_theory=100)$ccut_sample
ks_ccut_ci_df <- cbind(ks_ccut_df, dks=rej_KS(n=ks_ccut_df$n, m=100, signif=0.975))
# Creates a dataframe containing KS statistics over values from nvec and tvec. Notice the option M_nvec=NaN and size=1 to compare empirical and limit distributions (and not bootstrap and limit).

ggplot(ks_ccut_ci_df, aes(x=n, y=value, colour=t)) + geom_line(linewidth=1.2) +
  geom_ribbon(aes(ymin=pmax(value-dks,0), ymax=pmin(value+dks,1), fill=t), alpha=0.2, linetype=0) +
  coord_cartesian(xlim=range(nvec), ylim=c(0,.85)) + labs(x="n", y="Kolmogorov-Smirnov statistic") + theme_minimal() +
  scale_color_manual(name=TeX("Cheeger Cut"), values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                     labels=unname(TeX(c("$1\\leq t \\,< \\,\\sqrt{2}$", "$\\sqrt{2}\\leq t \\,<\\;2\\,$", "$2\\leq t \\,< \\,\\sqrt{5}$",
                                         "$\\sqrt{5}\\leq t \\,< \\,\\sqrt{8}$", "$\\sqrt{8} \\leq t \\;\\;\\;\\;$"))),
                     guide = guide_legend(nrow=2, byrow=TRUE, override.aes = list(linetype = rep("solid", 5), shape = rep(NA, 5)))) +
  scale_fill_manual(name=TeX("Cheeger Cut"), values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                    labels=unname(TeX(c("$1\\leq t \\,< \\,\\sqrt{2}$", "$\\sqrt{2}\\leq t \\,<\\;2\\,$", "$2\\leq t \\,< \\,\\sqrt{5}$",
                                        "$\\sqrt{5}\\leq t \\,< \\,\\sqrt{8}$", "$\\sqrt{8} \\leq t \\;\\;\\;\\;$")))) +
  theme(legend.position = c(0.85,0.8), legend.box.background = element_rect(fill="white",colour = "black")) +
  scale_x_continuous(breaks=nvec) + guides(colour=guide_legend(ncol=1), fill=guide_legend(ncol=1))
# PLots the KS statistic of the CCut empirical distribution against its limit over values for n and t.

ks_ccut_df2 <- ks_statistic_over_n_t(nvec, tvec, matrix(rep(1/9,9),3), M_nvec=NaN, BigM=100, size=1, size_theory=100)$ccut_sample
ks_ccut_ci_df2 <- cbind(ks_ccut_df2, dks=rej_KS(n=ks_ccut_df2$n, m=100, signif=0.975))
# Creates a dataframe containing KS statistics over values from nvec and tvec for the uniform distribution as above. Can then be plotted as above.


# Quantile-quantile plots for a bimodal and the uniform distribution ----


eps <- 0.4
nvec <- c(50,200,500,1000,10000)

prob_binom <- rep(c(sqrt(1+eps^2+1.5*eps), eps, 1), each=3)
g_bimodal <- generate_grid_graph(matrix(prob_binom/sum(prob_binom), 3), t=1)
qq_smpls_binom <- generate_qq_plots_xcut_limit_distribution(2000, g_bimodal, nvec)
# Generates QQ plots of the XCut statistics against their respective limiting distributions for a bimodal distribution.

g_unif <- generate_grid_graph(matrix(rep(1/9,9), 3), t=1)
qq_smpls_unif <- generate_qq_plots_xcut_limit_distribution(2000, g_unif, nvec)
# Generates QQ plots as above for the uniform distribution.

qq_smpls_binom$ccut_plot
qq_smpls_unif$ccut_plot
# Displays and contrasts the CCut QQ plots of the bimodal and uniform distributions. Similarly, QQ plots for the other cuts can be displayed.


# Kolmogorov-Smirnov distance between the bootstrap XCut distribution and its theoretical limit ----


nvec <- round(10^seq(4,9,0.5))
eps <- 0.4
lambda <- sqrt(1 + 1.5*eps + eps^2)
pmat <- matrix(rep(c(1,eps,lambda),each=3),3) %>% (function(x)x/sum(x))
g_true <- generate_grid_graph(pmat)
g_true_coords <- cbind(V(g_true)$x, V(g_true)$y)
tvec <- sort(unique(c(outer(1:dim(g_true_coords)[1], 1:dim(g_true_coords)[1], Vectorize(function(i,j)sqrt(sum((g_true_coords[i,]-g_true_coords[j,])^2)))))))[-1]
# Initializes the necessary parameters, including all possible neighbourhood distances tvec.

ks_ccut_array <- ks_statistic_over_n_t(nvec, tvec, pmat, round(sqrt(nvec)), 100, 100, 2000)$ccut_sample
ks_ccut_array2 <- ks_statistic_over_n_t(nvec, tvec, pmat, nvec, 100, 100, 2000)$ccut_sample
# Generates instances of KS tests of the bootstrap distribution against the theoretical limit, first for M_n = sqrt(n), then for M_n = n.
#save(ks_ccut_array, file="ks_ccut_over_n_t_Mneqn.RData")

ks_ccut_means <- data.frame(value = mapply(function(n,t)median(ks_ccut_array$value[round(ks_ccut_array$n)==n & ks_ccut_array$t==t]), rep(nvec,each=length(tvec)), tvec),
                            n = rep(nvec,each=length(tvec)), t = factor(tvec), dks=rej_KS(n=100,m=2000,signif=0.975))
ks_ccut_means2 <- data.frame(value = mapply(function(n,t)median(ks_ccut_array2$value[round(ks_ccut_array2$n)==n & ks_ccut_array2$t==t]), rep(nvec,each=length(tvec)), tvec),
                             n = rep(nvec,each=length(tvec)), t = factor(tvec), dks=rej_KS(n=100,m=2000,signif=0.975))
# Transforms the output of the KS test function slightly to ready it for plotting, adding confidence band diameter dks.

ggplot(ks_ccut_means, aes(x=n, y=value, colour=t)) + geom_line(size=1.2) +
  geom_ribbon(aes(ymin=pmax(value-dks,0), ymax=pmin(value+dks,1), fill=t), alpha=0.2, linetype=0) +
  coord_cartesian(xlim=range(nvec), ylim=c(0,.85)) + labs(x="n", y="Kolmogorov-Smirnov statistic") + theme_minimal() +
  scale_color_manual(name=TeX("Cheeger Cut bootstrap for $M\\,=\\,\\sqrt{n}$"), values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                     labels=unname(TeX(c("$1\\leq t \\,< \\,\\sqrt{2}$", "$\\sqrt{2}\\leq t \\,<\\;2\\,$", "$2\\leq t \\,< \\,\\sqrt{5}$",
                                         "$\\sqrt{5}\\leq t \\,< \\,\\sqrt{8}$", "$\\sqrt{8} \\leq t \\;\\;\\;\\;$"))),
                     guide = guide_legend(nrow=2, byrow=TRUE, override.aes = list(linetype = rep("solid", 5), shape = rep(NA, 5)))) +
  scale_fill_manual(name=TeX("Cheeger Cut bootstrap for $M\\,=\\,\\sqrt{n}$"), values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                    labels=unname(TeX(c("$1\\leq t \\,< \\,\\sqrt{2}$", "$\\sqrt{2}\\leq t \\,<\\;2\\,$", "$2\\leq t \\,< \\,\\sqrt{5}$",
                                        "$\\sqrt{5}\\leq t \\,< \\,\\sqrt{8}$", "$\\sqrt{8} \\leq t \\;\\;\\;\\;$")))) +
  theme(legend.position = "none", legend.box.background = element_rect(fill="white",colour = "black")) +
  scale_x_continuous(trans="log10", breaks=nvec[c(1,3,5,7,9,11)], labels=trans_format("log10", math_format(10^.x)))

ggplot(ks_ccut_means2, aes(x=n, y=value, colour=t)) + geom_line(size=1.2) +
  geom_ribbon(aes(ymin=pmax(value-dks,0), ymax=pmin(value+dks,1), fill=t), alpha=0.2, linetype=0) +
  coord_cartesian(xlim=range(nvec), ylim=c(0,.85)) + labs(x="n", y="Kolmogorov-Smirnov statistic") + theme_minimal() +
  scale_color_manual(name=TeX("Cheeger Cut bootstrap for $M\\,=\\,n$"), values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                     labels=unname(TeX(c("$1\\leq t \\,< \\,\\sqrt{2}$", "$\\sqrt{2}\\leq t \\,<\\;2\\,$", "$2\\leq t \\,< \\,\\sqrt{5}$",
                                         "$\\sqrt{5}\\leq t \\,< \\,\\sqrt{8}$", "$\\sqrt{8} \\leq t \\;\\;\\;\\;$"))),
                     guide = guide_legend(nrow=2, byrow=TRUE, override.aes = list(linetype = rep("solid", 5), shape = rep(NA, 5)))) +
  scale_fill_manual(name=TeX("Cheeger Cut bootstrap for $M\\,=\\,n$"), values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                    labels=unname(TeX(c("$1\\leq t \\,< \\,\\sqrt{2}$", "$\\sqrt{2}\\leq t \\,<\\;2\\,$", "$2\\leq t \\,< \\,\\sqrt{5}$",
                                        "$\\sqrt{5}\\leq t \\,< \\,\\sqrt{8}$", "$\\sqrt{8} \\leq t \\;\\;\\;\\;$")))) +
  theme(legend.position = c(0.65,0.85), legend.box.background = element_rect(fill="white",colour = "black"), legend.title=element_blank()) +
  scale_x_continuous(trans="log10", breaks=nvec[c(1,3,5,7,9,11)], labels=trans_format("log10", math_format(10^.x)))
# Plots the KS test value (with confidence band) against the vector nvec of sample sizes n, first for bootstrap sample size M_n = sqrt(n), then M_n = n.


# Rejection probability of an asymptotic test for XCuts ----


alph <- 0.05
nvec <- c(1000,5000,10000)
epsvec <- seq(0.1,2.5,0.1)

test_rejection_over_n_eps <- function(g_lst=c(), g_unif=generate_grid_graph(matrix(rep(1/16,16),4)),
                                      nvec, epsvec, BigN, BigN_limit, alph=0.05, progress_bar=TRUE){
  
  if(progress_bar) progress <- txtProgressBar(min = 1, max = length(nvec)*length(epsvec), style = 3)
  
  tr_df <- c()
  
  limit_sample <- apply(rnorm_xcut_limit(BigN_limit, g_unif, partitions=xcut_exact(g_unif)$ncut)$ncut_sample, 1, min)
  limit_quantile <- quantile(limit_sample, probs=c(alph/2, 1-alph/2))
  # Generates a sample from the limit distribution and uses it to compute the quantiles for the two-sided test.
  
  for(e_id in 1:length(epsvec)){
    for(n_id in 1:length(nvec)){
      
      if(length(pmat_lst) > 0) g_eps <- g_lst[[n_id + (e_id-1)*length(nvec)]]
      else{
        pvec <- c(rep(1,4), rep(epsvec[e_id],8), rep(1,4))
        g_eps <- generate_grid_graph(matrix(pvec/sum(pvec),4))
      }
      # Creates the graph to compare the uniform distribution one to; takes pmat_lst if it is supplied.
      
      empirical_sample <- c(rnorm_xcut_empirical(BigN, g_eps, nvec[n_id], g_truth=g_unif, progress_bar=FALSE, degree_up_to_date=FALSE)$ncut_sample)
      # Generates BigN samples from the XCut statistic where the true underlying graph is assumed to be generated from the uniform distribution.
      
      tr_df <- rbind(tr_df, data.frame(prob=sum(empirical_sample < limit_quantile[1] | empirical_sample > limit_quantile[2], na.rm=TRUE)/BigN, eps=epsvec[e_id], n=nvec[n_id]))
      # Calculates the empirical rejection probability and collects it in tr_df
      
      if(progress_bar) setTxtProgressBar(progress, n_id + (e_id - 1)*length(nvec))
      
    }
  }
  
  close(progress)
  
  return(tr_df)
  
}
# Empirical rejection probability of a two-sided asymptotic test of the XCut value of a given distribution against the uniform.
# Input: A list g_lst (of length length(nvec)*length(epsvec)) of graphs whose (n_id+(eps_id-1)*length(nvec))-th entry has been generated using nvec[n_id] and epsvec[eps_id]
#         (if g_lst is empty (the default), then a simple 4x4 bi-/unimodal distribution as in the paper is used), graph g_unif to test against (default: 4x4 uniform distribution),
#         vector nvec of samples sizes (to iterate over), vector epsvec of probability distribution parameters (to iterate over), number of XCut statistic computation iterations
#         to use for rejection probability estimation BigN, size BigN_limit of the sample used to estimate the limit quantile, significance level alph,
#         boolean progress_bar indicating of whether to draw a progress bar.
# Output: A data frame listing the empirical rejection probability with the corresponding values n and eps from nvec and epsvec, respectively, that were used to set up the test.

test_rej <- test_rejection_over_n_eps(nvec=nvec, epsvec=epsvec, BigN=5, BigN_limit=100, alph=alph)

colvec <- c("green4", "lawngreen", "darkgoldenrod1", "blue", "orangered1", "orange1", "yellow2", "mediumpurple2")
# Defines the colors for plotting the different values in nvec.

ggplot(test_rej, aes(x=eps, y=prob, colour=factor(n))) + geom_line(linewidth=1.2) + geom_abline(slope=0, intercept=alph, colour="red") +
  scale_x_continuous(breaks=seq(0,2.5,0.5)) + scale_y_continuous(breaks=seq(0,1,0.2)) + theme_minimal() +
  labs(x=unname(TeX("$\\epsilon$")) , y="Rejection probability") +
  scale_color_manual(name="Sample size", values=colvec[1:length(nvec)], labels=paste("n =",nvec)) +
  theme(legend.position.inside = c(0.2,0.81), legend.box="horizontal",
        legend.box.background = element_rect(color = "black",fill="white"), legend.text=element_text(size=14),
        text = element_text(size=rel(4)), axis.text=element_text(size=14), panel.background = element_rect(fill="gray98", color="white")) +
  guides(fill=guide_legend(order=1), color=guide_legend(order=2))
# Plots the rejection probability against epsvec for different sample sizes from nvec.

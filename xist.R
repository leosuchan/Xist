############################################################################################################################
# R code supplementary to the paper "A scalable clustering algorithm to approximate graph cuts" by Suchan, Li, Munk (2023) #
############################################################################################################################
# Part I: Function definitions


install.packages(c("igraph", "ggplot2", "scales", "raster", "rasterVis", "latticeExtra", "RANN", "dbscan", "sClust", "RColorBrewer"))
# Required packages:
# igraph: Required for essentially everything; igraph supplies the graph structure and many useful functions.
# ggplot: Required for plotting only.
# raster: Required for loading and cropping the confocal raster images.
# rasterVis and latticeExtra: Required for plotting raster images with overlaid colors.
# RANN: Used for its implementation of the k-nearest neighbor algorithm.
# dbscan: Used for its implementation of the DBSCAN algorithm.
# sClust: Used for its implementation of the standard spectral clustering algorithm as stated e.g. in von Luxburg (2007).
# RColorBrewer: Used for vanity only.

library(igraph)
# Implements the graph structure and provides an algorithm to compute st-MinCuts.
library(ggplot2)
# Used for most of the plotting.
library(scales)
# Used for adjusting the plot scaling.
library(raster)
# Enables reading of .tif data such as the cell image dataset used in the paper..

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
  # Computes the degree (i.e. the sum of the weights of all neighboring vertices) for each vertex.
  
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
  
  if(is.null(E(g)$weight)) E(g)$weight <- 1
  
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
  
  if(is.null(E(g)$weight)) E(g)$weight <- 1
  
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
  
  if(is.null(E(g)$weight)) E(g)$weight <- 1
  
  if(length(S) == length(V(g))) S <- which(S)
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
  
  if(is.null(E(g)$weight)) E(g)$weight <- 1
  
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
  
  if(is.null(E(g)$weight)) E(g)$weight <- 1
  
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

plot_partition_overlay <- function(sample_raster, partition, colorvec){
  
  k <- length(unique(c(partition)))
  pixels <- dim(sample_raster)[1:2]
  M <- dim(partition)
  if(length(M) != 2) simpleError("partition must be a 2d matrix")
  # Determines number of partitions k, image dimension pixels and grid dimension M.
  
  rasterVis::levelplot(flip(t(sample_raster), 2), margin = FALSE, colorkey=FALSE, par.settings = list(axis.line=list(col='transparent')),
                       scales = list(draw=FALSE), col.regions = colorRampPalette(RColorBrewer::brewer.pal(9, 'Greys'))) +
    # Plots a greyscale image of sample_raster to add the colors onto.
    latticeExtra::layer(lapply(1:k, function(i){
      
      partitionvec <- which(partition==i)
      # For each color, this selects all vertex indices that will be colored in one color.
      
      lwidth <- 0.5 * pixels / M
      # Computes half of the width (in each dimension) of one square.
      
      locations <- rbind((partitionvec-1) %/% M[2] + 0.5, (partitionvec-1) %% M[2] + 0.5) * pixels / M
      # Computes the locations of the centers of the rectangles corresponding to the indices in partitionvec,
      # given that it is based on a pixels[1] x pixels[2] grid that is discretized onto an M[1] x M[2] sized grid.
      
      locs_quad_offset <- cbind(-lwidth, c(-1,1)*lwidth, lwidth, c(1,-1)*lwidth, -lwidth)
      loc_poly <- lapply(1:length(partitionvec), function(j)sp::Polygon(t(locations[,j] - locs_quad_offset), hole=FALSE))
      names(loc_poly) <- paste0("square", partitionvec)
      buffer_locations <- sp::SpatialPolygons(list(polys = sp::Polygons(loc_poly, ID=i)))
      # Constructs a list of squares (through sp::Polygon) of appropriate size at each location of partitionvec (given by locations).
      
      return(sp::sp.polygons(buffer_locations, fill=colorvec[i], col=NA, alpha=0.5))
      # Returns the appropriately colored polygons with slight transparency.
      
    }), data=list(k=k,partition=partition,M=M,pixels=pixels,colorvec=colorvec))
  # Plots a layer of k colors via the "latticeExtra" layer function whose handling is slightly awkward.
  
}
# Plots the raster image sample_raster together with a colored partition overlay.
# Input: Raster image sample_raster, matrix of discretized partition assignments partition, vector of colors colorvec
# Output: None; plots the graph.

manual_multiway_spec_clustering <- function(g, k){
  
  if(is.null(E(g)$weight)) E(g)$weight <- 1
  
  L_mat_norm <- diag(rep(0,length(V(g))))
  edgemat <- ends(g, E(g))
  L_mat_norm[rbind(edgemat, cbind(edgemat[,2],edgemat[,1]))] <- - E(g)$weight
  L_mat_norm <- - L_mat_norm / rowSums(L_mat_norm)
  L_mat_norm[cbind(V(g), V(g))] <- 1
  # Computes the normalized graph Laplacian L_mat_norm = I - D^-1 W for weight matrix W and degree matrix W.
  
  L_mat_norm_eigen <- eigen(L_mat_norm)
  # Eigenvalues and eigenvectors of the Laplacians.
  
  U_mat_norm <- L_mat_norm_eigen$vectors[,(ncol(L_mat_norm_eigen$vectors)-k+1):ncol(L_mat_norm_eigen$vectors)]
  # Selects the matrix consisting of the k eigenvectors corresponding to the k smallest eigenvalues of L_mat_norm.
  
  set.seed(8472)
  # Fixes seed to allow for replicability since kmeans uses the non-deterministic Hartigan-Wong (1979) algorithm.
  
  return(lapply(2:k, function(i)kmeans(U_mat_norm[,(k-i+1):k], centers=i)$cluster))
  # Applies the kmeans algorithm to partition the entries of U_mat_norm (corresponding to the observations) into 2,...,k clusters.
  
}
# Clusters the graph g using normalized spectral clustering according to Shi and Malik (2000); algorithm due to von Luxburg (2007).
# Input: igraph graph g, number of clusters k.
# Output: A vector of integers 1 or 2 representing the cluster the nodes are allocated to.

classification_rate <- function(vec, truth){
  
  if(length(vec) != length(truth)) simpleError("vec and truth differ in length")
  # vec and truth have to be of the same length to allow comparison.
  
  vec_unique <- unique(vec)
  vec_intersect <- sapply(vec_unique, function(i)table(factor(truth[vec == i], levels=1:2)))
  # For each distinct entry of vec, checks how many classified entries are in truth either 1 or 2 (the only levels of truth).
  
  m1_idx <- which.max(vec_intersect[1,])
  m2_idx <- which.max(vec_intersect[2,])
  # Selects both entries that most often are classified as 1 and 2.
  
  if(m1_idx == m2_idx){
    # If both indices are equal, their corresponding sum should not be returned since vec doesn't distinguish 1 and 2 in this case.
    
    m1_2 <- ifelse(length(vec_unique) > 1, max(vec_intersect[1,-m1_idx]), 0)
    m2_2 <- ifelse(length(vec_unique) > 1, max(vec_intersect[2,-m2_idx]), 0)
    # Additionally selects the second best values (guaranteed to be different from m1_idx and m2_idx).
    # Note that this also handles the fringe case of vec only containing one distinct label.
    
    return(unname(max(vec_intersect[1,m1_idx] + m2_2, vec_intersect[2,m2_idx] + m1_2)) / length(vec))
    # Returns the best match of either m1_idx and m2_2 or m2_idx and m1_2.
    
  } else return(unname(vec_intersect[1,m1_idx] + vec_intersect[2,m2_idx]) / length(vec))
  # Simply returns the best possible sum of the number of entries classified into two categories.
  
}
# Computes the classification rate of vec w.r.t. the true classification truth.
# Input: Two vectors vec and truth, with truth only containing 1's or 2's.
# Output: Classification rate, i.e. the fraction of entries correctly classified under the best possible permutation of entries of vec.

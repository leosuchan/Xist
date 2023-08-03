####################################################################################################
# R code supplementary to the Master's thesis "Limits of graph cuts on finite grids" by Leo Suchan #
####################################################################################################

install.packages(c("igraph", "ggplot2", "raster", "mvtnorm", "tidyr", "RColorBrewer", "rasterVis"))
# Required packages. mvtnorm is only needed once to draw samples from a multivariate normal distribution.
# From tidyr, only the pivot_longer function is used for some data formatting. RColorBrewer is used for vanity only.
# The levelplot function from rasterVis is needed to plot NIH3T3 samples.

### PART I: FUNCTION DEFINITIONS ###

library(igraph)
# Implements the graph structure and provides an algorithm to compute st-MinCuts.
library(ggplot2)
# Used for most of the plotting.
library(raster)
# Enables reading of .tif data such as the NIH3T3 dataset used in the thesis.

# NOTE: Descriptions of the following self-implemented functions are given immediately *after* the function.
# The inner workings of the functions are sparsely commented to explain the way these functions operate.
# It is recommended to use an IDE that allows folding of {...} environments to be able to fold all function
# definitions for clarity and easier future referencing.

generate_grid_graph <- function(sample, t=1, sigm=NA, weights_by_sample=TRUE){
  
  if(t < 1) simpleError("t must be greater or equal than 1")
  if(length(dim(sample)) != 2) simpleError("sample must be 2-dimensional")
  ndim <- dim(sample)
  n <- ndim[1] * ndim[2]
  
  edgevec <- weightvec <- x1coordvec <- x2coordvec <- labelvec <- c()
  degreevec <- rep(0, n)
  
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
              degreevec[x1 + (x2-1)*ndim[2]] <- degreevec[x1 + (x2-1)*ndim[2]] + current_weight
              degreevec[y1 + (y2-1)*ndim[2]] <- degreevec[y1 + (y2-1)*ndim[2]] + current_weight
              
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
  
  return(g)
  
}
# TODO description


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
  
  
  for(x in 1:n){
    
    # <- (x+1):n
    
    if(!is.na(nearest_neighbors) && nearest_neighbors >= 1){
    
      for(i in 2:dim(knn_sample$nn.idx)[2]){
        # Ignore the first entry of knn_sample since this is always the point i itself
        
        i_idx <- knn_sample$nn.idx[x,i]
        
        if(sum((sample[x,] - sample[i_idx,])^2) > t){
          # Only additionally connect nearest neighbors if they would remain unconnected otherwise.
        
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
        
        if(sum((sample[x,] - sample[y,])^2) <= t){
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
  #print(edge_mat)
  
#  edge_dupes <- duplicated(edge_mat, nmax=2)
#  edge_mat <- edge_mat[-edge_dupes,]
#  weightvec <- weightvec[-edge_dupes]
#  # Detects duplicate edges (from nearest neighbors and t-radius).
  
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
  
  V(g)$degree <- degreevec #sapply(1:n, function(i)sum(weightvec[edge_mat[,1] == i || edge_mat[,2] == i]))
  # Computes the degree (i.e. the sum of the weights of all neighbouring vertices) for each vertex.
  
  return(g)
  # Returns the graph g
}
# TODO description

xcut_via_stmincut <- function(g, progress_bar=TRUE){
  
  ij_min_m <-ij_min_r <- ij_min_n <- ij_min_c <- Inf
  ij_min_mcut <- ij_min_rcut <- ij_min_ncut <- ij_min_ccut <- c()
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = length(V(g))*(length(V(g))+1)/2, style = 3)
  
  for(i in 1:(length(V(g))-1)){
    
    for(j in (i+1):length(V(g))){
      
      ij_mincut <- min_cut(g, source = i, target = j, capacity = E(g)$weight, value.only = FALSE)$partition1
      # Computes a partitions attaining the ij-MinCut.
      
      S_ij <- ij_mincut #rep(1,length(V(g))) %>% (function(a){a[ij_mincut] <- 2; a})
      
      ij_xcut <- xcut_partition_S(S_ij, g)
      # Computes the XCut of the partition S_ij.
      
      if(!is.na(ij_xcut$mcut_min) & ij_xcut$mcut_min < ij_min_m){
        ij_min_m <- ij_xcut$mcut_min
        ij_min_mcut <- S_ij
      }
      
      if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
        ij_min_n <- ij_xcut$ncut_min
        ij_min_ncut <- S_ij
      }
      
      if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
        ij_min_c <- ij_xcut$ccut_min
        ij_min_ccut <- S_ij
      }
      
      if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
        ij_min_r <- ij_xcut$rcut_min
        ij_min_rcut <- S_ij
      }
      # Updates the best XCut and the partition that attains it
      
      if(progress_bar) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
      # Updates progress bar.
    }
  }
  
  if(progress_bar) close(progress)
  
  return(list(mcut = ij_min_mcut, mcut_min = ij_min_m,
              rcut = ij_min_rcut, rcut_min = ij_min_r,
              ncut = ij_min_ncut, ncut_min = ij_min_n,
              ccut = ij_min_ccut, ccut_min = ij_min_c))#, bla=bla))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}
# TODO description

xvst_on_local_maxima <- function(g, locmax_by_degree=FALSE, progress_bar=TRUE){
  
  mcvec <- c()
  
  ij_min_m <- ij_min_r <- ij_min_n <- ij_min_c <- Inf
  ij_min_mcut <- ij_min_rcut <- ij_min_ncut <- ij_min_ccut <- c()
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  n <- length(V(g))
  
  if(locmax_by_degree) locmaxvec <- which(sapply(1:n, function(i)all(sapply(V(g)[.nei(i)], function(j)V(g)$degree[i] >= V(g)$degree[j]))))
  else locmaxvec <- which(sapply(1:n, function(i)all(sapply(V(g)[.nei(i)], function(j)V(g)$sample[i] >= V(g)$sample[j]))))
  # Computes the set of local maxima of g.
  print(locmaxvec)
  #locmaxvec <- 1:n
  
  N <- length(locmaxvec)
  
  if(N == 1) return(list(mcut = c(), mcut_min = Inf, rcut = c(), rcut_min = Inf, ncut = c(), ncut_min = Inf, ccut = c(), ccut_min = Inf))
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = N*(N+1)/2, style = 3)
  # Setting up the progress bar.
  
  for(i_idx in 1:(N-1)){
    
    for(j_idx in (i_idx+1):N){
      
      ij_mc <- min_cut(g, source = locmaxvec[i_idx], target = locmaxvec[j_idx], capacity = E(g)$weight, value.only = FALSE)
      ij_mincut <- ij_mc$partition1
      mcvec <- c(mcvec, ij_mc$value)
      # Computes a partitions attaining the ij-MinCut.
      
      S_ij <- ij_mincut#rep(1,n) %>% (function(a){a[ij_mincut] <- 2; a})
      
      ij_xcut <- xcut_partition_S(S_ij, g)
      # Computes the XCut of the partition S_ij.
      
      if(!is.na(ij_xcut$mcut_min) & ij_xcut$mcut_min < ij_min_m){
        ij_min_m <- ij_xcut$mcut_min
        ij_min_mcut <- S_ij
      }
      
      if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
        ij_min_n <- ij_xcut$ncut_min
        ij_min_ncut <- S_ij
      }
      
      if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
        ij_min_c <- ij_xcut$ccut_min
        ij_min_ccut <- S_ij
      }
      
      if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
        ij_min_r <- ij_xcut$rcut_min
        ij_min_rcut <- S_ij
      }
      # Updates the best XCut and the partition that attains it.
      
      if(progress_bar) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
      # Updates progress bar.
      
    }
  }
  
  if(progress_bar) close(progress)
  
  return(list(mcut = ij_min_mcut, mcut_min = ij_min_m,
              rcut = ij_min_rcut, rcut_min = ij_min_r,
              ncut = ij_min_ncut, ncut_min = ij_min_n,
              ccut = ij_min_ccut, ccut_min = ij_min_c))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}
# TODO description

xcut_exact <- function(g, normalize=FALSE, progress_bar=TRUE){
  
  mcut_min <- rcut_min <- ncut_min <- ccut_min <- Inf
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
  
  if(n_eff > 30){
    print("Error: More than 30 nodes.")
    return(NULL)
  }
  # Hard-coded check whether the number of nodes isn't too large in order to keep computations manageable.
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = 2^(n_eff-1)-1, style = 3)
  
  for(k in 1:(2^(n_eff-1)-1)){
    
    p <- vvec[intToBits(k)[1:n_eff]==1]
    # Converts a number k into a partition p by binary encoding.
    
    xcut_p <- xcut_partition_S(p, g)
    # Computes the XCut of partition p.
    #    print(k)
    #    print(p)
    #    print(xcut_p)
    
    if(xcut_p$mcut_min < mcut_min){
      mcut_min <- xcut_p$mcut_min
      mcut <- p
    }
    if(xcut_p$rcut_min < rcut_min){
      rcut_min <- xcut_p$rcut_min
      rcut <- p
    }
    if(xcut_p$ncut_min < ncut_min){
      ncut_min <- xcut_p$ncut_min
      ncut <- p
    }
    if(xcut_p$ccut_min < ccut_min){
      ccut_min <- xcut_p$ccut_min
      ccut <- p
    }
    # Updates the optimal XCuts and the partition attaining them.
    
    if(progress_bar) setTxtProgressBar(progress, k)
    
  }
  return(list(mcut=mcut, rcut=rcut, ncut=ncut, ccut=ccut,
              mcut_min=mcut_min, rcut_min=rcut_min, ncut_min=ncut_min, ccut_min=ccut_min))
}
# TODO description

xcut_partition_S <- function(S, g, normalize=FALSE){
  
  S_C <- V(g)[-S]
  
  if(F){ #!any(is.null(V(g)$degree))){
    
    S_vol <- sum(V(g)[S]$degree)
    S_C_vol <- sum(V(g)[S_C]$degree)
    
  } else {
    
    S_vol <- sum(E(g)[.inc(S)]$weight) + sum(E(g)[S %--% S]$weight)
    S_C_vol <- sum(E(g)[.inc(S_C)]$weight) + sum(E(g)[S_C %--% S_C]$weight)
    
  }
  # Computing the volumes of S and S_C (also computable using sum(E(g)[.inc(S)]$weight) if vertex degree attribute is missing)
  
  if(normalize) edgesum <- sum(E(g)$weight)
  
  mcut_min <- sum(E(g)[S %--% S_C]$weight) / ifelse(normalize, edgesum, 1)
  rcut_min <- mcut_min / (length(S) * length(S_C)) * ifelse(normalize, length(V(g)), 1)
  ncut_min <- mcut_min / (S_vol * S_C_vol) * ifelse(normalize, (S_vol + S_C_vol)^2, 1)
  ccut_min <- mcut_min / min(S_vol, S_C_vol) * ifelse(normalize, (S_vol + S_C_vol), 1)
  # Computes the various XCuts of S.
  
  return(list(mcut_min=mcut_min, rcut_min=rcut_min, ncut_min=ncut_min, ccut_min=ccut_min))
  # Returns a list of the computed XCut values of S.
}
# TODO description

xist <- function(g, allcuts=TRUE, progress_bar=TRUE){
  
  ij_min_r <- ij_min_n <- ij_min_c <- Inf
  ij_min_rcut <- ij_min_ncut <- ij_min_ccut <- c()
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  if(locmax_by_degree) locmaxvec <- which(sapply(1:length(V(g$graph)),
                                                 function(i)all(sapply(V(g)[.nei(i)], function(j)V(g$graph)$degree[i] >= V(g$graph)$degree[j]))))
  else locmaxvec <- which(sapply(1:length(V(g$graph)),
                                 function(i)all(sapply(V(g)[.nei(i)], function(j)V(g$graph)$sample[i] >= V(g$graph)$sample[j]))))
  locmaxvec <- sort(locmaxvec, decreasing=TRUE)
  # Computes and sorts the set of local maxima of g.
  
  g_merge <- g
  # Denote the graph where all local maxima will be merged by g_merge
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = length(locmaxvec)-1, style = 3)
  # Setting up the progress bar.
  
  while(length(locmaxvec) > 1){
    #    print(paste(locmaxvec[1],"|",locmaxvec[i]))
    #     cat("distmat:",distmat[1,],"\n")
    
    ij_cut <- min_cut(g, source = locmaxvec[i_idx], target = locmaxvec[j_idx], capacity = E(g_merge)$weight, value.only = FALSE)$partition1
    # Computes all partitions attaining the ij-MinCut.
    
    if(locmaxvec[1] %in% S_ij) S_ij  <- (V(g)[-S_ij]) #S_ij <- union(S_ij, already_merged)#c(S_ij, locmaxvec[2:(i-1)])
    else S_ij <- (V(g_directed)[S_ij])$tempnr#S_ij <- setdiff(S_ij, already_merged)
    #     print(S_ij)
    
    ij_xcut <- xcut_partition_S(S_ij, g)
    # Computes the XCut of the partition S_ij.
    
    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
      ij_min_n <- ij_xcut$ncut_min
      ij_min_ncut <- S_ij
      #        print(paste(locmaxvec[1],"|",locmaxvec[i]))
      #        print(S_ij)
    }
    if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
      ij_min_c <- ij_xcut$ccut_min
      ij_min_ccut <- S_ij
    }
    if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
      ij_min_r <- ij_xcut$rcut_min
      ij_min_rcut <- S_ij
    }
    # Updates the best XCut and the partition that attains it.
  
    
    if(length(locmaxvec) > 2){
      # We now merge vertices locmaxvec[1] and locmaxvec[i] by combining them into one (namely vertex locmaxvec[1]).
      
      mergelist <- (1:length(V(g_directed))) %>% (function(a){a[locmaxvec[2]] <- locmaxvec[1];a})
      #       cat("mergelist:",mergelist,"\n")
      g_directed <- contract(g_directed, mergelist, vertex.attr.comb=list(degree="sum", sample="sum", "first"))
      # Contracts all vertices in vertices_to_merge into one (namely locmaxvec[1]).
      
      g_directed <- simplify(g_directed, edge.attr.comb=list(capacity="sum", weight="sum", width="sum", distance="min", "first"))
      # Collapses all resulting multiple edges into one by summing their weights etc., and deletes loops, thus finishing the merging.
      
      if(locmaxvec[2] == length(V(g_directed))) #g_directed <- delete.vertices(g_directed, locmaxvec[mindist_idx])
      
      #for(i in 2:length(locmaxvec))print(distmat[c(1,mindist_idx),i])
      #    print(distmat)
      
      distmat[1,2:length(locmaxvec)] <- sapply(2:length(locmaxvec), function(i)min(distmat[c(1,mindist_idx),i]))
      distmat <- distmat[-mindist_idx,-mindist_idx]
      # Updates distmat to account for vertices_to_merge being merged into locmaxvec[1] (=vertices_to_merge[1])
      
      #    already_merged <- c(already_merged, locmax[mindist_idx])
      # Updates the already merged vertices.
    }
    locmaxvec <- locmaxvec[-mindist_idx] %>% (function(a){a - (a > mindist_idx)})
    #distmat <- distmat[,-mindist_idx]
    # Prevents locmaxvec[mindist_idx] from being merged again.
    
    if(progress_bar) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
    # Updates progress bar.
    
    #print(i)
    #print(j)
    #print(S_ij)
    #print(ij_xcut$ncut_min)
    
  }
  if(progress_bar) close(progress)
  #    bla <- c()
  #    for(i in 1:(length(V(g$graph))-1)){
  #      for(j in (i+1):length(V(g$graph))){
  #        ij_cut <- stMincuts(g_directed, i, j)
  #        # Computes all partitions attaining the ij-MinCut.
  #        
  # #       for(S_ij in ij_cut$partition1s){
  #    ij_xcut <- xcut_partition_S(S_ij, g)
  #    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min == ij_min_n){
  #      bla <- c(bla, i, j)
  #  }#}}
  #}
  return(list(rcut = ij_min_rcut, rcut_min = ij_min_r,
              ncut = ij_min_ncut, ncut_min = ij_min_n,
              ccut = ij_min_ccut, ccut_min = ij_min_c))#, bla=bla))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}

plot_graph_partition <- function(g, partition=NULL, vertex_scale=400/sqrt(length(V(g))), axes=FALSE, edge_scale=20,
                                 vertex_offset=8/log(length(V(g))), edge_offset=1, main=NULL,
                                 xlim=range(V(g)$x), ylim=range(V(g)$y)){
  #V(g$graph)[unlist(partition)]$color <- 'blue'
  #V(g$graph)[-unlist(partition)]$color <- 'green3'
  if(!is.null(partition)){
    if(length(partition) < length(V(g))){
      V(g)$color <- 'green3'
      #V(g$graph)[which(as.numeric(V(g$graph)$label) %in% partition)]$color <- 'blue'
      V(g)[unlist(partition)]$color <- 'blue'
    } else if(length(partition) == length(V(g))){
      if(is.numeric(partition)){
        allcolors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
        V(g)$color <- c("green3", "blue", sample(allcolors, length(unique(partition))-2))[partition]
        #colorRampPalette(c("red", "green3"))
      }
      else V(g)$color <- partition
    }
  }
  edge_weights <- E(g)$weight
  vertex_weights <- V(g)$sample
  #vertex_weights <- V(g$graph)$smmm
  
  par(cex=0.7, mai=c(0,0.1,0,0.1))
  #plot(g$graph, edge.width = E(g$graph)$weight, vertex.size=vertex_size, axes=TRUE)
  if(axes){
    plot(g, rescale=FALSE, axes=TRUE, xlim=xlim, ylim=ylim, asp = 0, vertex.label=NA,
         vertex.size = vertex_size, main=main)
  }else {
    plot(g, vertex.label=NA, edge.width = edge_scale*edge_weights/max(edge_weights)+edge_offset,
         #vertex.size = vertex_scale*vertex_weights/1.5+15, main=main)
         vertex.size = vertex_scale*vertex_weights/max(vertex_weights)+vertex_offset, main=main)
  }
}
# TODO description


xist_gusfield <- function(g, locmax_by_degree=FALSE, normalize=FALSE, progress_bar=TRUE){
  
  n <- length(V(g))
  
  ij_min_m <- ij_min_r <- ij_min_n <- ij_min_c <- Inf
  ij_min_mcut <- ij_min_rcut <- ij_min_ncut <- ij_min_ccut <- c()
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  if(locmax_by_degree){
    
    if(is.null(V(g)$degree)) V(g)$degree <- sapply(1:n, function(i)sum(E(g)[.inc(i)]$weight))
    
#    locmaxvec <- rep(TRUE,n)
    #edgemat <- ends(g, 1:length(E(g)))
#    for(i in 1:length(E(g))){
#      
##      oneedge <- ends(g,i)
#      if(V(g)[oneedge[1]]$degree > V(g)[oneedge[2]]$degree) locmaxvec[oneedge[2]] <- FALSE
#      else if(V(g)[oneedge[2]]$degree > V(g)[oneedge[1]]$degree) locmaxvec[oneedge[1]] <- FALSE
#      
#    }
    
#    locmaxvec <- which(locmaxvec)
    
    locmaxvec <- which(sapply(1:n, function(i)all(sapply(V(g)[.nei(i)], function(j)V(g)$degree[i] >= V(g)$degree[j]))))
  
  } else{
    
    if(is.null(V(g)$sample)) V(g)$sample <- 1
    
    locmaxvec <- which(sapply(1:n, function(i)all(sapply(V(g)[.nei(i)], function(j)V(g)$sample[i] >= V(g)$sample[j]))))
    
  }
  # Computes the set of local maxima of g.
  
  #locmaxvec <- 1:n
  N <- length(locmaxvec)
  #print(N)
  eqvec <- rep(1, N)
  #mcvec <- rep(NA,N-1)
  
  if(N == 1) return(list(mcut = c(), mcut_min = Inf, rcut = c(), rcut_min = Inf, ncut = c(), ncut_min = Inf, ccut = c(), ccut_min = Inf))
  
  if(progress_bar) progress <- txtProgressBar(min = 1, max = N-1, style = 3)
  
  for(s in 2:N){
    
    t <- eqvec[s]
    
    ij_cut <- min_cut(g, source = locmaxvec[s], target = locmaxvec[t], capacity = E(g)$weight, value.only = FALSE)
    #ij_cut <- min_cut(g, source = s, target = t, capacity = E(g)$weight, value.only = FALSE)
    S_ij <- ij_cut$partition1
    #mcvec[s-1] <- ij_cut$value
    
    if(!(locmaxvec[s] %in% S_ij)) S_ij  <- V(g)[-S_ij]
    #if(!(s %in% S_ij)) S_ij  <- V(g)[-S_ij]
    
    #print(paste(s,"|",t))
    #print(length(S_ij))
    
    for(i in s:N) if(locmaxvec[i] %in% S_ij && eqvec[i] == t) eqvec[i] <- s
    #for(i in s:n){ if(i %in% S_ij && eqvec[i] == t) eqvec[i] <- s }
#    if(s < 6){
#      print(S_ij)
      
#    print(head(eqvec))}
    
    ij_xcut <- xcut_partition_S(S_ij, g, normalize=normalize)
    # Computes the XCut of the partition S_ij.
    
    if(!is.na(ij_xcut$mcut_min) & ij_xcut$mcut_min < ij_min_m){
      ij_min_m <- ij_xcut$mcut_min
      ij_min_mcut <- S_ij
    }
    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
      ij_min_n <- ij_xcut$ncut_min
      ij_min_ncut <- S_ij
      #        print(paste(locmaxvec[1],"|",locmaxvec[i]))
      #        print(S_ij)
    }
    if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
      ij_min_c <- ij_xcut$ccut_min
      ij_min_ccut <- S_ij
    }
    if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
      ij_min_r <- ij_xcut$rcut_min
      ij_min_rcut <- S_ij
    }
    
    if(progress_bar) setTxtProgressBar(progress, s)
    # Updates progress bar.
    
  }
  
  if(progress_bar) close(progress)
  
  return(list(mcut = ij_min_mcut, mcut_min = ij_min_m,
              rcut = ij_min_rcut, rcut_min = ij_min_r,
              ncut = ij_min_ncut, ncut_min = ij_min_n,
              ccut = ij_min_ccut, ccut_min = ij_min_c))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
  
}

xist_gomoryhu <- function(g, locmax_by_degree=FALSE, progress_bar=TRUE){
  
  n <- length(V(g))
  
  ij_min_r <- ij_min_n <- ij_min_c <- Inf
  ij_min_rcut <- ij_min_ncut <- ij_min_ccut <- c()
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  if(locmax_by_degree) locmaxvec <- which(sapply(1:n, function(i)all(sapply(V(g)[.nei(i)], function(j)V(g)$degree[i] >= V(g)$degree[j]))))
  else locmaxvec <- which(sapply(1:n, function(i)all(sapply(V(g)[.nei(i)], function(j)V(g)$sample[i] >= V(g)$sample[j]))))
  # Computes the set of local maxima of g.
  
  print(locmaxvec)
  
  #locmaxvec <- 1:n
  N <- length(locmaxvec)
  #mcvec <- c()
  V(g)$realvertices <- lapply(1:n, list)
  g_mg <- g
  #xcuts_list <- list()
  
  #shift_a_back_by_b <- function(a,b) a + sapply(a, function(i)sum(b - 1:length(b) <= i))
  # Shifts a list (of vertices) a back to their original position in 1:length(b) assuming
  # they were a subset of the vector obtained by removing b from 1:length(b).
  
  if(progress_bar) progress <- txtProgressBar(min = 1, max = N-1, style = 3)
  
  
  cut_and_split_graph <- function(g_mg, lmv, xcuts_list){
    
    stcut <- min_cut(g_mg, source = lmv[1], target = lmv[2], capacity = E(g_mg)$weight, value.only = FALSE)
    #ij_cut <- min_cut(g, source = s, target = t, capacity = E(g)$weight, value.only = FALSE)
    S_st <- stcut$partition1
    
    #mcvec <- c(mcvec, stcut$value)
    
    #if(!(lmv[1] %in% S_st)) S_st  <- V(g)[-S_st]
    S_st_C <- V(g_mg)[-S_st]
    # Ensures that S_st contains lmv[1] and S_st_C contains lmv[2].
    
    S_real <- unlist(V(g_mg)[S_st]$realvertices)
    print(paste("Current lmv:",toString(lmv)))
    print(V(g_mg)[S_st]$realvertices[1:3])
    print("---------------")
    
    current_xcuts <- xcut_partition_S(S_real, g)
    xcuts_list <- list(mcut_min = c(xcuts_list$mcut_min, current_xcuts$mcut_min),
                       rcut_min = c(xcuts_list$rcut_min, current_xcuts$rcut_min),
                       ncut_min = c(xcuts_list$ncut_min, current_xcuts$ncut_min),
                       ccut_min = c(xcuts_list$ccut_min, current_xcuts$ccut_min),
                       partition = c(xcuts_list$partition, list(S_real)))
    # Computes the XCut of the partition S_ij.
    
    
    #S_st_real1 <- intersect(unlist(V(g_mg)[S_st]$realvertices), S)
    #S_st_real2 <- intersect(unlist(V(g_mg)[S_st_C]$realvertices), S)
    
    lmv_S1 <- which(S_st %in% lmv) + 1
    lmv_S2 <- which(S_st_C %in% lmv) + 1
    #lmv_S1 <- intersect(S_st, lmv)
    #lmv_S2 <- intersect(S_st_C, lmv)
    
    if(progress_bar) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
    # Updates progress bar.
    
    if(length(lmv_S1) > 1){
      
      mergelist1 <- rep(1,length(V(g_mg))) %>% (function(a){a[S_st] <- 2:(length(S_st)+1);a})
      #       cat("mergelist:",mergelist,"\n")
      g_mg1 <- contract(g_mg, mergelist1, vertex.attr.comb=list(degree="sum", sample="sum", realvertices="concat", "first"))
      # Contracts all vertices in vertices_to_merge into one (namely locmaxvec[1]).
      V(g_mg1)$realvertices <- lapply(V(g_mg1)$realvertices, unlist)
      
      g_mg1 <- simplify(g_mg1, edge.attr.comb=list(capacity="sum", weight="sum", width="sum", distance="min", "first"))
      # Collapses all resulting multiple edges into one by summing their weights etc., and deletes loops, thus finishing the merging.
      
      #if(1 %in% S_st_C) merged_vertices <- c(merged_vertices, )
      
      xcuts_list <- cut_and_split_graph(g_mg1, lmv_S1, xcuts_list)
      
    }
    
    if(length(lmv_S2) > 1){
      
      mergelist2 <- rep(1,length(V(g_mg))) %>% (function(a){a[S_st_C] <- 2:(length(S_st_C)+1);a})
      #       cat("mergelist:",mergelist,"\n")
      g_mg2 <- contract(g_mg, mergelist2, vertex.attr.comb=list(degree="sum", sample="sum", realvertices="concat", "first"))
      # Contracts all vertices in vertices_to_merge into one (namely locmaxvec[1]).
      V(g_mg2)$realvertices <- lapply(V(g_mg2)$realvertices, unlist)
      
      g_mg2 <- simplify(g_mg2, edge.attr.comb=list(capacity="sum", weight="sum", width="sum", distance="min", "first"))
      # Collapses all resulting multiple edges into one by summing their weights etc., and deletes loops, thus finishing the merging.
      
      #if(1 %in% S_st_C) merged_vertices <- c(merged_vertices, )
      
      xcuts_list <- cut_and_split_graph(g_mg2, lmv_S2, xcuts_list)
      
    }
    
    return(xcuts_list)
    
    #return(list(g1=, g2, xcuts = ij_xcut, bla=ij_cut$value))
    
  }
  
  
  xcuts_list <- cut_and_split_graph(g, locmaxvec, list())
  
  mc_idx <- which.min(xcuts_list$mcut_min)
  rc_idx <- which.min(xcuts_list$rcut_min)
  nc_idx <- which.min(xcuts_list$ncut_min)
  cc_idx <- which.min(xcuts_list$ccut_min)
  
  #print(xcuts_list)
  print(length(xcuts_list))
  
#  for(ij_xcut in xcuts_list){
#    print(ij_xcut[3]$ncut_min)
  
#    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
#      ij_min_n <- ij_xcut$ncut_min
#      ij_min_ncut <- ij_xcut$partition
#      #        print(paste(locmaxvec[1],"|",locmaxvec[i]))
#      #        print(S_ij)
#    }
#    if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
#      ij_min_c <- ij_xcut$ccut_min
#      ij_min_ccut <- ij_xcut$partition
#    }
#    if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
#      ij_min_r <- ij_xcut$rcut_min
#      ij_min_rcut <- ij_xcut$partition
#    }
  
#  }
  
  if(progress_bar) close(progress)
  
  return(list(mcut = xcuts_list$partition[[mc_idx]], mcut_min = xcuts_list$mcut_min[mc_idx],
              rcut = xcuts_list$partition[[rc_idx]], rcut_min = xcuts_list$rcut_min[rc_idx],
              ncut = xcuts_list$partition[[nc_idx]], ncut_min = xcuts_list$ncut_min[nc_idx],
              ccut = xcuts_list$partition[[cc_idx]], ccut_min = xcuts_list$ccut_min[cc_idx], bla=c()))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
  
}


xcut_multipartition_S <- function(g, partitions, normalize=FALSE){
  
  partitions <- as.integer(factor(partitions))
  
  k <- max(partitions)
  
  if(k==1) return(list(mcut_min=Inf, rcut_min=Inf, ncut_min=Inf, ccut_min=Inf))
  
  mcut_min <- rcut_min <- ncut_min <- ccut_min <- 0
  
  for(i in 1:k){
    
    S <- which(partitions == i)
    
    S_C <- V(g)[-S]
    
    if(F){ #!any(is.null(V(g)$degree))){
      
      S_vol <- sum(V(g)[S]$degree)
      S_C_vol <- sum(V(g)[S_C]$degree)
      
    } else {
      
      S_vol <- sum(E(g)[.inc(S)]$weight) + sum(E(g)[S %--% S]$weight)
      S_C_vol <- sum(E(g)[.inc(S_C)]$weight) + sum(E(g)[S_C %--% S_C]$weight)
      
    }
    # Computing the volumes of S and S_C (also computable using sum(E(g)[.inc(S)]$weight) if vertex degree attribute is missing)
    
    if(normalize) edgesum <- sum(E(g)$weight)
    mc_S <- sum(E(g)[S %--% S_C]$weight) / ifelse(normalize, edgesum, 1)
    
    mcut_min <- mcut_min + mc_S
    rcut_min <- rcut_min + mc_S / (length(S) * length(S_C)) * ifelse(normalize, length(V(g)), 1)
    ncut_min <- ncut_min + mc_S / (S_vol * S_C_vol) * ifelse(normalize, (S_vol + S_C_vol)^2, 1)
    ccut_min <- ccut_min + mc_S / min(S_vol, S_C_vol) * ifelse(normalize, (S_vol + S_C_vol), 1)
    # Computes the various XCuts of S.
    
  }
  
  return(list(mcut_min=mcut_min/2, rcut_min=rcut_min/2, ncut_min=ncut_min/2, ccut_min=ccut_min/2))
  # Returns a list of the computed XCut values of S.
}
# TODO description


generate_tif_grid_graph <- function(sample_raster, M, t=1, sigm=NA, weights_by_sample=TRUE){
  
  if(length(M) > 2) simpleError("M must be a numeric vector of length 1 or 2")
  
  sample <- as.matrix(raster::aggregate(sample_raster, fact=dim(sample_raster)[1:2]/M, fun=sum, expand=FALSE))
  # Discretizes sample_raster to matrix of dimensions M.
  
  sample <- sample / sum(sample)
  # Normalizes sample.
  
  return(generate_grid_graph(sample=sample, t=t, sigm=sigm, weights_by_sample=weights_by_sample))
  # Returns the grid graph built from sample.
}
# TODO description
# Input:
# Output:

# Discretizes and computes the NCut of a raster image (here: of a confocal .tif cell sample).
# M needs to divide the number of pixels in each direction without remainder, and the image has to be quadratic.
# Input: raster image cf, grid size M, distance t.
# Output: List containing an extended graph object, NCut partition and discretized sample of the confocal sample


xvst_part_of_graph <- function(partition, g, colour, locmax_by_degree=locmax_by_degree, draw_plot=TRUE){
  g1 <- delete_vertices(g, V(g)[-partition])
  V(g1)$label <- as.character(1:length(partition))
  if(length(V(g1)) > 1) g1_xcut <- xist_gusfield(g1, locmax_by_degree=locmax_by_degree, progress_bar=FALSE)
  else g1_xcut <- list(mcut=c(),rcut=c(),ncut=c(),ccut=c(),mcut_min=Inf,rcut_min=Inf,ncut_min=Inf,ccut_min=Inf)
  print(paste("JUST CUT:",ncg1_xcut$ncut_min))
  real_stuff <- sort((as.numeric(V(g)$label))[partition])[g1_xcut$ncut]
  V(g)[real_stuff]$color <- colour
  if(draw_plot) plot(g, vertex.label=NA, edge.width = 10*edge_weights/max(edge_weights)+1, vertex.size = V(g)$vertex.size)
  return(list(g=g, partition1=real_stuff, partition2=setdiff(partition, real_stuff), xcut=ncg1_xcut))
}
# Divides the partition nodes of the graph g$graph in two using xcut_via_stMinCuts and colours them appropriately (for now, using NCut).
# Returns the coloured graph, and the two partitions, allowing them to be divided further by another call of the function.


xist_iterated <- function(g, k, cut="ncut", locmax_by_degree=FALSE, normalize=FALSE, colourvec=NULL){
  
  cutnr <-  switch(cut, mcut=1, rcut=2, ncut=3, ccut=4, custom=5)
  # cutnr denotes the number of the desired cut; the order is given by the output of xist_gusfield.
  # "custom" is a placeholder for custom cuts to be added to xist_gusfield.
  # Note that this code is bad practice and should really be done via if statements; however, I am lazy.
  
  #simpleError('cut' must be either one of 'rcut', 'ncut' or 'ccut'.)
  # Computes the appropriate XCut of the discretized confocal sample.
  labels_list <- list(V(g))
  partitions_list <- list(rep(1,length(V(g))))
  #p_glist <- list(V(g))
  #ncut_vec <- c()
  index_best <- 1
  if(is.null(colourvec)) colourvec <- sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)], k)
  #  V(g)$color <- colourvec[1]
  first_cut <- xist_gusfield(g, locmax_by_degree=locmax_by_degree, normalize=normalize, progress_bar=FALSE)
  #cut_value_vec <- c(first_cut[[2*cutnr]])
  cut_list <- list(first_cut)#list(first_cut[[2*cutnr-1]])
  glist <- list(g)
  for(i in 1:(k-1)){
    g_current <- glist[[index_best]]
    partition_current <- cut_list[[index_best]][[2*cutnr-1]]
    #print(labels_list[[index_best]])
    #print(partition_current)
    labels_list[[i+1]] <- labels_list[[index_best]][-partition_current]
    labels_list[[index_best]] <- labels_list[[index_best]][partition_current]
    glist[[i+1]] <- delete_vertices(g_current, V(g_current)[partition_current])
    glist[[index_best]] <- delete_vertices(g_current, V(g_current)[-partition_current])
    cut_list[[i+1]] <- xist_gusfield(glist[[i+1]], locmax_by_degree=locmax_by_degree, normalize=normalize, progress_bar=FALSE)
    cut_list[[index_best]] <- xist_gusfield(glist[[index_best]], locmax_by_degree=locmax_by_degree, normalize=normalize, progress_bar=FALSE)
    
    partitions_list[[i]] <- (partitions_list[[max(1,i-1)]] %>% (function(a){a[labels_list[[i+1]]] <- i+1;a}))
    
    index_best <- which.min(sapply(1:(i+1), function(j)cut_list[[j]][[2*cutnr]]))
    #    V(g)[latest_cut$partition1]$color <- colourvec[i+1] # Colours the original graph
    #labels_list[c(index_best, i+1)] <- list(part1_current, part2_current)
    #partitions_list[[i]] <- partitions
    #glist[i+1] <- delete_vertices(g_current, V(g_current)[-g_current$partition])
    #glist[index_best] <- delete_vertices(g, V(g)[g_current$partition])
    #glist[c(index_best, i+1)] <- list(xvst_part_of_graph(latest_cut$partition1, g, "black", draw_plot=FALSE), xvst_part_of_graph(latest_cut$partition2, g, "black", draw_plot=FALSE))
    #ncut_vec[c(index_best, i+1)] <- c(glist[[index_best]]$xcut$ncut_min, glist[[i+1]]$xcut$ncut_min)
    #index_best <- which.min(cut_vec)
    #print(paste("Current best index:",index_best))
    #print(sapply(cut_list, function(j)j[[2*cutnr-1]]))
    #print(sapply(cut_list, function(j)j[[2*cutnr]]))
    #print(sapply(glist, function(j)length(V(j))))
    #print(labels_list)
    #print("----------------------------------")
  }
  #for(i in 1:k) V(g)[labels_list[[i]]]$color <- colourvec[i]
  V(g)$color <- colourvec[partitions_list[[k-1]]]
  
  #V(g)$vertex.size <- 15*V(g)$sample/max(V(g)$sample) + 3
  #edge_weights <- E(g)$weight
  #  par(cex=0.7, mai=c(0.5,0.5,0.3,0.3))
  #  plot(g, vertex.label=NA, edge.width = 10*edge_weights/max(edge_weights)+1, vertex.size = V(g)$vertex.size)
  return(list(g=g, partitions=partitions_list))
}


M <- 16
t <- 1.5
# Defining grid size M and neighbourhood distance t.

setwd("./RExercises/PhD testing range/")
#setwd("../Schreibtisch/Masterarbeit/R files/Data")
#setwd("/path/to/NIH3T3/data")
# Setting the correct path to the NIH3T3 dataset

cf1_13 <- raster::raster("./Confocal_still_images/s_C001Tubulin13.tif")
# Uses the "raster" package to read the .tif file into a raster object.

g <- generate_tif_grid_graph(cf1_13, M=128, t=1.5)
#g_sample <- V(g)$sample

collist <- c("orangered1", "lawngreen", "yellow", "purple","green4", "darkgoldenrod2", "blue", "cyan", "brown3", "magenta2")

xc_g_7 <- xist_iterated(g, 10, normalize=TRUE, colourvec=collist)
plot_graph_partition(xc_g_7$g, collist[xc_g_7$partitions[[6]]])
sapply(1:k, function(i)sum(xc_g_7$partitions[[k-1]]==i))

spec_g_7 <- manual_multiway_spec_clustering(g, k=10)

a <- xcut_via_stmincut(g)
g1 <- delete.vertices(g, a$ncut)
g2 <- delete.vertices(g, V(g)[-a$ncut])
a1 <- xcut_via_stmincut(g1)
a2 <- xcut_via_stmincut(g2)


#cf1_13 <- raster::raster("./Confocal_still_images/s_C001Tubulin13.tif")

#cf1_13_xcut_iterated <- xvst_iterated(k=7, pt$g, collist)
#stop_time <- Sys.time()
#print(stop_time - start_time)
# Dividing the cell image cf1_13 into k partitions, colouring them according to collist.

colr <- colorRampPalette(RColorBrewer::brewer.pal(9, 'YlOrRd'))
rasterVis::levelplot(flip(t(cf1_13), 2),
                     margin = FALSE,
                     colorkey = F,#list(space='bottom', labels=list(at=seq(0,5000,length.out=6))),
                     par.settings = list(axis.line=list(col='transparent')),
                     scales = list(draw=FALSE),
                     col.regions = colr,
                     at = seq(0,5000,length.out=101))
# Plotting the cell sample cf1_13 with proper scaling and nice colouring.

put_square_at_partition_KEEEEEP <- function(partitionvec, pixels, M){
  
  if(length(M)==1) M <- c(M, M)
  if(length(pixels)==1) pixels <- c(pixels, pixels)
  # Input can be provided in a simplified form if the underlying image (and resulting grid) should be square.
  
  locations <- cbind(pixels[1] / M[1] * ((partitionvec-1) %/% M[2] + 0.5),
                     pixels[2] / M[2] * ((partitionvec-1) %% M[2] + 0.5))
  # Computes the locations of the centers of the rectangles corresponding to the indices in partitionvec,
  # given that it is based on a pixels[1] x pixels[2] grid that is discretized onto an M[1] x M[2] sized grid.
  
  return(rgeos::gBuffer(SpatialPoints(locations), width=0.5*pixels[1]/M[1], quadsegs=1, capStyle="SQUARE"))
}


plot_partition_overlay <- function(sample_raster, partition, colourvec){
  
  k <- length(unique(c(partition)))
  pixels <- dim(sample_raster)[1:2]
  M <- dim(partition)
  if(length(M) != 2) simpleError("partition must be a 2d matrix")
  
  #if(pixels[1] * pixels[2] != length(partition)) simpleError("partition does not fit filename image size")
  
  
  #if(length(M)==1) M <- c(M, M)
  #if(length(pixels)==1) pixels <- c(pixels, pixels)
  # Input can be provided in a simplified form if the underlying image (and resulting grid) should be square.
  
  colr2 <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Greys'))
  
  rasterVis::levelplot(flip(t(sample_raster), 2), margin = FALSE, colorkey=FALSE, par.settings = list(axis.line=list(col='transparent')),
                       scales = list(draw=FALSE), col.regions = colr2) +
    latticeExtra::layer(lapply(1:k, function(i){
    
      partitionvec <- which(partition==i)
      #print(partitionvec)
      
      locations <- cbind(pixels[1] / M[1] * ((partitionvec-1) %/% M[2] + 0.5),
                         pixels[2] / M[2] * ((partitionvec-1) %% M[2] + 0.5))
      # Computes the locations of the centers of the rectangles corresponding to the indices in partitionvec,
      # given that it is based on a pixels[1] x pixels[2] grid that is discretized onto an M[1] x M[2] sized grid.
      #print(i)
      #print(dim(locations))
      buffer_locations <- rgeos::gBuffer(SpatialPoints(locations), width=0.5*pixels[1]/M[1], quadsegs=1, capStyle="SQUARE")
      
      return(sp::sp.polygons(buffer_locations, fill=colourvec[i], col=NA, alpha=0.5))
    
  }), data=list(k=k,partition=partition,M=M,pixels=pixels,colourvec=colourvec))
  
}
# Returns SpatialPolygons for squares around the locations specified by the discretized locations given by partitionvec and hM to overlay on top of an image of hpixels x vpixels
# Currently only computes squares. WIP.

cl <- collist[c(7,8,4,1,2,3,6,3,5,9)]
cl <- collist[c(1,2,3,4,8,6,10,5,9,7)]
plot_partition_overlay(cf1_13, matrix(xc_g_7$partitions[[9]], 128), cl)
#plot_partition_overlay(cf1_13, matrix(spec_g_7[[7-1]], 128), cl)



uncoloured_squares_list <- lapply(cf1_13_xcut_iterated$partitions, function(p)put_square_at_partition(p, dim(cf1_13)[1], M))
colr2 <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Greys'))
rasterVis::levelplot(flip(t(cf1_13), 2),
                     margin = FALSE, colorkey=FALSE,
                     par.settings = list(axis.line=list(col='transparent')),
                     scales = list(draw=FALSE),
                     col.regions = colr2) +
  latticeExtra::layer(lapply(1:k, function(idx)sp.polygons(uncoloured_squares_list[[idx]], fill=collist[idx], col=NA, alpha=0.5)))
# Plots the image together with the partitions (which are coloured appropriately)




compute_weight_matrix <- function(g){
  wm_g <- matrix(rep(0,length(V(g))^2), length(V(g)))
  wm_g[ends(g, E(g))] <- E(g)$weight
  return(wm_g + t(wm_g))
}
wm_gg <- compute_weight_matrix(gg)
gg2 <- generate_tif_grid_graph("./Confocal_still_images/s_C001Tubulin13.tif", M=32, t=1.5)
g_spec_8 <- sClust::VonLuxburgSC(compute_weight_matrix(gg), K=5)
plot_graph_partition(gg, g_spec_8$cluster)
#wm_g <- outer(V(g), V(g), function(i,j)ifelse())
#sClust::checking.gram.similarityMatrix()
a_spec <- manual_multiway_spec_clustering(g, k=8)
plot_graph_partition(g, collist[a_spec[[8-1]]])
plot_partition_overlay(cf1_13, matrix(a_spec[[8-1]], 128), collist)

manual_multiway_spec_clustering <- function(g, k=7){
  
#  wm_g <- compute_weight_matrix(g)
#  deg_mat_inv <- diag(1/rowSums(wm_g))
#  L_mat_norm <- diag(length(V(g))) - deg_mat_inv %*% wm_g
  
  L_mat_norm <- diag(rep(0,length(V(g)))) #diag(length(V(g))) #matrix(rep(0,length(V(g))^2), length(V(g)))
  edgemat <- ends(g, E(g))
  L_mat_norm[rbind(edgemat, cbind(edgemat[,2],edgemat[,1]))] <- - E(g)$weight
  L_mat_norm <- - L_mat_norm / rowSums(L_mat_norm)
  L_mat_norm[cbind(V(g), V(g))] <- 1 #- rowSums(L_mat_norm)
  
#  print(L_mat_norm[rbind(c(1,2),c(1,17),c(2,3),c(5,5))])
#  print(L_mat_norm2[rbind(c(1,2),c(1,17),c(2,3),c(5,5))])
  # Computes the normalized graph Laplacian L_mat_norm = I - D^-1 W for weight matrix W and degree matrix W.
  
  #  L_mat_unnorm <- diag(rowSums((g$weight_mat + t(g$weight_mat)))) - (g$weight_mat + t(g$weight_mat))
  # Computes the inverse deg_mat_inv of the degree matrix D and both graph Laplacians L_mat_norm and L_mat_unnorm.
  
  try({
    L_mat_norm_eigen <- eigen(L_mat_norm)
    #    L_mat_unnorm_eigen <- eigen(L_mat_unnorm)
    # Eigenvalues and eigenvectors of the Laplacians.
    
    U_mat_norm <- L_mat_norm_eigen$vectors[,(ncol(L_mat_norm_eigen$vectors)-k+1):ncol(L_mat_norm_eigen$vectors)]
    #    U_mat_unnorm <- L_mat_unnorm_eigen$vectors[,(ncol(L_mat_unnorm_eigen$vectors)-k+1):ncol(L_mat_unnorm_eigen$vectors)]
    # Computes the matrix U (i.e. the k smallest eigenvectors) for the (Un)Normalized spectral clustering algorithm.
    #save(U_mat_norm, file="umatnorm1.RData")
    #print(dim(U_mat_norm))
    #print(lapply(2:k, function(i)kmeans(U_mat_norm[,(k-i+1):k], centers=i)$cluster))
    set.seed(8472)
    return(lapply(2:k, function(i)kmeans(U_mat_norm[,(k-i+1):k], centers=i)$cluster)) #return(c(norm=lapply(...)))
    #                  kmeans(U_mat_norm[,(k-1):k], centers=k-5, nstart=dim(U_mat_norm)[1]-5)$cluster,
    #                norm3 = kmeans(U_mat_norm[,(k-2):k], centers=k-4, nstart=dim(U_mat_norm)[1]-4)$cluster,
    #                norm4 = kmeans(U_mat_norm[,(k-3):k], centers=k-3, nstart=dim(U_mat_norm)[1]-3)$cluster,
    #                norm5 = kmeans(U_mat_norm[,(k-4):k], centers=k-2, nstart=dim(U_mat_norm)[1]-2)$cluster,
    #                norm6 = kmeans(U_mat_norm[,(k-5):k], centers=k-1, nstart=dim(U_mat_norm)[1]-1)$cluster,
    #                norm7 = kmeans(U_mat_norm[,(k-6):k], centers=k, nstart=dim(U_mat_norm)[1])$cluster))#,
    #                unnorm = kmeans(U_mat_unnorm, centers=k, nstart=dim(U_mat_unnorm)[1])$cluster))
  }, silent=TRUE)
}
# TODO Add text + combine with above spectral clustering function

# , nstart=dim(U_mat_norm)[1]+i-k




a3 <- xist_gomoryhu(g1, progress_bar=F)
#a3$ncut
plot_graph_partition(g1, a3$ncut, edge_scale=10, vertex_offset=5, vertex_scale=20)
#collist <- c("orangered1", "lawngreen", "yellow", "purple","green4", "darkgoldenrod2", "blue", "cyan", "brown3")
g2 <- g1
V(g2)$rvtest <- sapply(1:length(V(g2)), list)
V(g2)$rvtest[1:4]
g22 <- contract(g2, c(1,1,2:255), vertex.attr.comb=list(degree="sum",sample="sum",rvtest="concat","first"))
V(g22)$rvtest[1:4]
V(g22)$rvtest <- lapply(V(g22)$rvtest, unlist)
V(g22)$rvtest[1:4]
g222 <- simplify(g22, edge.attr.comb=list(capacity="sum", weight="sum", width="sum", distance="min", "first"))
V(g222)$rvtest[1:4]
g23 <- contract(g222, c(1,1,2:254), vertex.attr.comb=list(degree="sum",sample="sum",rvtest="concat","first"))
V(g23)$rvtest[1:4]
V(g23)$rvtest <- lapply(V(g23)$rvtest, unlist)
V(g23)$rvtest[1:4]


g1 <- discretize_and_xcut_confocal(cf1_13, M=16, t=1, cut="none")$g
a1 <- xist_gusfield(g1, progress_bar=F)
a1$bla
a2 <- xvst_on_local_maxima(g1, progress_bar=F)
a2$bla



if(F){

aa <- matrix(rep(NA,20),nrow=2)
for(n in 10*(1:10)){
  st <- Sys.time()
  replicate(20, generate_grid_graph(t=1.5, matrix(1:n^2,nrow=n)))
  et <- Sys.time()
  aa[1,n/10] <- (et-st)/20
  gc()
  st <- Sys.time()
  replicate(20, generate_graph_OLD(t=1.5, matrix(1:n^2,nrow=n)))
  et <- Sys.time()
  aa[2,n/10] <- (et-st)/20
  gc()
}

}



# NOTATION: As in the thesis, t denotes the neighbourhood distance of a graph and n the sample size of a sample X
# Unlike the thesis, M is used here to denote the grid size m from the thesis. I apologize for the confusion.


euclidean_distance <- function(a, b) return(sqrt(sum((a-b)^2)))
# Computes the euclidean distance between vectors a and b.

t_radial_kernel <- function(x,y,t) return(ifelse(sum((x-y)^2)<=t, 1, 0))
# Realizes the t-radial kernel function: Returns 1 if the vectors x and y are no more than t euclidean distance
# apart. Otherwise, it returns 0.

generate_graph_OLD <- function(t, sample_discretized, locmax_degree=FALSE){
  
  M <- dim(sample_discretized)[1]
  sample_discretized <- matrix(as.numeric(sample_discretized), nrow=M)
  # Determines grid size M. Assumes that sample_discretized is an M x M matrix. Also sets all values in
  # sample_discretized to be double precision (instead of integer) to avoid receiving NAs by multiplication (see below).

  all_squares <- rbind(rep(1:M, each=M), rep(1:M, M))
  str_squares <- sapply(1:M^2, toString)
  # Determines coordinates and labels for the nodes of the graph.
  
  adj_mat <- matrix(rep(0,M^4), M^2, M^2)
  edge_list <- c()
  capacity_list <- c()
  distance_list <- c()
  removed_vertices <- c()
  local_maximum <- rep(TRUE, M^2)
  # Initializing. adj_mat will be the weight (or adjacency) matrix of the generated graph, edge list will be a vector
  # containing all edges, capacity_list will contain the same weights for using min_cuts, and removed_vertices will
  # consist of all vertices that are removed due to not containing any observations, i.e. sample_discretized[i] = 0
  # for some vertex (index) i. local_maximum[i] indicates whether vertex i is a local maximum, i.e. it consists of more
  # observations than any of its neighbouring vertices.
  
  for(i in 1:(M^2-1)){
    
    if(!is.na(sample_discretized[all_squares[1,i], all_squares[2,i]]) &
              sample_discretized[all_squares[1,i], all_squares[2,i]] > 0){
      
      for(j in (i+1):M^2){
        
        if(sample_discretized[all_squares[1,j], all_squares[2,j]] > 0
           & euclidean_distance(all_squares[,i], all_squares[,j]) <= t){
        # Are nodes i and j connected?
          
          edge_list <- c(edge_list, i,j)
          # Contains all edges (i1, j1), (i2, j2), ... as a vector like so: {i1, j1, i2, j2, ...}. Format needed for igraph.
          
          distance_list <- c(distance_list, euclidean_distance(all_squares[,i], all_squares[,j]))
          
          capacity_list <- c(capacity_list, sample_discretized[all_squares[1,i], all_squares[2,i]] *
                                            sample_discretized[all_squares[1,j], all_squares[2,j]])
          adj_mat[i,j] <- sample_discretized[all_squares[1,i], all_squares[2,i]] *
                          sample_discretized[all_squares[1,j], all_squares[2,j]]
          # For an explanation of adj_mat and capacity_list see above. The weights are given by w_{ij} = Y_i * Y_j.
          
          if(sample_discretized[all_squares[1,i], all_squares[2,i]] < sample_discretized[all_squares[1,j], all_squares[2,j]])
            local_maximum[i] <- FALSE
          if(sample_discretized[all_squares[1,i], all_squares[2,i]] > sample_discretized[all_squares[1,j], all_squares[2,j]])
            local_maximum[j] <- FALSE
        }
      }
    }
    else removed_vertices <- c(removed_vertices, i)
    # If a vertex i satisfies Y_i = 0, it can be removed since all edges containing i have weight 0.
  }
  if(is.na(sample_discretized[all_squares[1,M^2], all_squares[2,M^2]]) |
     sample_discretized[all_squares[1,M^2], all_squares[2,M^2]] == 0)
    removed_vertices <- c(removed_vertices, M^2)
  # The last vertex (M, M) was not covered by the for loop. This is rectified here.
  
  g <- graph(edges=edge_list, n=M^2, directed = FALSE)
  V(g)$label = 1:M^2
  g <- delete_vertices(g, removed_vertices)
  # Creates a graph with M^2 vertices and edges from edge_list, then removes redundant nodes.
  
  E(g)$capacity <- capacity_list
  # Assigns capacities (=weights) of the edges from edge_vector. This is needed for igraph::min_cut.
  
  E(g)$distance <- distance_list
  
  E(g)$weight <- 100 * capacity_list / sum(capacity_list)
  E(g)$width <- 10*E(g)$weight/max(E(g)$weight)+1
  # This is used solely for plotting purposes. To balance edge thickness, their capacities are normalized and scaled.
  
  if(!is.null(removed_vertices)){
    V(g)$x <- all_squares[1,-removed_vertices]
    V(g)$y <- all_squares[2,-removed_vertices]
    V(g)$label <- str_squares[-removed_vertices]
    V(g)$sample <- as.numeric(t(sample_discretized))[-removed_vertices]
    V(g)$local_maximum <- local_maximum[-removed_vertices]
    V(g)$degree <- sapply((1:M^2)[-removed_vertices], function(i)sum(adj_mat[i,-removed_vertices]))
    adj_mat <- adj_mat[-removed_vertices, -removed_vertices]
  }
  else{
    V(g)$x <- all_squares[1,]
    V(g)$y <- all_squares[2,]
    V(g)$label <- str_squares
    V(g)$sample <- as.numeric(t(sample_discretized))
    #V(g)$size <- 
    V(g)$local_maximum <- local_maximum
    V(g)$degree <- sapply(1:M^2, function(i)sum(adj_mat[i,]))
  }
  # Assigning coordinates and labels to all vertices but those removed due to not containing any observation.
  
  adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
  # Since it is known that adj_mat is symmetric, this line copies all weights from the upper triangular part of the matrix
  # to the bottom one, ensuring symmetry.
  
#  if(locmax_degree){
#    degreevec <- sapply(1:M^2, function(i)sum(adj_mat[i,]))
#    for(i in 1:M^2) V(g)$degree[i] <- ifelse(all(adj_mat[i,]==0 | degreevec[i] >= degreevec), degreevec[i], 0)
#    #print(V(g)$local_maximum)
#  }
  
  return(list(graph=g, weight_mat=adj_mat, removed_vertices=removed_vertices))
  # Returns a list consisting of the graph g, associated weight matrix weight_mat
  # and a list of vertices removed due to not containing any observations
}
# Generates an extended graph object containing a graph constructed using t-radius neighbourhood,
# associated weight matrix with weights w_ij = Y_i * Y_j for observations Y_i and Y_j if i~j,
# and vector removed_vertices containing vertices i removed due to Y_i = 0.
# Input: Distance t and square matrix sample_discretized (likely containing a discretized sample).
# Output: list(graph, weight_mat, removed_vertices)

generate_sample_graph_OLD <- function(t, raw_sample){
  
  n <- dim(raw_sample)[1]
  # Determines sample size n. Assumes that sample_discretized is an n x 2 matrix.
  
  adj_mat <- matrix(rep(0,n^2), n, n)
  edge_list <- c()
  # Initializing. adj_mat will be the weight (or adjacency) matrix of the generated graph, edge list will be a vector
  # containing all edges and capacity_list will contain the same weights for using min_cuts
  
  for(i in 1:(n-1)){
    
      for(j in (i+1):n){
        
        if(euclidean_distance(raw_sample[i,], raw_sample[j,]) <= t){
          # Are nodes i and j connected?
          
          edge_list <- c(edge_list, i,j)
          # Contains all edges (i1, j1), (i2, j2), ... as a vector like so: {i1, j1, i2, j2, ...}. Format needed for igraph.
          
          adj_mat[i,j] <- 1
          # The weights are given by w_{ij} = 1 for the t-neighbourhood graph.

        }
        
      }
  }
  
  g <- graph(edges=edge_list, n=n, directed = FALSE)
  # Creates a graph with n vertices and edges from edge_list.
  
  E(g)$capacity <- rep(1, length(edge_list)/2)
  # Assigns capacities (=weights) of the edges from edge_list. This is needed for igraph::min_cut.
  
  E(g)$weight <- rep(100/n, length(edge_list)/2)
  E(g)$width <- rep(10, length(edge_list)/2)
  # This is used solely for plotting purposes. To balance edge thickness, their capacities are normalized and scaled.
  
  V(g)$x <- raw_sample[,1]
  V(g)$y <- raw_sample[,2]
  V(g)$label <- sapply(1:n, toString)
  V(g)$sample <- rep(1, n)
  # Assigning coordinates and labels to all vertices.
  
  adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
  # Since it is known that adj_mat is symmetric, this line copies all weights from the upper triangular part of the matrix
  # to the bottom one, ensuring symmetry.
  
  return(list(graph=g, weight_mat=adj_mat, removed_vertices=NULL))
  # Returns a list consisting of the graph g, associated weight matrix weight_mat
  # and a list of vertices removed due to not containing any observations
}
# Generates an extended graph object containing a graph constructed using t-radius neighbourhood,
# associated weight matrix with weights w_ij = 1 for observations X_i and X_j if i~j.
# Input: Distance t and a matrix with two columns containing the two-dim. sample.
# Output: list(graph, weight_mat, removed_vertices=NULL)

manual_spec_clustering <- function(g, method="kmeans"){
  
  if(is.null(g$removed_vertices)){
    deg_mat_inv <- diag(1/rowSums((g$weight_mat + t(g$weight_mat))))
    L_mat_norm <- diag(length(V(g$graph))) - deg_mat_inv %*% (g$weight_mat + t(g$weight_mat))
    L_mat_unnorm <- diag(rowSums((g$weight_mat + t(g$weight_mat)))) - (g$weight_mat + t(g$weight_mat))
  } else{
    deg_mat_inv <- diag(1 / rowSums((g$weight_mat + t(g$weight_mat))[-g$removed_vertices, -g$removed_vertices]))
    L_mat_norm <- diag(length(V(g$graph))) -
                  deg_mat_inv %*% (g$weight_mat + t(g$weight_mat))[-g$removed_vertices, -g$removed_vertices]
    L_mat_unnorm <- diag(rowSums((g$weight_mat + t(g$weight_mat))[-g$removed_vertices, -g$removed_vertices])) -
                    (g$weight_mat + t(g$weight_mat))[-g$removed_vertices, -g$removed_vertices]
  }
  # Computes the inverse deg_mat_inv of the degree matrix D and both graph Laplacians L_mat_norm and L_mat_unnorm.
  
  try({
    L_mat_norm_eigen <- eigen(L_mat_norm)
    L_mat_unnorm_eigen <- eigen(L_mat_unnorm)
    # Eigenvalues and eigenvectors of the Laplacians.
    
    U_mat_norm <- L_mat_norm_eigen$vectors[,(ncol(L_mat_norm_eigen$vectors)-1):ncol(L_mat_norm_eigen$vectors)]
    U_mat_unnorm <- L_mat_unnorm_eigen$vectors[,(ncol(L_mat_unnorm_eigen$vectors)-1):ncol(L_mat_unnorm_eigen$vectors)]
    # Computes the matrix U (i.e. the two smallest eigenvectors) for the (Un)Normalized spectral clustering algorithm.
    
    if(method == "naive"){
      return(list(norm=ifelse(U_mat_norm[,1] >= 0, 1, 2), unnorm=ifelse(U_mat_unnorm[,1] >= 0, 1, 2)))
    }
    if(method == "mean"){
      return(list(norm=ifelse(U_mat_norm[,1] >= mean(U_mat_norm[,1]), 1, 2),
                  unnorm=ifelse(U_mat_unnorm[,1] >= mean(U_mat_unnorm[,1]), 1, 2)))
    }
    if(method == "kmeans"){
      return(list(norm = kmeans(U_mat_norm, centers=2, nstart=dim(U_mat_norm)[1])$cluster,
                  unnorm = kmeans(U_mat_unnorm, centers=2, nstart=dim(U_mat_unnorm)[1])$cluster))
    }else{
      print("Please use either one of 'naive', 'mean' or 'kmeans' as clustering method")
      return(NULL)
    }
  }, silent=TRUE)
}
# Clusters the extended graph object g with both unnormalized and normalized spectral clustering.
# Offers three methods for clustering the components of the first (two) eigenvector(s):
# - "naive": Clusters the components of the second eigenvector according to their sign.
# - "mean": Clusters the components of the second eigenvector according to their relation to its mean.
# - "kmeans": Uses stats::kmeans (clustering) to cluster the components of the first and second eigenvector.
# Input: Extended graph object g and a clustering method 'method'.
# Output: A vector of integers 1 or 2 representing the cluster the nodes are allocated to.

xcut_via_stmincut_OLD <- function(g){
  
  ij_min_r <- Inf
  ij_min_n <- Inf
  ij_min_c <- Inf
  ij_min_rcut <- NULL
  ij_min_ncut <- NULL
  ij_min_ccut <- NULL
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
#  g$weight_mat[g$removed_vertices, ] <- 0
#  g$weight_mat[, g$removed_vertices] <- 0
  # Making sure removed vertices do not contribute to the cut
  
  g_directed <- as.directed(g$graph, mode="mutual")
  # The st-MinCut function igraph::stMincuts only accepts directed graphs.
#  bla <- c()
  for(i in 1:(length(V(g$graph))-1)){
    for(j in (i+1):length(V(g$graph))){
      
      ij_cut <- stMincuts(g_directed, i, j)
      # Computes all partitions attaining the ij-MinCut.
      
      for(S_ij in ij_cut$partition1s){
      # Goes through all partitions attaining the ij-MinCut.
      
 #       S_ij_C <- V(g$graph)[-S_ij]
        # Divides the vertices into the partition S_ij that attained the ij-MinCut ij_cut.
        # Selects the last partition attaining the ij-MinCut. This is to obtain a larger partition.
#        print(S_ij)
#        print(S_ij_C)
        
        ij_xcut <- xcut_partition_S_OLD(S_ij, g)
        # Computes the XCut of the partition S_ij.
        
#        if(isTRUE(all.equal(ij_xcut$ncut_min, ij_min_n))){
#          print(paste(i,j))
#          print(S_ij)
#          print(ij_min_n)
#        }
        
        if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
          ij_min_n <- ij_xcut$ncut_min
          ij_min_ncut <- S_ij
        }
        if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
          ij_min_c <- ij_xcut$ccut_min
          ij_min_ccut <- S_ij
        }
        if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
          ij_min_r <- ij_xcut$rcut_min
          ij_min_rcut <- S_ij
        }
        # Updates the best XCut and the partition that attains it
      }
    }}
#    bla <- c()
#    for(i in 1:(length(V(g$graph))-1)){
#      for(j in (i+1):length(V(g$graph))){
#        ij_cut <- stMincuts(g_directed, i, j)
#        # Computes all partitions attaining the ij-MinCut.
#        
# #       for(S_ij in ij_cut$partition1s){
#    ij_xcut <- xcut_partition_S(S_ij, g)
#    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min == ij_min_n){
#      bla <- c(bla, i, j)
  #  }#}}
  #}
  return(list(rcut = ij_min_rcut, rcut_min = ij_min_r,
              ncut = ij_min_ncut, ncut_min = ij_min_n,
              ccut = ij_min_ccut, ccut_min = ij_min_c))#, bla=bla))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}
partitionToInt <- function(kvec) sum(2^(kvec-1))
intToPartition <- function(k, len=M^2) (1:len)[intToBits(k)[1:len]==1]

xcut_via_stmincut2 <- function(g){

  m <- length(V(g$graph))
  g_directed <- as.directed(g$graph, mode="mutual")
  # The st-MinCut function igraph::stMincuts only accepts directed graphs.
  
  #ij_cuts <- unique(unlist(lapply(1:(m-1), function(i)unlist(lapply((i+1):m, function(j)vapply(stMincuts(g_directed, i, j)$partition1s, partitionToInt, 1))))))
  ij_cuts <- unique(unlist(lapply(1:(m-1), function(i)unique(unlist(lapply((i+1):m, function(j)lapply(stMincuts(g_directed, i, j)$partition1s, function(x){ attributes(x) <- NULL; x })), recursive=FALSE))), recursive=FALSE))
  #ij_cuts <- unique(sapply(2:m, function(j)sapply(stMincuts(g_directed, 1, j)$partition1s, partitionToInt)))
  all_cuts <- sapply(ij_cuts, function(partition){ x <- xcut_partition_S(partition, g); attributes(x) <- NULL; unlist(x) } )
  #print(which(is.na(ij_cuts)))
  #print(ij_cuts)
  #print(all_cuts)
  #xcut_indices <- max.col(structure(vapply(all_cuts, "-", 1), dim=dim(all_cuts)))
  xcut_indices <- vapply(1:4, function(x)which.min(all_cuts[x,]), 1)
  #print(xcut_indices)
  #print(all_cuts)
  #return(all_cuts)
  return(list(mcut = ij_cuts[[xcut_indices[1]]], mcut_min = all_cuts[1, xcut_indices[1]],
              rcut = ij_cuts[[xcut_indices[2]]], rcut_min = all_cuts[2, xcut_indices[2]],
              ncut = ij_cuts[[xcut_indices[3]]], ncut_min = all_cuts[3, xcut_indices[3]],
              ccut = ij_cuts[[xcut_indices[4]]], ccut_min = all_cuts[4, xcut_indices[4]]))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}
# Realizes the XCut-via-st-MinCuts algorithm. Input: Extended graph object returned from generate_graph
# (or generate_cockroach or generate_r_graph). Output: List of XCut values with corresponding partitions.

st <- Sys.time()
rc_vec_st_5_1 <- lapply(t_vec2, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec5, nrow=M2))))
et <- Sys.time()
print(et-st)

xcut_via_stmincut_on_local_maxima_OLD <- function(g, allcuts=TRUE, progress_bar=TRUE){
  
  ij_min_r <- Inf
  ij_min_n <- Inf
  ij_min_c <- Inf
  ij_min_rcut <- NULL
  ij_min_ncut <- NULL
  ij_min_ccut <- NULL
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  g_directed <- as.directed(g$graph, mode="mutual")
  # The st-MinCut function igraph::stMincuts only accepts directed graphs.
  #  bla <- c()
  
  locmaxvec <- sort(which(V(g$graph)$local_maximum))
  # Computes the set of local maxima of g.
  #print(locmaxvec)
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = length(locmaxvec)*(length(locmaxvec)+1)/2, style = 3)
  # Setting up the progress bar.
  
  #for(i in 1:(length(V(g$graph))-1)){
  #  for(j in (i+1):length(V(g$graph))){
  #    if(V(g$graph)$local_maximum[i] & V(g$graph)$local_maximum[j]){
  for(i in locmaxvec){
    for(j in locmaxvec[locmaxvec > i]){
        
      if(allcuts) ij_cut <- stMincuts(g_directed, i, j)
      # Computes all partitions attaining the ij-MinCut.
      
      if(length(ij_cut$partition1s) == 0){
      # stMincuts doesn't return a partition if the MinCut value is zero.
        
        ij_cut <- min_cut(g$graph, source=i, target=j, value.only=FALSE)
        # In that case we use igraph::min_cut to determine the st-MinCut partition.
        
        ij_cut$partition1s <- list(ij_cut$partition1)
        # Stores the partition as a one-element list for the for loop below.
        
      }
      
      for(S_ij in ij_cut$partition1s){
        # Goes through all partitions attaining the ij-MinCut.
        
        ij_xcut <- xcut_partition_S(S_ij, g)
        # Computes the XCut of the partition S_ij.
        
        if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
          ij_min_n <- ij_xcut$ncut_min
          ij_min_ncut <- S_ij
          #print(paste(i,"|",j))
          #print(S_ij)
        }
        if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
          ij_min_c <- ij_xcut$ccut_min
          ij_min_ccut <- S_ij
        }
        if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
          ij_min_r <- ij_xcut$rcut_min
          ij_min_rcut <- S_ij
        }
        # Updates the best XCut and the partition that attains it.
      }
      if(progress_bar) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
      # Updates progress bar.
      
      #print(i)
      #print(j)
      #print(S_ij)
      #print(ij_xcut$ncut_min)
        
    }
  }
  if(progress_bar) close(progress)
  #    bla <- c()
  #    for(i in 1:(length(V(g$graph))-1)){
  #      for(j in (i+1):length(V(g$graph))){
  #        ij_cut <- stMincuts(g_directed, i, j)
  #        # Computes all partitions attaining the ij-MinCut.
  #        
  # #       for(S_ij in ij_cut$partition1s){
  #    ij_xcut <- xcut_partition_S(S_ij, g)
  #    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min == ij_min_n){
  #      bla <- c(bla, i, j)
  #  }#}}
  #}
  return(list(rcut = ij_min_rcut, rcut_min = ij_min_r,
              ncut = ij_min_ncut, ncut_min = ij_min_n,
              ccut = ij_min_ccut, ccut_min = ij_min_c))#, bla=bla))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}
# Realizes the XCut-via-st-MinCuts algorithm, restricted to pairs of "local maxima" (i.e. vertices with more observations than all their neigbours).
# Input: Extended graph object returned from generate_graph (or generate_cockroach or generate_r_graph).
# Output: List of XCut values with corresponding partitions.

ncut_via_stmincut_CURRENTLY_BROKEN <- function(g){
  
  ij_min_n <- Inf
  ij_min_ncut <- NULL
  partitionvec <- 1:length(V(g$graph))
  # partitionvec stores the partitions that have been tested. For example, if an ij-MinCut separates nodes k and l, partitionvec[k] == partitionvec[l].
  
#  edgemat <- as_edgelist(g$graph)
  g_directed <- as.directed(g$graph, mode="mutual")
  # The st-MinCut function igraph::stMincuts only accepts directed graphs.
#  i <- 1
#  j <- 2
  k <- 1
  i <- which.max(V(g$graph)$sample)
  i_before_j <- 0
  while(length(V(g_directed)) > 1){
    #print(paste("Length:",length(V(g_directed))))
#    i1 <- which.max(V(g_directed)$sample)
#    print(V(g_directed)[i1]$sample)
    #print(V(g_directed)[i]$sample)
    j <- which.max(V(g_directed)$sample[-i])
    if(j>=i) j <- j+1
    #print(V(g_directed)[j]$sample)
#    print(paste0("i1=",i1," | i=",i))
#    k <- k+1
#  for(i in 1:(length(V(g$graph))-1)){
#    for(j in (i+1):length(V(g$graph))){
      #print(j)
#      g_directed <- as.directed(g$graph, mode="mutual")
      ij_cut <- stMincuts(g_directed, i, j)
      # Computes all partitions attaining the ij-MinCut.
#      return(g_directed)
      for(S_ij in ij_cut$partition1s){
        # Goes through all partitions attaining the ij-MinCut.
        
#        S_ij_C <- V(g$graph)[-S_ij]
        # Divides the vertices into the partition S_ij that attained the ij-MinCut ij_cut.
        S_ij_retransformed <- sort(unlist(sapply(S_ij, function(x)which(partitionvec==x))))
        
        
        ij_xcut <- xcut_partition_S(S_ij_retransformed, g)
        # Computes the XCut of the partition S_ij.
        #print(ij_xcut$ncut_min)
        
        # Transforms S_ij into its original partition by untangling all combined vertices.
        
        if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
          ij_min_n <- ij_xcut$ncut_min
          ij_min_ncut <- S_ij_retransformed
        }
        # Updates the best XCut and the partition that attains it
        
#        g$weight_mat[i,] <- g$weight_mat[i,] + g$weight_mat[j,]
#        g$weight_mat[,i] <- g$weight_mat[,i] + g$weight_mat[,j]
#        g$weight_mat[i,i] <- 0
#        g$weight_mat <- g$weight_mat[-j,-j]
#        g$graph <- delete_vertices(g$graph, j)
#        V(g$graph)$label <- as.character(1:(length(V(g$graph))-k+1))
#        g_directed <- contract(g_directed, c(1,1:(length(V(g_directed))-1)), vertex.attr.comb=list(label="first", "ignore"))
        #if(i > j) i_before_j <- i_before_j + 1
        if(j > 1 & j < length(V(g_directed)))
          g_directed <- contract(g_directed, c(1:(j-1), i, j:(length(V(g_directed))-1)), vertex.attr.comb=list(label="first", sample="sum", "ignore"))
        if(j == 1)
          g_directed <- contract(g_directed, c(i, 1:(length(V(g_directed))-1)), vertex.attr.comb=list(label="first", sample="sum", "ignore"))
        if(j == length(V(g_directed)))
          g_directed <- contract(g_directed, c(1:(length(V(g_directed))-1), i), vertex.attr.comb=list(label="first", sample="sum", "ignore"))
        g_directed <- simplify(g_directed, edge.attr.comb=list(capacity="sum", weight="sum", "ignore"))
        # Combining vertices i and j into one (namely i)
        partitionvec[partitionvec == j] <- -42
        partitionvec[partitionvec > j] <- partitionvec[partitionvec > j] - 1
        partitionvec[partitionvec == -42] <- ifelse(i > j, i-1, i)
        if(i > j) i <- i-1
#        print(V(g_directed))
        # Adding a reference to every index referring to j that these nodes now are combined into i.
        #print(partitionvec)
#        print(length(V(g_directed))-1)
#      }
    }
  }
  return(list(ncut = ij_min_ncut, ncut_min = ij_min_n))
  # Returns a list consisting of the NCut value and the partition attaining it
}
# Realizes the NCut-via-st-MinCuts algorithm. Input: Extended graph object returned from generate_graph
# (or generate_cockroach or generate_r_graph). Output: List of XCut values with corresponding partitions.
# This function is a vastly improved version of xcut_via_stMinCuts for NCut only
# NOTE: THIS CURRENTLY DOES NOT WORK (AND PROBABLY NEVER WILL). DO NOT USE!

xcut_exact_OLD <- function(g){
  
  mcut_min <- Inf
  rcut_min <- Inf
  ncut_min <- Inf
  ccut_min <- Inf
  # Keeps track of the minimum respective XCut. Again, apologies for ugliness. Should be efficient, though.
  
  nvert <- length(V(g$graph))
  if(nvert > 30){
    print("Error: More than 30 nodes.")
    return(NULL)
  }
  # Hard-coded check whether the number of nodes isn't too large in order to keep computations manageable.
  
  for(k in 1:(2^(nvert-1)-1)){
    
    #if(k%%5000==0) print(paste0(k,"/",2^(nvert-1)))
    # If one wants some indication as to how long the algorithm will take.
#    if(length(g$removed_vertices) > 0)
#      p <- ((1:dim(g$weight_mat)[1])[-g$removed_vertices])[intToBits(k)[1:nvert]==1]
#    else
      p <- (1:dim(g$weight_mat)[1])[intToBits(k)[1:nvert]==1]
    # Converts a number k into a partition p by binary encoding.
    
    xcut_p <- xcut_partition_S(p, g)
    # Computes the XCut of partition p.
#    print(k)
#    print(p)
#    print(xcut_p)
    
    if(xcut_p$mcut_min < mcut_min){
      mcut_min <- xcut_p$mcut_min
      mcut <- p
    }
    if(xcut_p$rcut_min < rcut_min){
      rcut_min <- xcut_p$rcut_min
      rcut <- p
    }
    if(xcut_p$ncut_min < ncut_min){
      ncut_min <- xcut_p$ncut_min
      ncut <- p
    }
    if(xcut_p$ccut_min < ccut_min){
      ccut_min <- xcut_p$ccut_min
      ccut <- p
    }
    # Updates the optimal XCuts and the partition attaining them.
  }
  return(list(mcut=mcut, rcut=rcut, ncut=ncut, ccut=ccut,
              mcut_min=mcut_min, rcut_min=rcut_min, ncut_min=ncut_min, ccut_min=ccut_min))
}
# Computes the exact XCut value of a graph by brute force. Only works for small graphs (i.e. with <25 nodes).
# Input: Extended graph object g. Output: List of XCut values and partitions attaining them (includes MinCut).

plot_graph_partition_OLD <- function(g, partition, vertex_scale=30, axes=FALSE, edge_scale=20,
                                 vertex_offset=30, edge_offset=1, xlim=c(1,M), ylim=c(1,M), main=NULL){
  #V(g$graph)[unlist(partition)]$color <- 'blue'
  #V(g$graph)[-unlist(partition)]$color <- 'green3'
  
  if(length(partition) < length(V(g$graph))){
    V(g$graph)$color <- 'green3'
    #V(g$graph)[which(as.numeric(V(g$graph)$label) %in% partition)]$color <- 'blue'
    V(g$graph)[unlist(partition)]$color <- 'blue'
  } else if(length(partition) == length(V(g$graph))){
    allcolors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    V(g$graph)$color <- c("green3", "blue", sample(allcolors, length(unique(partition))-2))[partition]
    #colorRampPalette(c("red", "green3"))
  }
  edge_weights <- E(g$graph)$weight
  vertex_weights <- V(g$graph)$sample
  #vertex_weights <- V(g$graph)$smmm
  par(cex=0.7, mai=c(0,0.1,0,0.1))
  #plot(g$graph, edge.width = E(g$graph)$weight, vertex.size=vertex_size, axes=TRUE)
  if(axes){
    plot(g$graph, rescale=FALSE, axes=TRUE, xlim=xlim, ylim=ylim, asp = 0, vertex.label=NA,
         vertex.size = vertex_size, main=main)
  }else {
    plot(g$graph, vertex.label=NA, edge.width = edge_scale*edge_weights/max(edge_weights)+edge_offset,
         #vertex.size = vertex_scale*vertex_weights/1.5+15, main=main)
         vertex.size = vertex_scale*vertex_weights/max(vertex_weights)+vertex_offset, main=main)
  }
}
# Plots a partition 'partition' of the graph g$graph.
# Input: extended graph object g, several graphical parameters. Output: None.

generate_cockroach <- function(k){

  cr_coords <- rbind(cbind(1:(2*k-1), 2:(2*k)), cbind((2*k+1):(4*k-1), (2*k+2):(4*k)), cbind((k+1):(2*k), (3*k+1):(4*k)))
  # Gives a matrix containing the coordinates of the cockroach. rbind binds coordinates from both rows of the cockroach.
  
  cr_capacity <- rep(1, length(cr_coords[,1]))
  # Every edge of the cockroach graph has weight 1.
  
  cr_weight_mat <- matrix(rep(0, (4*k)^2), 4*k, 4*k)
  cr_weight_mat[cr_coords] <- 1
  cr_weight_mat[lower.tri(cr_weight_mat)] <- t(cr_weight_mat)[lower.tri(cr_weight_mat)]
  # Define weight matrix of the cockroach.
  
  g <- graph(edges=as.vector(t(cr_coords)), n=4*k, directed = FALSE)
  # Generate the graph from the edges 
  
  V(g)$x <- c(seq(1,2*k), seq(1,2*k))
  V(g)$y <- c(rep(1,2*k), rep(2,2*k))
  V(g)$label = 1:(4*k)
  # Assign coordinates and labels for plotting
  
  E(g)$capacity <- rep(1, length(E(g)))
  E(g)$weight <- 100*cr_capacity/sum(cr_capacity)
  # E(g)$capacity is only for using igraph::min_cut, $weight is for plotting purposes.
  
  return(list(graph=g, weight_mat=cr_weight_mat, removed_vertices=NULL))
  # Return value analogous to that of generate_graph.
}
# Generates a cockroach graph, see Guattery, Miller (1993). The integer k determines the length of the graph.

xcut_partition_S_OLD <- function(partition, g, with_removed=FALSE){
  
  if(with_removed)
    S <- partition - vapply(partition, function(x)sum(x > removed_vertices), 1)
  else
    S <- partition
  S_C <- (1:length(V(g$graph)))[-S]
  #S <- intersect(partition, as.numeric(V(g$graph)$label))
  #S_C <- intersect((1:dim(g$weight_mat)[1])[-partition], as.numeric(V(g$graph)$label))
#  print(S)
#  print(S_C)
  # Selecting the correct vertices of g$graph to determine the partition S and its complement S_C.
  # If with_removed == TRUE, the indices of the partition include removed_vertices in their counting.
  
  S_vol <- sum(g$weight_mat[S,])
  S_C_vol <- sum(g$weight_mat[S_C,])
  # Computing the volumes of S and S_C.
  
  mcut_min <- sum(g$weight_mat[S, S_C])
  rcut_min <- sum(g$weight_mat[S, S_C]) / (length(S) * length(S_C))
  ncut_min <- sum(g$weight_mat[S, S_C]) / (S_vol * S_C_vol)
  ccut_min <- sum(g$weight_mat[S, S_C]) / min(S_vol, S_C_vol)
  # Computes the various XCuts of S.
  
  return(list(mcut_min=mcut_min, rcut_min=rcut_min, ncut_min=ncut_min, ccut_min=ccut_min))
  # Returns a list of the computed XCut values of S.
}
# Computes the XCut values of partition 'partition' w.r.t. g$graph and its weight matrix.
# Input: Partition 'partition', i.e. a subset of the powerset of 1:length(V(g$graph)), and extended graph object g.
# with_removed determines whether partition includes g$removed_vertices in its counting.

xcut_cov <- function(p1, partitionS, gS, partitionR, gR){

  p <- as.vector(t(matrix(p1, sqrt(length(p1)))))
  # Re-ordering the probability vector to match it with the nodes from the graph.
  
  S <- intersect(partitionS, as.numeric(V(gS$graph)$label))
  S_C <- intersect((1:dim(gS$weight_mat)[1])[-partitionS], as.numeric(V(gS$graph)$label))
  # Computes the complement of S
  
  q_S <- sapply(1:length(p), function(x)sum(p[S[gS$weight_mat[S,x]>0 | gS$weight_mat[x,S]>0]]))
  q_S_C <- sapply(1:length(p), function(x)sum(p[S_C[gS$weight_mat[S_C,x]>0 | gS$weight_mat[x,S_C]>0]]))
  vol_S <- sum(p*q_S)
  vol_S_C <- sum(p*q_S_C)
  # Computes q and the volume of S and S_C. q_S from the thesis corresponds to q_S[S_C] (and q_S^C to q_S_C[S]).
  # q_S[S] and q_S_C[S_C] also serve as substitutes for b, and q_S + q_S_C for a (also from the thesis).
  
  if(vol_S < vol_S_C){
    vol_S_min <- vol_S
    S_min <- S
    S_C_min <- S_C
    q_S_min <- q_S
    q_S_C_min <- q_S_C
  }
  else{
    vol_S_min <- vol_S_C
    S_min <- S_C
    S_C_min <- S
    q_S_min <- q_S_C
    q_S_C_min <- q_S
  }
  # Determines the partition (either S or S_C) with minimal volume for Cheeger Cut.
  # NOTE: THIS WILL NOT WORK IF THE VOLUMES OF S AND S_C ARE EQUAL! THIS CASE IS COVERED BY cov_ccut_equal_vol.
  
  g_xc_S <- xcut_partition_S(S, gS)
  g_nc_S <- g_xc_S$ncut
  g_cc_S <- g_xc_S$ccut
  # Computes the NCut and Cheeger Cut values of S.
  
  R <- intersect(partitionR, as.numeric(V(gR$graph)$label))
  R_C <- intersect((1:dim(gR$weight_mat)[1])[-partitionR], as.numeric(V(gR$graph)$label))
  q_R <- sapply(1:length(p), function(x)sum(p[R[gR$weight_mat[R,x]>0 | gR$weight_mat[x,R]>0]]))
  q_R_C <- sapply(1:length(p), function(x)sum(p[R_C[gR$weight_mat[R_C,x]>0 | gR$weight_mat[x,R_C]>0]]))
  vol_R <- sum(p*q_R)
  vol_R_C <- sum(p*q_R_C)
  if(vol_R<=vol_R_C){
    vol_R_min <- vol_R
    R_min <- R
    R_C_min <- R_C
    q_R_min <- q_R
    q_R_C_min <- q_R_C
  }
  else{
    vol_R_min <- vol_R_C
    R_min <- R_C
    R_C_min <- R
    q_R_min <- q_R_C
    q_R_C_min <- q_R
  }
  g_xc_R <- xcut_partition_S(R, gR)
  g_nc_R <- g_xc_R$ncut
  g_cc_R <- g_xc_R$ccut
  # Repeat the same computations for the second partition R.
  
  SnR <- intersect(S, R)
  S_CnR <- intersect(S_C, R)
  SnR_C <- intersect(S, R_C)
  S_CnR_C <- intersect(S_C, R_C)
  S_minnR_min <- intersect(S_min, R_min)
  S_C_minnR_min <- intersect(S_C_min, R_min)
  S_minnR_C_min <- intersect(S_min, R_C_min)
  S_C_minnR_C_min <- intersect(S_C_min, R_C_min)
  # These intersections are necessary to account for the partitions each node belongs to.
  
  #print(sum(p[SnR]*q_S_C[SnR]*q_R_C[SnR]))
  #print(sum(p[SnR_C]*q_S_C[SnR_C]*q_R[SnR_C]))
  #print(sum(p[S_CnR]*q_S[S_CnR]*q_R_C[S_CnR]))
  #print(sum(p[S_CnR_C]*q_S[S_CnR_C]*q_R[S_CnR_C]))
  #print(p[SnR_C]*q_S_C[SnR_C]*q_R[SnR_C])
  
  return(list(
    mcut_cov = sum(p[SnR]*q_S_C[SnR]*q_R_C[SnR]) + sum(p[SnR_C]*q_S_C[SnR_C]*q_R[SnR_C]) + sum(p[S_CnR]*q_S[S_CnR]*q_R_C[S_CnR]) + sum(p[S_CnR_C]*q_S[S_CnR_C]*q_R[S_CnR_C]) - 4*g_xc_S$mcut*g_xc_R$mcut,
    rcut_cov = (sum(p[SnR]*q_S_C[SnR]*q_R_C[SnR]) + sum(p[SnR_C]*q_S_C[SnR_C]*q_R[SnR_C]) + sum(p[S_CnR]*q_S[S_CnR]*q_R_C[S_CnR]) + sum(p[S_CnR_C]*q_S[S_CnR_C]*q_R[S_CnR_C]))/(length(S)*length(S_C)*length(R)*length(R_C)) - 4*g_xc_S$rcut*g_xc_R$rcut,
    ncut_cov = (sum(p[SnR] * (q_S_C[SnR] - g_nc_S*(q_S_C[SnR]*(vol_S - vol_S_C) + 2*vol_S_C*(q_S[SnR] + q_S_C[SnR])))
                    * (q_R_C[SnR] - g_nc_R*(q_R_C[SnR]*(vol_R - vol_R_C) + 2*vol_R_C*(q_R[SnR] + q_R_C[SnR]))))
                + sum(p[SnR_C] * (q_S_C[SnR_C] - g_nc_S*(q_S_C[SnR_C]*(vol_S - vol_S_C) + 2*vol_S_C*(q_S[SnR_C] + q_S_C[SnR_C])))
                      * (q_R[SnR_C] - g_nc_R*(q_R[SnR_C]*(vol_R_C - vol_R) + 2*vol_R*(q_R[SnR_C] + q_R_C[SnR_C]))))
                + sum(p[S_CnR] * (q_S[S_CnR] - g_nc_S*(q_S[S_CnR]*(vol_S_C - vol_S) + 2*vol_S*(q_S[S_CnR] + q_S_C[S_CnR])))
                      * (q_R_C[S_CnR] - g_nc_R*(q_R_C[S_CnR]*(vol_R - vol_R_C) + 2*vol_R_C*(q_R[S_CnR] + q_R_C[S_CnR]))))
                + sum(p[S_CnR_C] * (q_S[S_CnR_C] - g_nc_S*(q_S[S_CnR_C]*(vol_S_C - vol_S) + 2*vol_S*(q_S[S_CnR_C] + q_S_C[S_CnR_C])))
                      * (q_R[S_CnR_C] - g_nc_R*(q_R[S_CnR_C]*(vol_R_C - vol_R) + 2*vol_R*(q_R[S_CnR_C] + q_R_C[S_CnR_C])))))
    / (vol_S * vol_S_C * vol_R * vol_R_C) - 4*g_nc_S*g_nc_R,
    ccut_cov = (sum(p[S_minnR_min]*(q_S_C_min[S_minnR_min] - g_cc_S*(q_S_C_min[S_minnR_min]+2*q_S_min[S_minnR_min]))*(q_R_C_min[S_minnR_min]-g_cc_R*(q_R_C_min[S_minnR_min]+2*q_R_min[S_minnR_min])))
                + sum(p[S_minnR_C_min]*(q_S_C_min[S_minnR_C_min] - g_cc_S*(q_S_C_min[S_minnR_C_min]+2*q_S_min[S_minnR_C_min]))*(q_R_min[S_minnR_C_min]-g_cc_R*(q_R_C_min[S_minnR_C_min])))
                + sum(p[S_C_minnR_min]*(q_S_min[S_C_minnR_min] - g_cc_S*(q_S_C_min[S_C_minnR_min]))*(q_R_C_min[S_C_minnR_min]-g_cc_R*(q_R_C_min[S_C_minnR_min]+2*q_R_min[S_C_minnR_min])))
                + sum(p[S_C_minnR_C_min]*(q_S_min[S_C_minnR_C_min] - g_cc_S*(q_S_C_min[S_C_minnR_C_min]))*(q_R_min[S_C_minnR_C_min]-g_cc_R*(q_R_C_min[S_C_minnR_C_min]))))/(vol_S_min*vol_R_min)))
  # I apologize. This is correct, though not pleasant to look at. It is the case-sensitive implementation of
  # the covariance between the asymptotic distributions of sqrt(n) * (\hat{XC(S)} - XC(S)) and the same for R.
}
# Computes the covariance between the asymptotic distributions of sqrt(n) * (\hat{XC(S)} - XC(S)) and
# sqrt(n) * (\hat{XC(R)} - XC(R)) for partitions S and R of the graphs gS$graph and gR$graph.
# gS and gR are extended graph objects built from a sample from the distribution defined by the probability vector p.
# Input: Probability vector p underlying the empirical graphs gS$graph and gR$graph, partitions S and R, and
# extended graph objects gS and gR. Output: A list of values of the covariances mentioned above for XCut.
# NOTE: THIS WILL RETURN THE WRONG CHEEGER CUT COVARIANCE IF VOL(S)==VOL(S^C) OR VOL(R)==VOL(R^C).
# IN THIS CASE, USE ccut_cov_equal_vol.

ccut_cov_equal_vol <- function(p1, S, gS, R, gR){
  
  p <- as.vector(t(matrix(p1, sqrt(length(p1)))))
  # Re-ordering the probability vector to match it with the nodes from the graph.
  
  S_C <- (1:length(p))[-S]
  q_S <- sapply(1:length(p), function(x)sum(p[S[gS$weight_mat[S,x]>0 | gS$weight_mat[x,S]>0]]))
  q_S_C <- sapply(1:length(p), function(x)sum(p[S_C[gS$weight_mat[S_C,x]>0 | gS$weight_mat[x,S_C]>0]]))
  vol_S <- sum(p*q_S)
  vol_S_C <- sum(p*q_S_C)
  vol_S_C <- sum(p*q_S_C)
  # Computing the various statistics around p and gS$graph with associated weight matrix gS$weight_mat.
  
  Sigma_mult <- diag(p) - outer(p,p)
  # The covariance matrix Sigma of a multinomial distribution with probability vector p.
  
  Z_sample_ccut_cov <- function(){
  # We need to define a function that generates a multinomial samples Z and then returns a sample of
  # the limiting covariance of the CCut statistic.
      
    Z <- as.vector(mvtnorm::rmvnorm(1, sigma=Sigma_mult))
    # Generating a sample Z following a multivariate normal distribution (see  thesis, Thm. 3.12).
    
    aZ_S <- sapply(1:length(p), function(x)sum(Z[gS$weight_mat[,x]>0 | gS$weight_mat[x,]>0]))
    # Compute the statistic 'a' w.r.t. Z (see Thm. 3.12).
    
    minequalS <- sum(Z[S]*(q_S[S]+q_S_C[S])) + sum(p[S]*aZ_S[S])
    minequalS_C <- sum(Z[S_C]*(q_S[S_C]+q_S_C[S_C])) + sum(p[S_C]*aZ_S[S_C])
    # Compute the statistics to check whether S or S_C is the one stored as \tilde{S}.
    
    if(vol_S < vol_S_C | (vol_S == vol_S_C & minequalS <= minequalS_C)){
      vol_S_min <- vol_S
      S_min <- S
      S_C_min <- S_C
      q_S_min <- q_S
      q_S_C_min <- q_S_C
    }
    else{
      vol_S_min <- vol_S_C
      S_min <- S_C
      S_C_min <- S
      q_S_min <- q_S_C
      q_S_C_min <- q_S
    }
    # Determines \tilde{S}.
    
    g_xc_S <- xcut_partition_S(S, gS)
    g_nc_S <- g_xc_S$ncut
    g_cc_S <- g_xc_S$ccut
    # Computing the XCut.
    
    R_C <- (1:length(p))[-R]
    q_R <- sapply(1:length(p), function(x)sum(p[R[gR$weight_mat[R,x]>0 | gR$weight_mat[x,R]>0]]))
    q_R_C <- sapply(1:length(p), function(x)sum(p[R_C[gR$weight_mat[R_C,x]>0 | gR$weight_mat[x,R_C]>0]]))
    vol_R <- sum(p*q_R)
    vol_R_C <- sum(p*q_R_C)
    aZ_R <- sapply(1:length(p), function(x)sum(Z[gR$weight_mat[,x]>0 | gR$weight_mat[x,]>0]))
    minequalR <- sum(Z[R]*(q_R[R]+q_R_C[R])) + sum(p[R]*aZ_R[R])
    minequalR_C <- sum(Z[R_C]*(q_R[R_C]+q_R_C[R_C])) + sum(p[R_C]*aZ_R[R_C])
    if(vol_R < vol_R_C | (vol_R == vol_R_C & minequalR <= minequalR_C)){
      vol_R_min <- vol_R
      R_min <- R
      R_C_min <- R_C
      q_R_min <- q_R
      q_R_C_min <- q_R_C
    }
    else{
      vol_R_min <- vol_R_C
      R_min <- R_C
      R_C_min <- R
      q_R_min <- q_R_C
      q_R_C_min <- q_R
    }
    g_cc_R <- xcut_partition_S(R, gR)$ccut
    # Doing the same for R.
    
    xcut_cov_sum_component <- function(x, aTbF){
      if(aTbF & x%in%S_min) return((q_S_C_min[x] - g_cc_S*(q_S_C_min[x]+2*q_S_min[x]))/vol_S_min)
      if(aTbF & x%in%S_C_min) return((q_S_min[x] - g_cc_S*(q_S_C_min[x]))/vol_S_min)
      if(!aTbF & x%in%R_min) return((q_R_C_min[x] - g_cc_R*(q_R_C_min[x]+2*q_R_min[x]))/vol_R_min)
      if(!aTbF & x%in%R_C_min) return((q_R_min[x] - g_cc_R*(q_R_C_min[x]))/vol_R_min)
    }
    # Defining the components of the sum to compute the asymptotic distribution of the XCut statistic.
    
    return(list(ESR = outer(1:length(p), 1:length(p), Vectorize(function(i,j)
              xcut_cov_sum_component(i, TRUE) * xcut_cov_sum_component(j, FALSE) * Z[i] * Z[j])),
            ES = sum(Z*sapply(1:length(p), function(i)xcut_cov_sum_component(i, TRUE))),
            ER = sum(Z*sapply(1:length(p), function(j)xcut_cov_sum_component(j, FALSE)))))
    # ESR is the "empirical expectation" of the asymptotic statistic of S times that of R.
    # ES (resp. ER) is the "empirical expectation" of the asymptotic statistic of S (resp. R).
  }
  
  Z_cov_sample <- replicate(1000, Z_sample_ccut_cov())
  return(mean(Z_cov_sample[1,]) - mean(Z_cov_sample[2,] * mean(Z_cov_sample[3,])))
  # Returns the (estimated) asymptotic covariance of the Cheeger Cut statistic.
}
# Estimates the asymptotic covariance of the Cheeger Cut statistics sqrt(n) * (\hat{XC(S)} - XC(S)) and
# sqrt(n) * (\hat{XC(R)} - XC(R)) under consideration that vol(S)==vol(S^C) or vol(R)==vol(R^C) may be possible.
# Input and Output: See xcut_cov.

xvstlocmax_fast <- function(g, allcuts=TRUE, progress_bar=TRUE, mergeByDegree=TRUE){
  
  ij_min_r <- Inf
  ij_min_n <- Inf
  ij_min_c <- Inf
  ij_min_rcut <- NULL
  ij_min_ncut <- NULL
  ij_min_ccut <- NULL
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  g_directed <- as.directed(g$graph, mode="mutual")
  # The st-MinCut function igraph::stMincuts only accepts directed graphs.
  
  #locmaxvec <- sort(which(V(g$graph)$local_maximum))
  #locmaxvec <- order(V(g$graph)$local_maximum, decreasing=TRUE)
  #locmaxvec <- locmaxvec[V(g$graph)$local_maximum[locmaxvec] > 0]
  locmaxvec_idxs <- which(V(g$graph)$local_maximum)
  if(mergeByDegree) locmaxvec <- locmaxvec_idxs[order(V(g$graph)$degree[locmaxvec_idxs], decreasing=TRUE)]
  else{
    locmaxdeg_idx <- which.max(V(g$graph)$degree[locmaxvec_idxs])
    #distvec <- c(distances(g$graph, v=locmaxvec_idxs[locmaxdeg_idx], to=locmaxvec_idxs))
    #locmaxvec <- locmaxvec_idxs[order(distvec, decreasing=FALSE)]
    locmaxvec <- locmaxvec_idxs[order(V(g$graph)$sample[locmaxvec_idxs], decreasing=TRUE)]
  }
  # Computes the set of local maxima of g and sorts it by degree.
  print(locmaxvec)
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = length(locmaxvec)-1, style = 3)
  # Setting up the progress bar.
  
  if(length(locmaxvec) < 2) return(xcut_via_stmincut_OLD(g))
  
  for(i in 2:length(locmaxvec)){
#    print(paste(locmaxvec[1],"|",locmaxvec[i]))
    if(allcuts) ij_cut <- stMincuts(g_directed, locmaxvec[1], locmaxvec[i])
    # Computes all partitions attaining the ij-MinCut.
    
    if(length(ij_cut$partition1s) == 0){
      # stMincuts doesn't return a partition if the MinCut value is zero.
      
      ij_cut <- min_cut(g$graph, source=locmaxvec[1], target=locmaxvec[i], value.only=FALSE)
      # In that case we use igraph::min_cut to determine the st-MinCut partition.
      
      ij_cut$partition1s <- list(ij_cut$partition1)
      # Stores the partition as a one-element list for the for loop below.
      
    }
    
    for(S_ij in ij_cut$partition1s){
      # Goes through all partitions attaining the ij-MinCut.
      
#      shift_a_back_by_b <- function(a,b) a + sapply(a, function(i)sum(b - 1:length(b) <= i))
      # Shifts a list (of vertices) a back to their original position in 1:m assuming
      # they were a subset of the vector obtained by removing b from 1:m.
      
      #print(length(S_ij))
      #print(locmaxvec[1:(i-1)])
      if(locmaxvec[1] %in% S_ij) S_ij <- base::union(S_ij, locmaxvec[1:(i-1)])#c(S_ij, locmaxvec[2:(i-1)])
      else S_ij <- base::setdiff(S_ij, locmaxvec[1:(i-1)])
#      else if(i > 2) S_ij <- S_ij#shift_a_back_by_b(S_ij, locmaxvec[2:(i-1)])
#      # Recover the original partition (before any merging, e.g. if i==2) of S_ij by shifting it back.
#      #print(locmaxvec[2:(i-1)])
      #print(S_ij)
      
      ij_xcut <- xcut_partition_S_OLD(S_ij, g)
      # Computes the XCut of the partition S_ij.
      
      #print(paste(locmaxvec[1],"|",locmaxvec[i]))
      #print(S_ij)
      #print(ij_xcut$ncut_min)
      #print("----")
      
      if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
        ij_min_n <- ij_xcut$ncut_min
        ij_min_ncut <- S_ij
#        print(paste(locmaxvec[1],"|",locmaxvec[i]))
#        print(S_ij)
      }
      if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
        ij_min_c <- ij_xcut$ccut_min
        ij_min_ccut <- S_ij
      }
      if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
        ij_min_r <- ij_xcut$rcut_min
        ij_min_rcut <- S_ij
      }
      # Updates the best XCut and the partition that attains it.
    }
    
    # We now merge vertices locmaxvec[1] and locmaxvec[i] by combining them into one (namely vertex locmaxvec[1]).
    if(i < length(locmaxvec)){
      merged_indices <- 1:length(V(g_directed))
      #print(locmaxvec)
      locmaxvec[c(1,i)] <- c(max(locmaxvec[i],locmaxvec[1]), min(locmaxvec[i],locmaxvec[1]))
      #print(locmaxvec)
      merged_indices[locmaxvec[i]] <- locmaxvec[1]
      #print(i)
      #print(merged_indices)
      #print(length(V(g_directed)))
      #print(V(g_directed)$degree)
      
      g_directed <- contract(g_directed, merged_indices, vertex.attr.comb=list(degree="sum","last"))
#    E(g_directed)[E(g_directed) == locmaxvec[i]] <- locmaxvec[1]
#    # Redirects all edges involving locmaxvec[i] to locmaxvec[1].
      #print(V(g_directed)$degree)
    
      g_directed <- simplify(g_directed, edge.attr.comb=list(capacity="sum", weight="sum", width="sum"))
      # Collapses multiple edges into one by summing their weights etc., and deletes loops.
    }
    
#    
#    g_directed <- delete.vertices(g_directed, locmaxvec[i])
#    # Deletes locmaxvec[i] to complete the merging process.
    
    if(progress_bar) setTxtProgressBar(progress, i-1)
    # Updates progress bar.
    
    #print(i)
    #print(j)
    #print(S_ij)
    #print(ij_xcut$ncut_min)
    
  }
  if(progress_bar) close(progress)
  #    bla <- c()
  #    for(i in 1:(length(V(g$graph))-1)){
  #      for(j in (i+1):length(V(g$graph))){
  #        ij_cut <- stMincuts(g_directed, i, j)
  #        # Computes all partitions attaining the ij-MinCut.
  #        
  # #       for(S_ij in ij_cut$partition1s){
  #    ij_xcut <- xcut_partition_S(S_ij, g)
  #    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min == ij_min_n){
  #      bla <- c(bla, i, j)
  #  }#}}
  #}
  return(list(rcut = ij_min_rcut, rcut_min = ij_min_r,
              ncut = ij_min_ncut, ncut_min = ij_min_n,
              ccut = ij_min_ccut, ccut_min = ij_min_c))#, bla=bla))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}
# Realizes the XCut-via-st-MinCuts algorithm, restricted to pairs of "local maxima" (i.e. vertices with more observations than all their neigbours).
# Input: Extended graph object returned from generate_graph (or generate_cockroach or generate_r_graph).
# Output: List of XCut values with corresponding partitions.

xvstlocmax_merge <- function(g, allcuts=TRUE, progress_bar=TRUE){
  
  ij_min_r <- Inf
  ij_min_n <- Inf
  ij_min_c <- Inf
  ij_min_rcut <- NULL
  ij_min_ncut <- NULL
  ij_min_ccut <- NULL
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  g_directed <- as.directed(g$graph, mode="mutual")
  # The st-MinCut function igraph::stMincuts only accepts directed graphs.
  
  locmaxvec <- which(V(g$graph)$local_maximum)
  locmaxdeg_idx <- which.max(V(g$graph)$degree[locmaxvec])
  locmaxvec <- c(locmaxvec[locmaxdeg_idx], locmaxvec[-locmaxdeg_idx])
  # Computes the local maxima of g and the local maximum with highest degree, putting the latter at position 1.
  
  already_merged <- c(locmaxvec[1])
  # Vector to keep track of vertices that have already merged with locmaxvec[1].
  
  distmat <- distances(g$graph, weights=E(g$graph)$distance, to=locmaxvec, mode="all")
  
  disconnected_locmax <- which(is.infinite(distmat[locmaxvec[1],]))
  if(length(disconnected_locmax) > 0){
    distmat[locmaxvec[1],disconnected_locmax] <- 0
    #warning("Graph is disconnected. All cut values will be zero.")
    locmaxvec <- locmaxvec[c(1,disconnected_locmax[1])]
    distmat <- distmat[,c(1,disconnected_locmax[1])]
    # If the graph is disconnected, perform only one MinCut
    # (between locmaxvec[1] and locmaxvec[disconnected_locmax[1]]), ensuring a cut value of 0.
    #return(xvstlocmax_fast(g=g, allcuts=allcuts, progress_bar=progress_bar))
  }
  distmat[cbind(locmaxvec,1:length(locmaxvec))] <- Inf
  # Computing the matrix of shortest path distances between every vertex and locmaxvec, assigning Inf to "loops".
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = length(locmaxvec)-1, style = 3)
  # Setting up the progress bar.
  
  
  while(length(locmaxvec) > 1){
    #    print(paste(locmaxvec[1],"|",locmaxvec[i]))
#    cat("distmat:",distmat[locmaxvec[1],],"\n")
    
    #mindist_idx <- which.min(distances(g_directed, v=locmaxvec[1], to=locmaxvec[-1])) + 1
    mindist_idx <- which.min(distmat[locmaxvec[1],])
    # Determines the locmax index of the vertex attaining the minimal distance to locmaxvec[1].
    
#    cat("locmaxvec:",locmaxvec,"\n")
#    cat("mindist_idx:",mindist_idx,"\n")
    
    if(allcuts) ij_cut <- stMincuts(g_directed, locmaxvec[1], locmaxvec[mindist_idx])
    # Computes all partitions attaining the ij-MinCut.
    
    if(length(ij_cut$partition1s) == 0){
      # stMincuts doesn't return a partition if the MinCut value is zero.
      
      ij_cut <- min_cut(g$graph, source=locmaxvec[1], target=locmaxvec[mindist_idx], value.only=FALSE)
      # In that case we use igraph::min_cut to determine the st-MinCut partition.
      
      ij_cut$partition1s <- list(ij_cut$partition1)
      # Stores the partition as a one-element list for the for loop below.
      
    }
    
    for(S_ij in ij_cut$partition1s){
      # Goes through all partitions attaining the ij-MinCut.
      
      #print(S_ij)
      if(locmaxvec[1] %in% S_ij) S_ij <- union(S_ij, already_merged)#c(S_ij, locmaxvec[2:(i-1)])
      else S_ij <- setdiff(S_ij, already_merged)
      #print(S_ij)
      
      ij_xcut <- xcut_partition_S(S_ij, g)
      # Computes the XCut of the partition S_ij.
      
      if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
        ij_min_n <- ij_xcut$ncut_min
        ij_min_ncut <- S_ij
        #        print(paste(locmaxvec[1],"|",locmaxvec[i]))
        #        print(S_ij)
      }
      if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
        ij_min_c <- ij_xcut$ccut_min
        ij_min_ccut <- S_ij
      }
      if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
        ij_min_r <- ij_xcut$rcut_min
        ij_min_rcut <- S_ij
      }
      # Updates the best XCut and the partition that attains it.
    }
    
    
    if(length(locmaxvec) > 2){
    # We now merge vertices locmaxvec[1] and locmaxvec[i] by combining them into one (namely vertex locmaxvec[1]).
  
      vertices_to_merge <- shortest_paths(g_directed, from=locmaxvec[1], to=locmaxvec[mindist_idx],
                                          weights=E(g_directed)$distance, mode="all", output="vpath")$vpath[[1]]
      # Determines the vertices constituting the shortest path between locmaxvec[1] and locmaxvec[mindist_idx] for merging.
#      cat("vertices_to_merge:",vertices_to_merge,"\n")
      
      mergelist <- (1:length(V(g_directed))) %>% (function(a){a[vertices_to_merge] <- locmaxvec[1];a})
  #    cat("mergelist:",mergelist,"\n")
      g_directed <- contract(g_directed, mergelist, vertex.attr.comb=list(degree="sum","ignore"))
      # Contracts all vertices in vertices_to_merge into one (namely locmaxvec[1]).
      
      g_directed <- simplify(g_directed, edge.attr.comb=list(capacity="sum", weight="sum", width="sum", distance="min"))
      # Collapses all resulting multiple edges into one by summing their weights etc., and deletes loops, thus finishing the merging.
      
      distmat[locmaxvec[1],2:length(locmaxvec)] <- sapply(2:length(locmaxvec), function(i)min(distmat[vertices_to_merge,i]))
      # Updates distmat to account for vertices_to_merge being merged into locmaxvec[1] (=vertices_to_merge[1])
      
      already_merged <- c(already_merged, vertices_to_merge[-1])
      # Updates the already merged vertices.
    }
    locmaxvec <- locmaxvec[-mindist_idx]
    distmat <- distmat[,-mindist_idx]
    # Prevents locmaxvec[mindist_idx] from being merged again.
    
    if(progress_bar) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
    # Updates progress bar.
    
    #print(i)
    #print(j)
    #print(S_ij)
    #print(ij_xcut$ncut_min)
    
  }
  if(progress_bar) close(progress)
  #    bla <- c()
  #    for(i in 1:(length(V(g$graph))-1)){
  #      for(j in (i+1):length(V(g$graph))){
  #        ij_cut <- stMincuts(g_directed, i, j)
  #        # Computes all partitions attaining the ij-MinCut.
  #        
  # #       for(S_ij in ij_cut$partition1s){
  #    ij_xcut <- xcut_partition_S(S_ij, g)
  #    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min == ij_min_n){
  #      bla <- c(bla, i, j)
  #  }#}}
  #}
  return(list(rcut = ij_min_rcut, rcut_min = ij_min_r,
              ncut = ij_min_ncut, ncut_min = ij_min_n,
              ccut = ij_min_ccut, ccut_min = ij_min_c))#, bla=bla))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}
# Realizes the XCut-via-st-MinCuts algorithm, restricted to pairs of "local maxima" (i.e. vertices with more observations than all their neigbours).
# Input: Extended graph object returned from generate_graph (or generate_cockroach or generate_r_graph).
# Output: List of XCut values with corresponding partitions.

xist_OLD <- function(g, allcuts=TRUE, progress_bar=TRUE){
  
  ij_min_r <- Inf
  ij_min_n <- Inf
  ij_min_c <- Inf
  ij_min_rcut <- NULL
  ij_min_ncut <- NULL
  ij_min_ccut <- NULL
  # I apologize for this ugliness. This is where the respective Xvst minimum and associated partition are stored.
  
  g_directed <- as.directed(g$graph, mode="mutual")
  # The st-MinCut function igraph::stMincuts only accepts directed graphs.
  V(g_directed)$tempnr <- 1:length(V(g_directed))
  
  print(dim(g$weight_mat))
  print(length(V(g$graph)))
  
  locmaxvec <- which(sapply(1:length(V(g$graph)), function(i)all(g$weight_mat[i,]==0 | V(g$graph)$sample[i] >= V(g$graph)$sample)))
  #locmaxvec <- which(V(g$graph)$local_maximum)
  locmaxdeg_idx <- which.max(V(g$graph)$sample[locmaxvec])
  locmaxvec <- c(locmaxvec[locmaxdeg_idx], locmaxvec[-locmaxdeg_idx])
  # Computes the local maxima of g and the local maximum with highest degree, putting the latter at position 1.
  
#  already_merged <- c(locmaxvec[1])
  # Vector to keep track of vertices that have already merged with locmaxvec[1].
  
  #distmat <- distances(g$graph, weights=E(g$graph)$distance, v=locmaxvec, to=locmaxvec, mode="all")
  distmat <- distances(g$graph, weights=E(g$graph)$capacity, v=locmaxvec, to=locmaxvec, mode="all")
  #distmat <- t(sapply(1:length(locmaxvec), function(i)order(V(g$graph)$degree[locmaxvec], decreasing=TRUE)))
  #distmat <- matrix(length(locmaxvec)+1 - as.numeric(factor(V(g$graph)$degree[locmaxvec],labels=1:length(locmaxvec))),
  #                  byrow=T, nrow=locmaxvec, ncol=length(locmaxvec))
  #distmat <- sapply(1:length(locmaxvec),function(i)rep(i,length(locmaxvec)))
  print(paste0("Graph length: ",length(V(g$graph))," | local maxima: ",length(locmaxvec)))
  #print(distmat[1,])
  
  disconnected_locmax <- which(is.infinite(distmat[1,]))
  if(length(disconnected_locmax) > 0){
    distmat[1,disconnected_locmax] <- 0
    warning("Graph is disconnected. All cut values will be zero.")
    locmaxvec <- locmaxvec[c(1,disconnected_locmax[1])]
    distmat <- distmat[,c(1,disconnected_locmax[1])]
    # If the graph is disconnected, perform only one MinCut
    # (between locmaxvec[1] and locmaxvec[disconnected_locmax[1]]), ensuring a cut value of 0.
    #return(xvstlocmax_fast(g=g, allcuts=allcuts, progress_bar=progress_bar))
  }
  distmat[cbind(1:length(locmaxvec),1:length(locmaxvec))] <- Inf
  # Computing the matrix of shortest path distances between every vertex and locmaxvec, assigning Inf to "loops".
  
  if(progress_bar) progress <- txtProgressBar(min = 0, max = length(locmaxvec)-1, style = 3)
  # Setting up the progress bar.
  
  
  while(length(locmaxvec) > 1){
    #    print(paste(locmaxvec[1],"|",locmaxvec[i]))
   #     cat("distmat:",distmat[1,],"\n")
    
    #mindist_idx <- which.min(distances(g_directed, v=locmaxvec[1], to=locmaxvec[-1])) + 1
    mindist_idx <- which.min(distmat[1,])
    # Determines the locmax index of the vertex attaining the minimal distance to locmaxvec[1].
    
     #   cat("locmaxvec:",locmaxvec,"\n")
    #    cat("mindist_idx:",mindist_idx,"\n")
    
    if(allcuts) ij_cut <- stMincuts(g_directed, locmaxvec[1], locmaxvec[mindist_idx])
    # Computes all partitions attaining the ij-MinCut.
    
    if(length(ij_cut$partition1s) == 0){
      # stMincuts doesn't return a partition if the MinCut value is zero.
      
      ij_cut <- min_cut(g$graph, source=locmaxvec[1], target=locmaxvec[mindist_idx], value.only=FALSE)
      # In that case we use igraph::min_cut to determine the st-MinCut partition.
      
      ij_cut$partition1s <- list(ij_cut$partition1)
      # Stores the partition as a one-element list for the for loop below.
      
    }
    
    for(S_ij in ij_cut$partition1s){
      # Goes through all partitions attaining the ij-MinCut.
      
      #print(S_ij)
      if(locmaxvec[1] %in% S_ij) S_ij  <- (V(g_directed)[-S_ij])$tempnr #S_ij <- union(S_ij, already_merged)#c(S_ij, locmaxvec[2:(i-1)])
      else S_ij <- (V(g_directed)[S_ij])$tempnr#S_ij <- setdiff(S_ij, already_merged)
 #     print(S_ij)
      
      ij_xcut <- xcut_partition_S(S_ij, g)
      # Computes the XCut of the partition S_ij.
      
      if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min < ij_min_n){
        ij_min_n <- ij_xcut$ncut_min
        ij_min_ncut <- S_ij
        #        print(paste(locmaxvec[1],"|",locmaxvec[i]))
        #        print(S_ij)
      }
      if(!is.na(ij_xcut$ccut_min) & ij_xcut$ccut_min < ij_min_c){
        ij_min_c <- ij_xcut$ccut_min
        ij_min_ccut <- S_ij
      }
      if(!is.na(ij_xcut$rcut_min) & ij_xcut$rcut_min < ij_min_r){
        ij_min_r <- ij_xcut$rcut_min
        ij_min_rcut <- S_ij
      }
      # Updates the best XCut and the partition that attains it.
    }
    
    
    if(length(locmaxvec) > 2){
      # We now merge vertices locmaxvec[1] and locmaxvec[i] by combining them into one (namely vertex locmaxvec[1]).
      
     # vertices_to_merge <- shortest_paths(g_directed, from=locmaxvec[1], to=locmaxvec[mindist_idx],
    #                                      weights=E(g_directed)$distance, mode="all", output="vpath")$vpath[[1]]
      # Determines the vertices constituting the shortest path between locmaxvec[1] and locmaxvec[mindist_idx] for merging.
      #      cat("vertices_to_merge:",vertices_to_merge,"\n")
      
      mergelist <- (1:length(V(g_directed))) %>% (function(a){a[locmaxvec[mindist_idx]] <- locmaxvec[1];a})
   #       cat("mergelist:",mergelist,"\n")
      g_directed <- contract(g_directed, mergelist, vertex.attr.comb=list(degree="sum", sample="sum", "first"))
      # Contracts all vertices in vertices_to_merge into one (namely locmaxvec[1]).
      
      g_directed <- simplify(g_directed, edge.attr.comb=list(capacity="sum", weight="sum", width="sum", distance="min", "first"))
      # Collapses all resulting multiple edges into one by summing their weights etc., and deletes loops, thus finishing the merging.
      
      if(locmaxvec[mindist_idx] < length(V(g_directed))) g_directed <- delete.vertices(g_directed, locmaxvec[mindist_idx])
      
      #for(i in 2:length(locmaxvec))print(distmat[c(1,mindist_idx),i])
  #    print(distmat)
      
      distmat[1,2:length(locmaxvec)] <- sapply(2:length(locmaxvec), function(i)min(distmat[c(1,mindist_idx),i]))
      distmat <- distmat[-mindist_idx,-mindist_idx]
      # Updates distmat to account for vertices_to_merge being merged into locmaxvec[1] (=vertices_to_merge[1])
      
  #    already_merged <- c(already_merged, locmax[mindist_idx])
      # Updates the already merged vertices.
    }
    locmaxvec <- locmaxvec[-mindist_idx] %>% (function(a){a - (a > mindist_idx)})
    #distmat <- distmat[,-mindist_idx]
    # Prevents locmaxvec[mindist_idx] from being merged again.
    
    if(progress_bar) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
    # Updates progress bar.
    
    #print(i)
    #print(j)
    #print(S_ij)
    #print(ij_xcut$ncut_min)
    
  }
  if(progress_bar) close(progress)
  #    bla <- c()
  #    for(i in 1:(length(V(g$graph))-1)){
  #      for(j in (i+1):length(V(g$graph))){
  #        ij_cut <- stMincuts(g_directed, i, j)
  #        # Computes all partitions attaining the ij-MinCut.
  #        
  # #       for(S_ij in ij_cut$partition1s){
  #    ij_xcut <- xcut_partition_S(S_ij, g)
  #    if(!is.na(ij_xcut$ncut_min) & ij_xcut$ncut_min == ij_min_n){
  #      bla <- c(bla, i, j)
  #  }#}}
  #}
  return(list(rcut = ij_min_rcut, rcut_min = ij_min_r,
              ncut = ij_min_ncut, ncut_min = ij_min_n,
              ccut = ij_min_ccut, ccut_min = ij_min_c))#, bla=bla))
  # Returns a list consisting of all balanced graph cuts and the partitions attaining them.
}
# TODO Add text

manual_multiway_spec_clustering_OLD <- function(g, k=7){
  
  if(is.null(g$removed_vertices)){
    deg_mat_inv <- diag(1/rowSums((g$weight_mat + t(g$weight_mat))))
    L_mat_norm <- diag(length(V(g$graph))) - deg_mat_inv %*% (g$weight_mat + t(g$weight_mat))
#    L_mat_unnorm <- diag(rowSums((g$weight_mat + t(g$weight_mat)))) - (g$weight_mat + t(g$weight_mat))
  } else{
    deg_mat_inv <- diag(1 / rowSums((g$weight_mat + t(g$weight_mat))[-g$removed_vertices, -g$removed_vertices]))
    L_mat_norm <- diag(length(V(g$graph))) -
      deg_mat_inv %*% (g$weight_mat + t(g$weight_mat))[-g$removed_vertices, -g$removed_vertices]
#    L_mat_unnorm <- diag(rowSums((g$weight_mat + t(g$weight_mat))[-g$removed_vertices, -g$removed_vertices])) -
#      (g$weight_mat + t(g$weight_mat))[-g$removed_vertices, -g$removed_vertices]
  }
  # Computes the inverse deg_mat_inv of the degree matrix D and both graph Laplacians L_mat_norm and L_mat_unnorm.
  
  try({
    L_mat_norm_eigen <- eigen(L_mat_norm)
#    L_mat_unnorm_eigen <- eigen(L_mat_unnorm)
    # Eigenvalues and eigenvectors of the Laplacians.
    
    U_mat_norm <- L_mat_norm_eigen$vectors[,(ncol(L_mat_norm_eigen$vectors)-k+1):ncol(L_mat_norm_eigen$vectors)]
#    U_mat_unnorm <- L_mat_unnorm_eigen$vectors[,(ncol(L_mat_unnorm_eigen$vectors)-k+1):ncol(L_mat_unnorm_eigen$vectors)]
    # Computes the matrix U (i.e. the k smallest eigenvectors) for the (Un)Normalized spectral clustering algorithm.
    #save(U_mat_norm, file="umatnorm1.RData")
    print(dim(U_mat_norm))
    return(c(a = U_mat_norm, norm = lapply(2:k, function(i)kmeans(U_mat_norm[,(k-i+1):k], centers=i, nstart=dim(U_mat_norm)[1]+k-i)$cluster)))
#                  kmeans(U_mat_norm[,(k-1):k], centers=k-5, nstart=dim(U_mat_norm)[1]-5)$cluster,
#                norm3 = kmeans(U_mat_norm[,(k-2):k], centers=k-4, nstart=dim(U_mat_norm)[1]-4)$cluster,
#                norm4 = kmeans(U_mat_norm[,(k-3):k], centers=k-3, nstart=dim(U_mat_norm)[1]-3)$cluster,
#                norm5 = kmeans(U_mat_norm[,(k-4):k], centers=k-2, nstart=dim(U_mat_norm)[1]-2)$cluster,
#                norm6 = kmeans(U_mat_norm[,(k-5):k], centers=k-1, nstart=dim(U_mat_norm)[1]-1)$cluster,
#                norm7 = kmeans(U_mat_norm[,(k-6):k], centers=k, nstart=dim(U_mat_norm)[1])$cluster))#,
#                unnorm = kmeans(U_mat_unnorm, centers=k, nstart=dim(U_mat_unnorm)[1])$cluster))
    }, silent=TRUE)
}
# TODO Add text + combine with above spectral clustering function


### <TESTING RANGE> #----------------

resseq <- (1e-10)*1.1^(1:250)
setwd("./RExercises/PhD testing range/")
sample_hep <- read.table("CA-HepPh.txt")
sample_hep <- read.table("Cit-HepPh.txt")
sample_hep <- read.csv("musae_squirrel_edges.csv") + 1
sample_hep <- read.csv("musae_facebook_edges.csv") + 1
sample_hep <- read.csv("./gemsec_facebook/artist_edges.csv") + 1
sample_hep <- read.table("email-Enron.txt") + 1
sample_hep <- read.csv("large_twitch_edges.csv") + 1
g_hep <- simplify(largest_component(graph(c(t(sample_hep)), directed=FALSE)))
length(V(g_hep))
length(E(g_hep))
E(g_hep)$weight <- 1
V(g_hep)$degree <- sapply(1:length(V(g_hep)), function(i)sum(E(g_hep)[.inc(i)]$weight))
st <- Sys.time()
xc_g_hep <- xist_gusfield(g_hep, locmax_by_degree=TRUE)
et <- Sys.time()
print(et - st)
print(xc_g_hep$ncut_min)
length(xc_g_hep$ncut)
leiden_g_hep <- cluster_leiden(g_hep, resolution_parameter=5e-6)
table(leiden_g_hep$membership)

resseq <- (2.9e-6)*(1.001)^(1:50)
resseq1 <- resseq[resseq > 1e-6 & resseq < 3e-6]
resleiden <- which.min(sapply(resseq1, function(res)xcut_multipartition_S(g_hep, cluster_leiden(g_hep, resolution_parameter=res)$membership)$ncut_min))
st <- Sys.time()
xcut_multipartition_S(g_hep, cluster_leiden(g_hep, resolution_parameter=resseq1[resleiden])$membership)$ncut_min
et <- Sys.time()
print(et - st)
table(cluster_leiden(g_hep, resolution_parameter=resseq1[resleiden])$membership)

xc_g_hep_42 <- xist_iterated(g_hep, k=42, normalize=TRUE)






setwd("./RExercises/PhD testing range/Confocal_still_images")
#cagrqc <- read.table("CA-GrQc.txt")
cagrqc0 <- read.table("CA-HepPh.txt")
cagrqc <- cagrqc0[sapply(1:dim(cagrqc0)[1], function(x)cagrqc0[x,1]!=cagrqc0[x,2]),]
#cagrqc <- read.table("email-Eu-core.txt")
ca_nodevec <- sort(unique(c(cagrqc$V1, cagrqc$V2)))
n <- length(ca_nodevec)
cagrqc_new <- cbind(sapply(cagrqc$V1, function(x)which(ca_nodevec==x)), sapply(cagrqc$V2, function(x)which(ca_nodevec==x)))
g_ca0 <- graph(edges=as.vector(t(cagrqc_new)), n=n, directed = FALSE)
g_ca0_ccs <- components(g_ca0)$membership
sort(table(g_ca0_ccs),decreasing=T)

ca_nodevec1 <- (1:n)[g_ca0_ccs == 1]
tempmap <- rep(NA, n)
for(i in 1:length(ca_nodevec1)) tempmap[ca_nodevec1[i]] <- i
ca_tempmat <- matrix(tempmap[cagrqc_new], ncol=2)
cagrqc_new2 <- ca_tempmat[sapply(1:dim(ca_tempmat)[1], function(i)!any(is.na(ca_tempmat[i,]))), ]

g_ca1 <- delete.vertices(g_ca0, which(g_ca0_ccs != 1))
#ca_nodevec2 <- sort(unique(c(cagrqc_new2)))
n <- length(ca_nodevec1)
#cagrqc_new2 <- cagrqc_new[sapply(1:dim(cagrqc_new)[1], function(i)all(g_ca0_ccs[cagrqc_new[i,]]==1)), ]
#cagrqc_new3 <- cbind(sapply(cagrqc_new2[,1], function(x)which(ca_nodevec2==x)), sapply(cagrqc_new2[,2], function(x)which(ca_nodevec2==x)))
g_ca1 <- graph(edges=as.vector(t(cagrqc_new2)), n=n, directed = FALSE)

E(g_ca1)$capacity <- 1
#V(g)$x <- sample[,1]
#V(g)$y <- sample[,2]
#n <- length(V(g_ca1))
V(g_ca1)$label <- sapply(1:n, toString)
#V(g_ca1)$label <- sapply(ca_nodevec, toString)
V(g_ca1)$sample <- rep(1, n)
g_ca1_weightmat <- matrix(rep(0,n^2),n)
for(i in 1:length(E(g_ca1))) g_ca1_weightmat[cagrqc_new2[i,1], cagrqc_new2[i,2]] <- 1
g_ca1_weightmat <- g_ca1_weightmat + t(g_ca1_weightmat)
for(i in 1:n) g_ca1_weightmat[i,i] <- 0
V(g_ca1)$degree <- sapply(1:n, function(i)sum(g_ca1_weightmat[,i]))
V(g_ca1)$local_maximum <- sapply(1:n, function(i)all(g_ca1_weightmat[i,]==0 | V(g_ca1)$degree[i] >= V(g_ca1)$degree))
g_ca2 <- list(graph=g_ca1, weight_mat=g_ca1_weightmat, removed_vertices=c())
#g_ca_xcut2 <- xvstlocmax_fast(g_ca)
a1 <- xvstlocmax_fast(g_ca2)#xvst_iterated(k=4, g=g_ca2)
b1 <- manual_spec_clustering(g_ca2)
delvec <- V(a1$graph)[-a1_xc$ncut]
#length(delvec)
#delvec <- (1:n)[-g_ca2_42$ncut]#g_ca2_42
a1 <- list(graph=delete.vertices(a1$graph, delvec), weight_mat=a1$weight_mat[-delvec, -delvec], removed_vertices=c())
a1_xc <- xvstlocmax_fast(a1)
a1_xc$ncut_min
length(a1_xc$ncut)/length(V(a1$graph))
gc()
st <- Sys.time()
#cc1 <- cluster_leiden(g_ca2$graph, resolution_parameter=1e-5)
a1_xc <- xvstlocmax_fast(g_ca2)
#b1 <- manual_spec_clustering(g_ca2)
et <- Sys.time()
et - st
c1 <- xvst_iterated(k=9, g_ca2, colourvec=collist)
b2 <- manual_multiway_spec_clustering(g_ca2, k=9)
cc1 <- cluster_leiden(g_ca2$graph, resolution_parameter=1e-5)
table(cc1$membership)


fb1 <- rbind(fb1, read.table("./gemsec_facebook_dataset/company_edges.csv", header=TRUE, sep=",") + max(fb1) + 1)
fb1 <- rbind(fb1, read.table("./gemsec_facebook_dataset/athletes_edges.csv", header=TRUE, sep=",") + max(fb1) + 1)
fb1 <- rbind(fb1, read.table("./gemsec_facebook_dataset/government_edges.csv", header=TRUE, sep=",") + max(fb1) + 1)
fb1 <- rbind(fb1, read.table("./gemsec_facebook_dataset/new_sites_edges.csv", header=TRUE, sep=",") + max(fb1) + 1)
fb1 <- rbind(fb1, read.table("./gemsec_facebook_dataset/politician_edges.csv", header=TRUE, sep=",") + max(fb1) + 1)
fb1 <- rbind(fb1, read.table("./gemsec_facebook_dataset/public_figure_edges.csv", header=TRUE, sep=",") + max(fb1) + 1)
fb1 <- rbind(fb1, read.table("./gemsec_facebook_dataset/tvshow_edges.csv", header=TRUE, sep=",") + max(fb1) + 1)




M <- 5
prob_vec <- c(0,0,0,0,0, 0,0,0,0,0, 3,2,0,1,1, 2,4,1,2,3, 0,2,1,3,1)
#prob_vec <- c(0,0,0,0,0, 0,0,0,0,0, 1,1,0,2,3, 3,2,1,4,2, 1,2,1,2,0)
#prob_vec <- c(0,0,2,3,2, 0,0,1,2,1, 0,0,0,2,2, 0,0,1,4,4, 0,0,4,1,0)
prob_vec <- prob_vec/sum(prob_vec)
gsd <- generate_graph(t=1, matrix(prob_vec, nrow=M))
which(V(gsd$graph)$local_maximum)
xc_gsd <- xvstlocmax_fast(gsd, progress_bar=F)
plot_graph_partition(gsd, xc_gsd)
stMincuts(as.directed(gsd$graph, mode="mutual"), 1,4)$value#partition1s[[1]]
xc_gsd_exact <- xcut_exact(gsd)

###

M <- 3
eps <- 0.4
lmbd <- sqrt(1+1.5*eps+eps^2)
tvarmat <- matrix(rep(NA,20), nrow=5, ncol=4)
#prob_vec <- rep(1/M^2, M^2)
#prob_vec <- c(1,1,1,1, 1,eps,eps,1, 1,eps,eps,1, 1,1,1,1)
#prob_vec <- c(rep(1, M), rep(eps, M), rep(1, (M-2)*M))
prob_vec <- c(eps,0,0, 2,2,eps, 0,eps,eps)
prob_vec <- prob_vec/sum(prob_vec)
prob_vec_list <- list(c(eps,0,0, 2,2,eps, eps,eps,eps), c(eps,eps,eps, 1,1,1, eps,eps,eps), c(1,1,1, eps,eps,eps, 1,1,1), rep(1, 9))
#prob_vec_list <- list(c(1,1,1, eps,eps,eps, lmbd,lmbd,lmbd), c(1,1,1, eps,eps,eps, lmbd,lmbd,lmbd), rep(1, 9), rep(1, 9))
smlst <- prob_vec_list
prob_vec_list <- lapply(prob_vec_list, function(x)x/sum(x))
for(i in 1:4){
tlist <- c(1, sqrt(2), 2, sqrt(5), sqrt(8))#, 3, sqrt(10), 4)
for(tind in 1:5){
  g1 <- generate_graph(tlist[tind], matrix(prob_vec_list[[i]], nrow=M))
  #mp <- xcut_exact(g1)
  smlst_trafo <- as.vector(t(matrix(smlst[[i]], M)))
  V(g1$graph)$smmm <- smlst_trafo[smlst_trafo!=0]
  #if(tind==1) print(V(g1$graph)$smmm)
  if(i==1) xmp <- xcut_cov(prob_vec_list[[i]], c(5,8), g1, c(5,8), g1)
  else xmp <- xcut_cov(prob_vec_list[[i]], c(1,4,7), g1, c(1,4,7), g1)
  tvarmat[tind,i] <- xmp$rcut_cov
  #print(paste("i=",i,"| t=",t,"partition=",toString(mp$rcut)))
  if(tind==1 & (i==1 | i==3)) plot_graph_partition(g1, c(1,4,7))
  if(tind==1 && i==2) plot_graph_partition(g1, c(1,2,4,5,7,8))
  if(tind==1 && i==4) plot_graph_partition(g1, c(1,2,4,5))
}
#plot(tlist, tlistvar, type="l")
}
tvarmat

tvar_df <- tidyr::pivot_longer(as.data.frame(rbind(rep(0,4),tvarmat,tvarmat[4,]), row.names=NULL),
                                         cols=1:4, names_to="cuts", values_to="values")
tvar_df_new <- cbind(tvar_df, x=rep(c(0,tlist,4), each=4))
tvar_df_new$xend <- rep(c(tlist,4,NA), each=4)
tvar_df_new[21,3] <- 4
tvar_df_new[17,4] <- 4
gg3 <- ggplot(tvar_df_new, aes(x=x, y=values, xend=xend, yend=values, colour=cuts)) +
  geom_vline(aes(xintercept=x), linetype=2, color="grey") +
  #scale_y_continuous(trans = scales::trans_new(name="almost_log", function(x)log10(x+1e-5), function(x)10^x-1e-5), breaks=c(1e-5,1e-4,1e-3)) + #ylim(c(0,0.0015)) +
  #geom_hline(yintercept = 0.042947, colour="blue", linetype="dashed", show.legend=FALSE) +
  geom_point() + geom_point(aes(x=xend, y=values), shape=1) + geom_segment(key_glyph="path") +
  coord_cartesian(xlim=c(0,3.7)) + labs(x="t", y="Limiting variance") + theme_minimal() +
  scale_color_manual(name="Distribution", values = c("green4", "lawngreen", "yellow", "orangered1"),
                     labels=c(paste("(a)", dQuote("Strange")), "(b) Unimodal", "(c) Bimodal", "(d) Uniform"),
                     guide = guide_legend(override.aes = list(linetype = rep("solid", 4), shape = rep(NA, 4)))) +
  theme(legend.position = "bottom", legend.box.background = element_rect(colour="black"))
  #theme(legend.position = c(0.22,0.6), legend.box.background = element_rect(colour = "black"))
gg3
ggsave("rcut_various_limvar.svg", plot=gg3, width=6, height=4.5)



for(pr in prob_vec_list){
  g3 <- generate_graph(1, matrix(pr, nrow=M))
  print(xcut_exact(g3)$rcut_min)
  print(xcut_exact(g3)$rcut)
  print(xcut_partition_S(c(1,4,7), g3)$rcut_min)
  print("----")
  #plot_graph_partition(g3, c(1,4,7))
}
abc <- 3
pvdummy <- c(1,abc,1, 1,0,1, 1,abc,1)
pvdummy <- pvdummy/sum(pvdummy)
for(tindx in c(1,1.5)){
  g41 <- generate_graph(tindx, matrix(pvdummy, nrow=M))
  print(xcut_exact_indexlist(g41)$rcut)
  print(xcut_cov(pvdummy, c(1,2,3,4), g41, c(1,2,3,4), g41)$rcut_cov)
}
function(abc)return



pv <- c(3,4,2,0,0, 1,3,1,0,0, 1,2,0,0,0, 2,3,1,0,0, 0,2,1,0,0)/26
pv_g <- generate_graph(1.5, matrix(pv, nrow=5))
pv_gc <- xcut_exact(pv_g)
plot_graph_partition(pv_g, pv_gc$ncut)

xmv <- c(0.6,0.7,1.5,2.7,2.2,3.2,2.7,3.6,1.1,1.8,0.5,1.3,0.5,1.6,7.1,9.3,5.0,5.4,6.4,7.6,6.7,8.5,6.6,7.5,4.5,8.4)
ymv <- c(2.2,5.2,3.3,2.6,3.4,4.3,1.2,3.1,4.5,4.8,1.0,2.4,3.5,1.7,1.4,1.5,2.7,3.7,2.3,2.5,3.3,3.7,4.4,4.7,4.2,2.5)
x2mv <- c(3.7,5.5,4.5,8.3,2.2,3.2,2.7,4.4,2.9,1.8,0.9,2.3,3.5,1.6,7.1,5.9,5.0,5.4,6.4,7.6,6.7,7.8,6.6,5.5,4.5,8.4)
y2mv <- c(2.2,1.8,1.3,4.7,3.4,4.3,1.2,3.1,3.0,4.8,2.7,2.4,3.5,1.7,1.6,3.0,2.5,3.8,2.3,2.5,3.3,3.7,4.4,4.7,4.2,2.5)

smplmat <- cbind(x2mv, y2mv)
gsmpl <- generate_sample_graph(1.5, smplmat)
#gsc <- xcut_exact(gsmpl)
gsc <- manual_spec_clustering(gsmpl)
plot_graph_partition(gsmpl, which(gsc$norm==1), vertex_scale=5, edge_scale=5)

##
discretize_sample <- function(smpl, M){
  #partition <- sapply(c(1,2), function(x)seq(a[x], b[x], length.out=M))
  smpl_df <- data.frame(x=smpl[,1], y=smpl[,2])
  partition <- funModeling::discretize_get_bins(data=smpl_df, input=c("x", "y"), n_bins=M)
  smpl_discretized <- funModeling::discretize_df(data=smpl_df, data_bins=partition, stringsAsFactors=TRUE)
  smpl_discretized <- sapply(smpl_discretized, as.numeric)
  sample_discretized <- as.matrix(table(smpl_discretized[,1], smpl_discretized[,2]))
  dimnames(sample_discretized) <- NULL
  return(list(sample_discretized=sample_discretized, bins=partition))
} # Discretizes smpl by dividing the data into M^2 parts (each dimension in M parts)

discretize_sample2 <- function(smpl, nbins=M, xlim, ylim){
  smpl <- smpl[smpl[,1]>=xlim[1] & smpl[,1]<=xlim[2] & smpl[,2]>=ylim[1] & smpl[,2]<=ylim[2],]
  smpl <- cbind(ceiling((smpl[,1]-xlim[1])/xlim[2]*nbins), floor((smpl[,2]-ylim[1])/ylim[2]*nbins))
  return(as.numeric(table(factor(smpl[,1]+smpl[,2]*nbins, levels=1:nbins^2))))
}
twodnormal_disc <- discretize_sample2(twodnormal, xlim=c(0,10), ylim=c(0,10))
print(matrix(twodnormal_disc, M))
ggplot(data.frame(x=twodnormal[,1], y=twodnormal[,2]), aes(x,y)) + stat_bin2d() + theme(legend.position = "none")

eps2 <- 0.1
M <- 3
p_vec <- c(eps2, 1, rep(eps2,5), 1, eps2)/(2+7*eps2)
p_vec2 <- c(rep(eps2,4), 1, rep(eps2,4))/(1+8*eps2)
p_vec25 <- rep(1/M^2,M^2)
M2 <- 9
p_vec3 <- c(rep(eps2,22), 1, rep(eps2,35), 1, rep(eps2,22))/(2+79*eps2)
p_vec35 <- c(rep(eps2,40), 1, rep(eps2,40))/(1+80*eps2)
p_vec4 <- rep(1/81,81)
twodnormal <- rbind(mvtnorm::rmvnorm(50000, mean=c(3,3), sigma=0.6*diag(2)), mvtnorm::rmvnorm(50000, mean=c(7,7), sigma=0.6*diag(2)))
onednormal <- mvtnorm::rmvnorm(100000, mean=c(5,5), sigma=4*diag(2))
p_vec6 <- discretize_sample2(onednormal, nbins=M2, xlim=c(0,10), ylim=c(0,10))
p_vec6 <- c(p_vec6)/sum(p_vec6)
p_vec5 <- discretize_sample2(twodnormal, nbins=M2, xlim=c(0,10), ylim=c(0,10))
p_vec5 <- c(p_vec5)/sum(p_vec5)

t_vec2 <- sort(unique(unlist(sapply(1:(M2-1), function(x)sqrt(x^2+(0:x)^2)))))

#rc_vec_ex <- sapply(t_vec, function(t)xcut_exact(generate_graph(t, matrix(p_vec, nrow=M)))$rcut_min)
rc_vec_st_1 <- sapply(t_vec, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec, nrow=M)))$ccut_min)
rc_vec_st_25 <- sapply(t_vec, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec25, nrow=M)))$ccut_min)
#rc_vec_ex_2 <- sapply(t_vec, function(t)xcut_exact(generate_graph(t, matrix(p_vec2, nrow=M)))$rcut_min)
rc_vec_st_2 <- sapply(t_vec, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec2, nrow=M)))$ccut_min)
rc_vec_st_3 <- sapply(t_vec2, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec3, nrow=M2)))$ncut_min)
rc_vec_st_35 <- sapply(t_vec2, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec35, nrow=M2)))$ncut_min)
rc_vec_st_4 <- sapply(t_vec2, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec4, nrow=M2)))$ncut_min)
rc_vec_st_3_1 <- lapply(t_vec, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec25, nrow=M))))
rc_vec_st_5_1 <- lapply(t_vec2, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec5, nrow=M2))))
rc_vec_st_6_1 <- lapply(t_vec2, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec6, nrow=M2))))
rc_vec_st_5 <- sapply(rc_vec_st_5_1, function(x)x$ccut_min)
rc_vec_st_6 <- sapply(rc_vec_st_6_1, function(x)x$ccut_min)
rc_mat <- cbind(rc_vec_st_3, rc_vec_st_35, rc_vec_st_4)#, rc_vec_st_3, rc_vec_st_4, rc_vec_st_5)
#rc_mat <- cbind(rc_vec_ex, rc_vec_st, rc_vec_ex_2, rc_vec_st_2)

nc_vec_ex <- sapply(t_vec, function(t)xcut_exact(generate_graph(t, matrix(p_vec, nrow=M)))$ncut_min)
nc_vec_st <- sapply(t_vec, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec, nrow=M)))$ncut_min)
nc_vec_ex_2 <- sapply(t_vec, function(t)xcut_exact(generate_graph(t, matrix(p_vec2, nrow=M)))$ncut_min)
nc_vec_st_2 <- sapply(t_vec, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec2, nrow=M)))$ncut_min)
nc_mat <- cbind(nc_vec_ex, nc_vec_st, nc_vec_ex_2, nc_vec_st_2)

cc_vec_ex <- sapply(t_vec, function(t)xcut_exact(generate_graph(t, matrix(p_vec, nrow=M)))$ccut_min)
cc_vec_st <- sapply(t_vec, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec, nrow=M)))$ccut_min)
cc_vec_ex_2 <- sapply(t_vec, function(t)xcut_exact(generate_graph(t, matrix(p_vec2, nrow=M)))$ccut_min)
cc_vec_st_2 <- sapply(t_vec, function(t)xcut_via_stmincut(generate_graph(t, matrix(p_vec2, nrow=M)))$ccut_min)
cc_mat <- cbind(cc_vec_ex, cc_vec_st, cc_vec_ex_2, cc_vec_st_2)


rcuts_df <- tidyr::pivot_longer(as.data.frame(rbind(rep(0,3), rc_mat, rc_mat[length(t_vec2),]), row.names=NULL), cols=1:3, names_to="cuts", values_to="values")
rcuts_df_new <- cbind(rcuts_df, x=rep(c(0,t_vec2,12), each=3))
rcuts_df_new$xend <- rep(c(t_vec2,12,NA), each=3)

ggplot(rcuts_df_new, aes(x=x, y=values, xend=xend, yend=values, colour=cuts)) +
  geom_vline(aes(xintercept=x), linetype=2, color="grey") + ylim(c(0,0.01)) +
  #geom_hline(yintercept = 0.042947, colour="blue", linetype="dashed", show.legend=FALSE) +
  geom_point() + geom_point(aes(x=xend, y=values), shape=1) + geom_segment(key_glyph="path") +
  coord_cartesian(xlim=c(0,3.7)) + labs(x="t", y="Ratio Cut value") + theme_minimal() +
  scale_color_manual(name="Ratio Cut", values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                     labels=c("Exact RCut (bimodal)", "Rvst (bimodal)", "Exact RCut (unimodal)", "Rvst (unimodal)"),
                     guide = guide_legend(override.aes = list(linetype = rep("solid", 4), shape = rep(NA, 4)))) +
  theme(legend.position = c(0.16,0.7), legend.box.background = element_rect(colour = "black"))

gg4 <- ggplot(rcuts_df_new, aes(x=x, y=values, xend=xend, yend=values, colour=cuts)) +
  #geom_vline(aes(xintercept=x), linetype=2, color="grey")+ 
  ylim(c(0,1)) +
  #geom_hline(yintercept = 0.042947, colour="blue", linetype="dashed", show.legend=FALSE) +
  #geom_point() + 
  #geom_point(aes(x=xend, y=values), shape=1) +
  geom_segment(key_glyph="path") +
  coord_cartesian(xlim=c(0,11.2)) + labs(x="t", y="Cheeger Cut value") + theme_minimal() +
  scale_color_manual(name="Cheeger Cut (9x9 grid)", values = c("green4", "orange", "blue"),
                     labels=c("Artif. bimodal", "Artif. unimodal", "Uniform"),
                     guide = guide_legend(override.aes = list(linetype = rep("solid", 3), shape = rep(NA, 3)))) +
  theme(legend.position = c(0.75,0.7), legend.box.background = element_rect(colour = "black"))
gg4
ggsave("ccut_overt_2_3.svg", plot=gg4, width=6, height=4.5)


#-----------
M <- 2
prob_vec <- rep(1/4,4)
prob_vec <- c(3.5/13, 2/13, 3/13, 4.5/13)
#prob_vec <- c(1/3,1/3,1/6,1/6)
g1 <- generate_graph(1, matrix(prob_vec, 2))
four_partitions <- list(c(1), c(2), c(3), c(4), c(1,2), c(1,3), c(1,4))
Sigma_mcut <- sapply(four_partitions, function(x)sapply(four_partitions, function(y)xcut_cov(prob_vec, x, g1, y, g1)$mcut_cov))
compute_tau_i <- function(){
  smpl <- mvtnorm::rmvnorm(1, sigma=Sigma_mcut)
  #return(ifelse(smpl[1] < smpl[2] & smpl[1] < smpl[6], 1, 0)) #tau_1
  return(ifelse(smpl[2] < smpl[1] & smpl[2] < smpl[6], 1, 0)) #tau_2
  #return(ifelse(smpl[6] < smpl[1] & smpl[6] < smpl[2], 1, 0)) #tau_3
}
sum(replicate(100000, compute_tau_i()))/100000

stpair_list <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
tau_mat <- rbind(c(3/8, 3/8, 0, 0, 0, 1/4, 0), c(3/8, 0, 3/8, 0, 1/4, 0, 0), c(1/2, 0, 0, 1/2, 0, 0, 0),
                 c(0, 1/2, 0, 1/2, 0, 0, 0), c(0, 3/8, 0, 3/8, 1/4, 0, 0), c(0, 0, 3/8, 3/8, 0, 1/4, 0))
Sigma_ncut <- sapply(four_partitions, function(x)sapply(four_partitions, function(y)xcut_cov(prob_vec, x, g1, y, g1)$ncut_cov))
M_ncut <- sapply(1:length(stpair_list), function(st)sapply(1:length(stpair_list), function(uv)
  sum(sapply(1:length(four_partitions), function(i)sum(sapply(1:length(four_partitions), function(j)
    sum(tau_mat[st, i] * tau_mat[uv, j] * Sigma_ncut[i, j])))))))
M_ncut

BigM <- 50000
n <- 100000000
stdistvec <- rep(NA, BigM)
altogethervec <- rep(NA, BigM)
atv2 <- rep(NA, BigM)
atv3 <- rep(NA, BigM)
whichpartvec <- c()
for(i in 1:BigM){
  sample_discretized <- matrix(rmultinom(1, n, prob_vec), 2)/n
  g2 <- generate_graph(1, sample_discretized)
  altv_cut <- xcut_via_stmincut(g2)
  altogethervec[i] <- sqrt(n) * (altv_cut$ncut_min - xcut_partition_S(altv_cut$ncut, g1)$ncut_min)
#  for(stpair in stpair_list){
#    stMincuts(as.directed(g2$graph, mode="mutual"), stpair[1], stpair[2])$partition1s
#  }
  atv2[i] <- n * (xcut_partition_S(c(1,2), g2)$ncut_min - xcut_partition_S(c(1,2), g1)$ncut_min)^2
  atv3[i] <- n * (xcut_partition_S(c(1,3), g2)$ncut_min - xcut_partition_S(c(1,3), g1)$ncut_min)
  stsample <- stMincuts(as.directed(g2$graph, mode="mutual"), 2, 4)$partition1s
  stsample_ncuts <- sapply(stsample, function(S)xcut_partition_S(S, g2)$ncut_min)
  #print(stsample[[which.min(stsample_ncuts)]])
  #print(min(stsample_ncuts))
  #print(xcut_partition_S(stsample[[which.min(stsample_ncuts)]], g1)$ncut_min)
  #whichpartvec <- c(whichpartvec, stsample)
  whichpartvec <- c(whichpartvec, do.call(function(x=altv_cut$ncut){#stsample[[which.min(stsample_ncuts)]]){
    attributes(x) <- NULL
    #ifelse(isTRUE(all.equal(x,c(1)))|isTRUE(all.equal(sort(x),c(2,3,4))), 1, ifelse(isTRUE(all.equal(x,c(2)))|isTRUE(all.equal(sort(x),c(1,3,4))), 2, ifelse(isTRUE(all.equal(sort(x),c(1,3)))|isTRUE(all.equal(sort(x),c(2,4))), 6, 99)))
    ifelse(isTRUE(all.equal(x,c(1)))|isTRUE(all.equal(sort(x),c(2,3,4))), 1, ifelse(isTRUE(all.equal(x,c(4)))|isTRUE(all.equal(sort(x),c(1,2,3))), 2, ifelse(isTRUE(all.equal(sort(x),c(1,2)))|isTRUE(all.equal(sort(x),c(3,4))), 6, 99)))
    },list()))
  #stdistvec[i] <- sqrt(n) * (min(stsample_ncuts) - xcut_partition_S(stsample[[which.min(stsample_ncuts)]], g1)$ncut_min)
  stdistvec[i] <- sqrt(n) * (min(stsample_ncuts) - xcut_partition_S(stsample[[which.min(stsample_ncuts)]], g1)$ncut_min)
}
hist(stdistvec, breaks=100)
hist(altogethervec, breaks=100)
hist(atv2, breaks=100)
hist(atv3, breaks=100)

stdistlimvec <- c()
stdistlimvec2 <- c()
stdistlimvec3 <- c()
stdistvec2 <- c()
altogetherlimvec <- c()
for(i in 1:BigM){
  p_tau <- runif(1)
  ncut_lim_sample <- mvtnorm::rmvnorm(1, sigma=Sigma_ncut)
  if(p_tau < 3/8){
    stdistlimvec[i] <- ncut_lim_sample[1]
  } else if(p_tau < 6/8){
    stdistlimvec[i] <- ncut_lim_sample[2]
  } else {
    stdistlimvec[i] <- ncut_lim_sample[6]
  }
  stdistlimvec2[i] <- min(ncut_lim_sample[1], ncut_lim_sample[4], ncut_lim_sample[5], ncut_lim_sample[6])
  altogetherlimvec[i] <- ncut_lim_sample[4]#min(ncut_lim_sample[5], ncut_lim_sample[6])
  sample_discretized <- matrix(rmultinom(1, n, prob_vec), 2)/n
  g2 <- generate_graph(1, sample_discretized)
  stdistvec2[i] <- sqrt(n) * (xcut_partition_S(c(1,3), g2)$ncut_min - xcut_partition_S(c(1,3), g1)$ncut_min)
  stdistlimvec3[i] <- ncut_lim_sample[6]
}
#hist(stdistlimvec, breaks=200)
hist(stdistlimvec2, breaks=100)
hist(stdistvec2, breaks=200)
#hist(stdistlimvec3, breaks=200)
hist(altogetherlimvec, breaks=200)

ncutlimcovmat <- cbind(c(5,-18,-11,-2), c(-18,5,-2,-11), c(-11,-2,5,-18), c(-2,-11,-18,5))
ncuttestlimvec <- rep(NA, BigM)
for(i in 1:BigM){
  sample_vec <- rmultinom(1, n, prob_vec)/n - prob_vec
  ncuttestlimvec[i] <- n*t(sample_vec)%*%ncutlimcovmat%*%sample_vec
}
hist(ncuttestlimvec, breaks=100)

##---CONFOCAL TESTING RANGE
t <- 1.5
#setwd("./RExercises/PhD testing range/Confocal_still_images")
timearray96 <- matrix(rep(NA,10*21), 10, 21)

for(i in 1:10){
for(j in c(1:21)){
  M <- c(8,9,12,14,18,21,24,28,36,42)[i]
  cfr1 <- raster::raster(paste0("./s_C001Tubulin", ifelse(j<17,j,j+1), ".tif"))
  cfr <- crop(cfr1, extent(cfr1, 9, 512, 1, 504))
  
#  pt <- discretize_and_xcut_confocal(cfr, M=M)
#  collist <- c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1", "purple")
#  V(pt$g$graph)[unlist(pt$partition)]$color <- collist[1]
#  V(pt$g$graph)[-unlist(pt$partition)]$color <- collist[2]
#  confocal_vertex_size <- t(pt$cf_disc)
#  V(pt$g$graph)$vertex.size <- 15*confocal_vertex_size[confocal_vertex_size>0]/max(confocal_vertex_size) + 3
#  edge_weights <- E(pt$g$graph)$weight
#  par(cex=0.7, mai=c(0.5,0.5,0.3,0.3))
#  plot(pt$g$graph, vertex.label=NA, edge.width = 10*edge_weights/max(edge_weights)+1, vertex.size = V(pt$g$graph)$vertex.size)
  
  #confocal_discretized <- as.matrix(aggregate(cfr, dim(cfr)[1]/M, sum, expand=FALSE)) %>% (function(a)a/sum(a))
  #gc <- generate_graph_OLD(t=1.5, confocal_discretized, locmax_degree=TRUE)
  gc <- generate_tif_grid_graph(cfr, M=M, t=1.5)
  gc <- delete.vertices(gc, which(V(gc)$sample == 0))
  
  #distmat <- as.dist(1/gc$weight_mat)
  wm <- compute_weight_matrix(gc)
#  k <- 2
  start_time <- Sys.time()
  #a <- dbscan::dbscan(distmat, eps=1e6)
  #a <- manual_spec_clustering(gc, method="kmeans")
  #a <- sClust::VonLuxburgSC(gc$weight_mat, K=2)
  a <- sClust::VonLuxburgSC(wm, K=2)
  #a <- xist_gusfield(gc, progress_bar = FALSE)
  #a <- xvstlocmax_fast(gc, progress_bar = FALSE, mergeByDegree=FALSE)
  #a <- cluster_leiden(gc$g, resolution_parameter = 1e-5)
  #tm_xvst <- system.time(xcut_via_stmincut(gc))
#  cfr_xcut_iterated <- xvst_iterated(k, pt$g, collist)
  stop_time <- Sys.time()
  timearray96[i,j] <- difftime(stop_time, start_time, units="secs")[[1]]
#  print(stop_time - start_time)
  #timearray21[i,j] <- as.numeric(tm_xvst[3])#stop_time - start_time
  print(paste0("(",i,", ",j,")"))
}}
load("xvst_algs_timearrays6.RData")
#save(timearray11, timearray21, timearray31, file="xvst_algs_timearrays3.RData")
#save(timearray11, timearray21, timearray31, timearray71, timearray81, timearray91, file="xvst_algs_timearrays6.RData")

ta_xvst <- data.frame(m=rep(NA,21*8), value=rep(NA,21*8), alg=rep("alg1",21*8))
ta_xvstlocmax <- data.frame(m=rep(NA,21*8), value=rep(NA,21*8), alg=rep("alg2",21*8))
ta_xvstlocmaxmerge <- data.frame(m=rep(NA,21*8), value=rep(NA,21*8), alg=rep("alg3",21*8))
ta_specclust <- data.frame(m=rep(NA,21*8), value=rep(NA,21*8), alg=rep("alg4",21*8))
ta_leiden <- data.frame(m=rep(NA,21*8), value=rep(NA,21*8), alg=rep("alg5",21*8))
ta_dbscan <- data.frame(m=rep(NA,21*8), value=rep(NA,21*8), alg=rep("alg6",21*8))
Mvec <- c(8,9,12,14,18,21,24,28,36,42)^2
for(i in 1:8){for(j in 1:21){
  ta_xvst[(i-1)*21+j,1:2] <- c(Mvec[i], timearray11[i,j])
  ta_xvstlocmax[(i-1)*21+j,1:2] <- c(Mvec[i], timearray95[i,j])
  ta_xvstlocmaxmerge[(i-1)*21+j,1:2] <- c(Mvec[i], timearray95[i,j])
  ta_dbscan[(i-1)*21+j,1:2] <- c(Mvec[i], timearray71[i,j])
  ta_specclust[(i-1)*21+j,1:2] <- c(Mvec[i], timearray81[i,j])
  ta_leiden[(i-1)*21+j,1:2] <- c(Mvec[i], timearray91[i,j])
}}
ta_arr <- rbind(ta_xvst, ta_xvstlocmax, ta_xvstlocmaxmerge, ta_specclust, ta_leiden, ta_dbscan)
library(scales)
library(latex2exp)

ggplot(ta_arr, aes(x=m, y=value, group=m)) +
  geom_boxplot(data=ta_xvst, aes(fill="Basic Xvst")) +
  geom_boxplot(data=ta_xvstlocmax, aes(fill="Spec. Clust.")) +
  geom_boxplot(data=ta_xvstlocmaxmerge, aes(fill="Xist")) +
  scale_x_continuous(trans='log2', breaks=Mvec,
                     labels=trans_format("sqrt", math_format(.x^2))) + theme_minimal() +
  scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x="Grid size m", y="Running time (sec.)") +
  #geom_abline(aes(intercept = -17.5, slope = 3, fill="asdf"), color="blue", linetype="dashed", size=.8) +
  #geom_abline(aes(intercept = -15, slope = 2.3, fill="sdfsd"), color="purple", linetype="dashed", size=.8) +
  geom_abline(data=data.frame(slope=c(3,2.3,1.6), intercept=c(-17.5,-15,-12), id=c("A", "B", "C")), linetype="dashed", size=.8,
              mapping=aes(slope=slope, intercept=intercept, color=id)) +
  scale_fill_manual(name="Algorithm",
                    #labels=c("Algorithm 4.1", "Algorithm 4.2", "Algorithm 4.3"),
                    values=c("blue", "purple", "red")) +
  scale_color_manual(name="Runtime",
                     labels=unname(TeX(c("$m^3\\;$", "$m^{2.3}$", "$m^{1.6}$"))),
                     values=c("A", "B", "C"),
                     guide=guide_legend(override.aes = list(color=c("blue", "purple", "red")))) +
  theme(legend.position = c(0.15,0.7), legend.box.background = element_rect(colour = "black",fill="white")) +
  #guides(fill=guide_legend(order=1), color=guide_legend(order=0))
  geom_abline(aes(intercept = -17.5, slope = 3), color="blue", linetype="dashed", size=.8) +
  geom_abline(aes(intercept = -15, slope = 2.3), color="purple", linetype="dashed", size=.8) +
  geom_abline(aes(intercept = -12, slope = 1.6), color="red", linetype="dashed", size=.8)

ta_xvst <- data.frame(m=rep(NA,21*10), value=rep(NA,21*10), alg=rep("alg1",21*10))
ta_xvstlocmax <- data.frame(m=rep(NA,21*10), value=rep(NA,21*10), alg=rep("alg2",21*10))
ta_xvstlocmaxmerge <- data.frame(m=rep(NA,21*10), value=rep(NA,21*10), alg=rep("alg3",21*10))
Mvec <- c(8,9,12,14,18,21,24,28,36,42)^2
for(i in 1:10){for(j in 1:21){
  ta_xvst[(i-1)*21+j,1:2] <- c(Mvec[i], timearray11[i,j])
  ta_xvstlocmax[(i-1)*21+j,1:2] <- c(Mvec[i], timearray21[i,j])
  ta_xvstlocmaxmerge[(i-1)*21+j,1:2] <- c(Mvec[i], timearray31[i,j])
}}
ta_arr <- rbind(ta_xvst, ta_xvstlocmax, ta_xvstlocmaxmerge)

ggplot(ta_arr, aes(x=m, y=value, group=m)) +
  geom_boxplot(data=ta_xvst, aes(fill="Alg. 4.1")) +
  geom_boxplot(data=ta_xvstlocmax, aes(fill="Alg. 4.3")) +
  geom_boxplot(data=ta_xvstlocmaxmerge, aes(fill="Alg. 4.4")) +
  scale_x_continuous(trans='log2', breaks=Mvec,
                     labels=trans_format("sqrt", math_format(.x^2))) + theme_minimal() +
  scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x="Grid size m", y="Running time (sec.)") + coord_cartesian(ylim=c(2^(-10),2^11)) +
  geom_abline(data=data.frame(slope=c(3.1,2.8,1.6), intercept=c(-19.8,-24.9,-17.5), id=c("blue", "purple", "red")), linetype="dashed", size=.8,
              mapping=aes(slope=slope, intercept=intercept, color=id)) +
  scale_fill_manual(name="Algorithm",
                    labels=unname(TeX(c("Xvst", "$Xvst^{loc}$", "$Xvst_{merge}^{loc}$"))),
                    values=c("blue", "purple", "red")) +
  scale_color_manual(name="Runtime",
                     labels=unname(TeX(c("$m^{3.1}\\;$", "$m^{2.8}$", "$m^{1.6}$"))),
                     values=c("blue", "purple", "red"),
                     guide=guide_legend(override.aes = list(color=c("blue", "purple", "red")))) +
  theme(legend.text.align = 0, legend.position = c(0.18,0.85), legend.box="horizontal", 
        legend.box.background = element_rect(colour = "black",fill="white")) +
  guides(fill=guide_legend(order=1), color=guide_legend(order=2))
#geom_abline(aes(intercept = -23.5, slope = 2.9), color="blue", linetype="dashed", size=.8) +
#geom_abline(aes(intercept = -17.5, slope = 1.5), color="purple", linetype="dashed", size=.8)

ggplot(ta_arr, aes(x=m, y=value, group=m)) +
  geom_boxplot(data=ta_xvst, aes(fill="Alg. 4.1")) +
  #geom_boxplot(data=ta_xvstlocmax, aes(fill="Alg. 4.3")) +
  geom_boxplot(data=ta_xvstlocmaxmerge, aes(fill="Alg. 4.4")) +
  scale_x_continuous(trans='log2', breaks=Mvec,
                     labels=trans_format("sqrt", math_format(.x^2))) + theme_minimal() +
  scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x="Number of vertices n", y="Running time (sec.)") + coord_cartesian(ylim=c(2^(-10),2^11)) +
  #geom_abline(data=data.frame(slope=c(3.1,1.6), intercept=c(-19.8,-17.5), id=c("blue", "red")), linetype="dashed", size=.8,
  geom_abline(data=data.frame(slope=c(3.1,2.1), intercept=c(-19.2,-20.0), id=c("blue", "red")), linetype="dashed", size=.8,
              mapping=aes(slope=slope, intercept=intercept, color=id)) +
  scale_fill_manual(name="Algorithm",
                    labels=c("Xvst", "Xist"),
                    values=c("blue", "red")) +
  scale_color_manual(name="Runtime",
                     labels=unname(TeX(c("$n^{3.1}\\;$", "$n^{2.1}$"))),
                     values=c("blue", "red"),
                     guide=guide_legend(override.aes = list(color=c("blue", "red")))) +
  theme(legend.text.align = 0, legend.position = c(0.75,0.14), legend.box="horizontal", 
        legend.box.background = element_rect(colour = "black",fill="white"), legend.text=element_text(size=14),
        text = element_text(size=rel(4)), axis.text=element_text(size=14), panel.background = element_rect(fill="gray98", colour="white")) +
  guides(fill=guide_legend(order=1), color=guide_legend(order=2))
#geom_abline(aes(intercept = -23.5, slope = 2.9), color="blue", linetype="dashed", size=.8) +
#geom_abline(aes(intercept = -17.5, slope = 1.5), color="purple", linetype="dashed", size=.8)

ggplot(ta_arr, aes(x=m, y=value, group=m)) +
  geom_boxplot(data=ta_xvst, aes(fill="A1")) +
  #geom_boxplot(data=ta_xvstlocmax, aes(fill="Alg. 4.3")) +
  geom_boxplot(data=ta_xvstlocmaxmerge, aes(fill="A2")) +
  geom_boxplot(data=ta_specclust, aes(fill="A3")) +
  geom_boxplot(data=ta_leiden, aes(fill="A4")) +
  geom_boxplot(data=ta_dbscan, aes(fill="A5")) +
  scale_x_continuous(trans='log2', breaks=Mvec,
                     labels=trans_format("sqrt", math_format(.x^2))) + theme_minimal() +
  scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x="Number of vertices n", y="Running time (sec.)") + coord_cartesian(ylim=c(2^(-10),2^11)) +
  #geom_abline(data=data.frame(slope=c(3.1,1.6), intercept=c(-19.8,-17.5), id=c("blue", "red")), linetype="dashed", size=.8,
  geom_abline(data=data.frame(slope=c(3.1,0.9,2.3,1.5,1.1), intercept=c(-19.2,-17.2,-21.7,-19,-11.3), id=c("blue", "orangered1", "orange1", "yellow2", "mediumpurple2")), linetype="dashed", size=.8,
              mapping=aes(slope=slope, intercept=intercept, color=id)) +
  scale_fill_manual(name="Algorithm",
                    labels=c("Xvst", "Xist", "Spec. Clust.", "Leiden", "DBSCAN"),
                    values=c("blue", "orangered1", "orange1", "yellow2", "mediumpurple2")) +
  scale_color_manual(name="Runtime",
                     labels=unname(TeX(c("$n^{3.1}\\;$", "$n^{1.1}\\;$", "$n^{2.3}\\;$", "$n^{0.9}\\;$", "$n^{1.5}$"))),
                     values=c("blue", "orangered1", "orange1", "yellow2", "mediumpurple2"),
                     guide=guide_legend(override.aes = list(color=c("blue", "orangered1", "orange1", "yellow2", "mediumpurple2")))) +
  theme(legend.text.align = 0, legend.position = c(0.2,0.81), legend.box="horizontal", 
        legend.box.background = element_rect(colour = "black",fill="white"), legend.text=element_text(size=14),
        text = element_text(size=rel(4)), axis.text=element_text(size=14), panel.background = element_rect(fill="gray98", colour="white")) +
  guides(fill=guide_legend(order=1), color=guide_legend(order=2))
#geom_abline(aes(intercept = -23.5, slope = 2.9), color="blue", linetype="dashed", size=.8) +
#geom_abline(aes(intercept = -17.5, slope = 1.5), color="purple", linetype="dashed", size=.8)



m <- 5
meanvec2 <- rep(NA, 20)
for(m in 2:20){
n <- 100000
BigM <- 100000
pvec <- rep(1/m^2, m^2)
#Msigm <- diag(pvec) - pvec%*%t(pvec)
#Z_smpl <- mvtnorm::rmvnorm(100000, sigma=Msigm)

exptlocmaxvec <- rep(NA, BigM)
for(i in 1:BigM){
  sample_discretized <- matrix(rmultinom(1, n, pvec), m)
  g1 <- generate_graph(1.5, sample_discretized)
  exptlocmaxvec[i] <- length(which(V(g1$graph)$local_maximum))
  #if((i*10)%%BigM == 0) print(i)
  if(i==BigM) print(i)
}
meanvec2[m] <- mean(exptlocmaxvec) # Estimated values: m=2: 1.336879, m=3: 2.54353, m=4: 4.15877, m=5: 6.18176
}
save(meanvec2, file="expected_nr_loc_max_2_20_t15.Rdata")
plot(2:20, meanvec[1:19], type="l")
lines(2:20, (2:20)^2*0.4223267/2 - (2:20)*0.1104657 - 0.2932241, type="l", col="green")
lines(2:19, meanvec2[1:18], type="l", col="red")
lines(2:20, (2:20)^2/8.5-1, type="l", col="blue")
lines(2:20, (2:20)*log(2:20)^2/2, col="yellow")

plot(log(2:20), log(meanvec[1:19]), type="l")
lines(log(2:20), log(meanvec2[1:19]), type="l", col="red")

### </TESTING RANGE>



### PART II: EXEMPLARY USE OF THE FUNCTIONS FROM PART I ###


# Cutting an NIH3T3 cell sample repeatedly using Xvst ---------------------

M <- 16
t <- 1.5
# Defining grid size M and neighbourhood distance t.

setwd("./RExercises/PhD testing range/")
#setwd("/path/to/NIH3T3/data")
# Setting the correct path to the NIH3T3 dataset

cf1_13 <- raster("./Confocal_still_images/s_C001Tubulin13.tif")
# Loads the confocal image of sample #13. The NIH3T3_Data folder contains ten confocal images in total.
# The interested reader is invited to load some of these images to explore the NIH3T3 dataset.

discretize_and_xcut_confocal <- function(cf, M=16, t=1.5, cut="ncut"){

  confocal_discretized <- as.matrix(aggregate(cf, dim(cf)[1]/M, sum, expand=FALSE))
  # Discretizes the raster image into an M x M matrix.
  
  n <- sum(confocal_discretized)
  confocal_discretized <- confocal_discretized/n
  # Normalizes the observations for convenience. The scaling does not matter with the confocal dataset.
  
  gc <- generate_grid_graph(confocal_discretized, t=t)
  # Generates a extended graph object out of the discretized confocal sample.
  
  if(cut=="rcut"){
    nvst_confocal <- xvstlocmax_fast(gc)$rcut
  }else if(cut=="ncut"){
    nvst_confocal <- xist_gomoryhu(gc)#xvstlocmax_fast(gc)$ncut#xist(gc)$ncut
  }else if(cut=="ccut"){
    nvst_confocal <- xvstlocmax_fast(gc)$ccut
  }else {
    print("Error: 'cut' must be either one of 'rcut', 'ncut' or 'ccut'.")
    return(list(g=gc, partition=c(), cf_disc=confocal_discretized))
  }
  # Computes the appropriate XCut of the discretized confocal sample.
  
  confocal_vertex_size <- t(confocal_discretized)
  plot_graph_partition(gc, nvst_confocal, edge_scale=10)#,
#                       vertex_size=15*confocal_vertex_size[confocal_vertex_size>0]/max(confocal_vertex_size) + 3)
  # Plots the discretized and cut confocal sample.
  
  return(list(g=gc, partition=nvst_confocal, cf_disc=confocal_discretized))
  # Returns an extended graph object, NCut partition and discretized sample of the confocal sample.
}
# Discretizes and computes the NCut of a raster image (here: of a confocal .tif cell sample).
# M needs to divide the number of pixels in each direction without remainder, and the image has to be quadratic.
# Input: raster image cf, grid size M, distance t.
# Output: List containing an extended graph object, NCut partition and discretized sample of the confocal sample

pt <- discretize_and_xcut_confocal(cf1_13, M=32, cut="ncut")
collist <- c("orangered1", "lawngreen", "yellow", "purple","green4", "darkgoldenrod2", "blue", "cyan", "brown3")
#collist <- c("lawngreen", "orangered1", "green4", "yellow", "darkgoldenrod2", "purple")
V(pt$g$graph)[unlist(pt$partition)]$color <- collist[1]
V(pt$g$graph)[-unlist(pt$partition)]$color <- collist[2]
confocal_vertex_size <- t(pt$cf_disc)
V(pt$g$graph)$vertex.size <- 15*confocal_vertex_size[confocal_vertex_size>0]/max(confocal_vertex_size) + 3
edge_weights <- E(pt$g$graph)$weight
par(cex=0.7, mai=c(0.5,0.5,0.3,0.3))
plot(pt$g$graph, vertex.label=NA, edge.width = 10*edge_weights/max(edge_weights)+1, vertex.size = V(pt$g$graph)$vertex.size)
# Plots the NCut partition of the discretized confocal sample pt.

xvst_part_of_graph <- function(partition, g, colour, draw_plot=TRUE){
  g1 <- delete_vertices(g$graph, V(g$graph)[-partition])
  ncg1 <- list(graph=g1, removed_vertices=g$removed_vertices, weight_mat=g$weight_mat[partition, partition])
  V(ncg1$graph)$label <- as.character(1:length(partition))
  if(length(V(ncg1$graph)) > 1) ncg1_xcut <- xvstlocmax_fast(ncg1,progress_bar=F,mergeByDegree=T)#xcut_via_stmincut_on_local_maxima(ncg1)
  else ncg1_xcut <- list(mcut=c(),rcut=c(),ncut=c(),ccut=c(),mcut_min=Inf,rcut_min=Inf,ncut_min=Inf,ccut_min=Inf)
  print(paste("JUST CUT:",ncg1_xcut$ncut_min))
  real_stuff <- (sort((as.numeric(V(g$graph)$label))[partition]))[ncg1_xcut$ncut]
  V(g$graph)[real_stuff]$color <- colour
  if(draw_plot) plot(g$graph, vertex.label=NA, edge.width = 10*edge_weights/max(edge_weights)+1, vertex.size = V(g$graph)$vertex.size)
  return(list(g=g, partition1=real_stuff, partition2=setdiff(partition, real_stuff), xcut=ncg1_xcut))
}
# Divides the partition nodes of the graph g$graph in two using xcut_via_stMinCuts and colours them appropriately (for now, using NCut).
# Returns the coloured graph, and the two partitions, allowing them to be divided further by another call of the function.

xvst_iterated <- function(k, g, colourvec=NULL){
  partitions <- list(1:length(V(g$graph)))
  partitions_list <- list()
  ncut_vec <- c()
  index_best <- 1
  if(is.null(colourvec)) colourvec <- c("green3", "blue", sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)], k-2))
#  V(g$graph)$color <- colourvec[1]
  glist <- list(xvst_part_of_graph(1:length(V(g$graph)), g, "black", draw_plot=FALSE))
  for(i in 1:(k-1)){
    latest_cut <- glist[[index_best]]
#    V(g$graph)[latest_cut$partition1]$color <- colourvec[i+1] # Colours the original graph
    partitions[c(index_best, i+1)] <- list(latest_cut$partition1, latest_cut$partition2)
    partitions_list[[i]] <- partitions
    glist[c(index_best, i+1)] <- list(xvst_part_of_graph(latest_cut$partition1, g, "black", draw_plot=FALSE), xvst_part_of_graph(latest_cut$partition2, g, "black", draw_plot=FALSE))
    ncut_vec[c(index_best, i+1)] <- c(glist[[index_best]]$xcut$ncut_min, glist[[i+1]]$xcut$ncut_min)
    index_best <- which.min(ncut_vec)
    print(ncut_vec)
    print(partitions[index_best])
  }
  for(i in 1:k) V(g$graph)[partitions[[i]]]$color <- colourvec[i]
  V(g$graph)$vertex.size <- 15*V(g$graph)$sample/max(V(g$graph)$sample) + 3
  edge_weights <- E(g$graph)$weight
#  par(cex=0.7, mai=c(0.5,0.5,0.3,0.3))
#  plot(g$graph, vertex.label=NA, edge.width = 10*edge_weights/max(edge_weights)+1, vertex.size = V(g$graph)$vertex.size)
  return(list(g=g, partitions=partitions, partitions_list=partitions_list))
}
# Divides the graph g$graph into k clusters using iterated application of the Xvst algorithm. The clusters are coloured using colourlist, a list of k colours.
# Returns the coloured extended graph object g as well as a list of all k partitions.

#k <- 7
#start_time <- Sys.time()
pt <- discretize_and_xcut_confocal(cf1_13, M=16, t=1.5, cut="none")
aa <- manual_multiway_spec_clustering(pt$g, k=7)
#V(pt$g$graph)$color <- collist[manual_multiway_spec_clustering(pt$g,k=6)$norm]
plot(pt$g$graph, vertex.label=NA, edge.width = 10*edge_weights/max(edge_weights)+1, vertex.size = V(pt$g$graph)$vertex.size)
cf1_13_xcut_iterated <- xvst_iterated(k=7, pt$g, collist)
#stop_time <- Sys.time()
#print(stop_time - start_time)
# Dividing the cell image cf1_13 into k partitions, colouring them according to collist.

colr <- colorRampPalette(RColorBrewer::brewer.pal(9, 'YlOrRd'))
rasterVis::levelplot(flip(t(cf1_13), 2),
                     margin = FALSE,
                     colorkey = F,#list(space='bottom', labels=list(at=seq(0,5000,length.out=6))),
                     par.settings = list(axis.line=list(col='transparent')),
                     scales = list(draw=FALSE),
                     col.regions = colr,
                     at = seq(0,5000,length.out=101))
# Plotting the cell sample cf1_13 with proper scaling and nice colouring.

put_square_at_partition <- function(partitionvec, hpixels, hM=M){
  vpixels <- hpixels
  vM <- hM
  # Dummy variables for future implementations with rectangular images. WIP.
  
  locations <- sapply(partitionvec, function(x)c(hpixels/hM*((x-1)%/%vM + 0.5), vpixels/vM*((x-1)%%vM + 0.5)))
  return(rgeos::gBuffer(SpatialPoints(t(locations)), width=0.5*hpixels/hM, quadsegs=1, capStyle="SQUARE"))
}
# Returns SpatialPolygons for squares around the locations specified by the discretized locations given by partitionvec and hM to overlay on top of an image of hpixels x vpixels
# Currently only computes squares. WIP.

uncoloured_squares_list <- lapply(cf1_13_xcut_iterated$partitions, function(p)put_square_at_partition(p, dim(cf1_13)[1], M))
# Returns a list of all squares (that belong to a partition) for every partition that cf1_13 has been divided into.
k=6
uncoloured_squares_list <- lapply(cf1_13_xcut_iterated$partitions_list[[k-1]], function(p)put_square_at_partition(p, dim(cf1_13)[1], M))

colr2 <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Greys'))
rasterVis::levelplot(flip(t(cf1_13), 2),
                     margin = FALSE, colorkey=FALSE,
                     par.settings = list(axis.line=list(col='transparent')),
                     scales = list(draw=FALSE),
                     col.regions = colr2) +
  latticeExtra::layer(lapply(1:k, function(idx)sp.polygons(uncoloured_squares_list[[idx]], fill=collist[idx], col=NA, alpha=0.5)))
# Plots the image together with the partitions (which are coloured appropriately)

cf1_13_cropped <- crop(cf1_13, extent(cf1_13, 9, 512, 1, 504))
# Cropping cf1_13 so that it can be divided more easily by different M without remainder. This is helpful for comparing
# how the cell sample is cut for different grid sizes M. Be aware of the dimensions of the sample you want to crop!

pt_M24 <- discretize_and_xcut_confocal(cf1_13_cropped, M = 12)
# Trying out a different grid size M. Note that the number of pixels in each direction is evenly divisible by 12.


# Demonstrating some basic properties of some cuts ------------------------

M <- 4
t <- 1
n <- 200
eps <- 0.4
# Defining the necessary parameters. eps corresponds to epsilon from Chapter 4 of the thesis.

#prob_vec <- c(rep(1/M^2, M^2))
prob_vec <- c(rep(1, M), rep(eps, M), rep(1, (M-2)*M))
prob_vec <- prob_vec / sum(prob_vec)
# Two probability vectors: one uniform, one with a slight "epsilon-gap" in between (indicated by two rows of eps).

g_true <- generate_graph(t, matrix(prob_vec, nrow = M))
# Computing the true underlying graph structure using the probability vector prob_vec

sample_discretized <- matrix(rmultinom(1, n, prob_vec), nrow = M)
g_sd <- generate_graph(t, sample_discretized)
gMC_mcut <- min_cut(g_sd$graph, value.only = FALSE)
plot_graph_partition(g_sd, gMC_mcut$partition1)
# Consider the uselessness of MinCut in practice.

xvst_g_sd <- xcut_via_stmincut(g_sd)
plot_graph_partition(g_sd, xvst_g_sd$ncut, main = "Nvst partition")
# Now cut the empirical graph built from sample_discretized using the Xvst algorithm.

xcut_g_sd <- xcut_exact(g_sd)
plot_graph_partition(g_sd, xcut_g_sd$rcut, main = "True Ratio Cut partition")
plot_graph_partition(g_sd, xcut_g_sd$ncut, main = "True NCut partition")
plot_graph_partition(g_sd, xcut_g_sd$ccut, main = "True Cheeger Cut partition")
# Computes and cut the exact XCut of sample_discretized. Compare the resulting plots.


# Comparing Rvst and unnormalized spectral clustering ---------------------

BigM <- 2000
t <- 1
nvec <- c(50, 100, 150, 200, 500)
rcut_sample <- matrix(rep(NA,BigM*2*length(nvec)), nrow=BigM)
# This matrix will contain the empirical cut values of BigM samples computed using Xvst and spectral clustering for
# several values of n contained within nvec.

progress <- txtProgressBar(min = 0, max = length(nvec)*BigM, style = 3)
# Sets up a little progress bar.

for(i in 1:length(nvec)){
  for(m in 1:BigM){
    sample_discretized <- matrix(rmultinom(1, nvec[i], prob_vec), nrow = M)
    g_smpl <- generate_graph(t, sample_discretized)
    # Preparing the extended graph object built from a random sample of the distribution defined by prob_vec

    rcut_sample_xvst <- xcut_via_stmincut(g_smpl)
    rcut_sample_spec_clust <- manual_spec_clustering(g_smpl, method="kmeans")
    # Computes the spectral clustering of the graph g_smpl$graph
  
    rcut_sample[m, i+length(nvec)] <- ifelse( inherits(rcut_sample_spec_clust, "try-error"), NA,
                                          xcut_partition_S(which(rcut_sample_spec_clust$unnorm == 1), g_smpl)$rcut_min)
    rcut_sample[m, i] <- rcut_sample_xvst$rcut_min
    # Fills the matrix ncut_sample
    
    setTxtProgressBar(progress, m + BigM*(i-1))
    # Updates the progress bar
  }
}
# Now ncut_sample is filled with ncut values for different n using two methods: ncut_sample[,1:5] computed using Xvst,
# and ncut_sample[,6:10] used normalized spectral clustering.

close(progress)
# Closes the progress bar.

g_true <- generate_graph(t, matrix(prob_vec, nrow=M))
g_true_cut <- xcut_exact(g_true)
plot_graph_partition(g_true, g_true_cut$ncut, vertex_size = 30)

tm <- tidyr::pivot_longer(as.data.frame(rcut_sample), cols=c(4,9), names_to="cuts", values_to="values")
vlines_rcut <- data.frame(xint = c(g_true_cut$rcut_min,
                                   mean(rcut_sample[,4])/(n*(n-1)),
                                   mean(rcut_sample[,9])/(n*(n-1))), grp = letters[1:3])
# Prepare the data for plotting: tm contains the RCut samples for n=200 computed by the Xvst algorithm
# (column 4 of rcut_sample) and Unnormalized spectral clustering (column 9 of rcut_sample).
# vlines_rcut contains the true rcut value as well as the empirical mean of both samples.

ggplot(tm, aes(x = values/(n*(n-1)), y=(..count..)/sum(..count..), fill=cuts, group=cuts)) +
  geom_histogram(position="identity", bins=40, lty="blank", alpha=0.5) +
  geom_vline(data=vlines_rcut, aes(xintercept=xint, colour=grp), linetype="dashed", size=0.8) +
  labs(x="Ratio Cut value", y="Empirical density") +
  scale_x_continuous(breaks=c(seq(0.00005,0.00025,0.0001)), limits=c(0.00003,0.00035)) + theme_minimal() +
  scale_color_manual(name="Ratio Cut value",
                     labels=c("True Ratio Cut value", "Empirical mean Rvst", "Empirical mean spec. clust."),
                     values=c("red", "blue", "purple")) +
  scale_fill_manual(name="Ratio Cut bimodal distribution",
                    labels=c("RCut via st-MinCuts", "Unnorm. spectral clustering"),
                    values=c("blue", "purple")) +
  theme(legend.position = c(0.77,0.66), legend.box.background = element_rect(colour = "black")) +
  guides(fill=guide_legend(order=1), color=guide_legend(order=0))
# Comparing the Rvst algorithm and unnormalized spectral clustering. Notice the scaling of 'values'
# by a factor of n*(n-1) in order to make it unbiased.


# Where does spectral clustering fail? ------------------------------------

k <- 4
cockroach_king <- generate_cockroach(k)
cockroach_cut <- xcut_via_stmincut(cockroach_king)
# Create and cut a cockroach graph.

plot_graph_partition(cockroach_king, cockroach_cut$rcut,
                     edge_scale=10, xlim=c(0,2*k+1), ylim=c(0,3), main="RatioCut via st-MinCuts")
plot_graph_partition(cockroach_king, cockroach_cut$ncut,
                     edge_scale=10, xlim=c(0,2*k+1), ylim=c(0,3), main="NCut via st-MinCuts")
plot_graph_partition(cockroach_king, cockroach_cut$ccut,
                     edge_scale=10, xlim=c(0,2*k+1), ylim=c(0,3), main="CheegerCut via st-MinCuts")
# Xvst yields a reasonable partition (although the Cheeger Cut partition is not particularly good).
# Still, for all three cuts, Xvst agrees with the true XCut partition as can be seen below.

cockroach_cut_spectral <- manual_spec_clustering(cockroach_king, method = "naive")
plot_graph_partition(cockroach_king, which(cockroach_cut_spectral$unnorm==1),
                     edge_scale=10, xlim=c(0,2*k+1), ylim=c(0,3), main="Simple Spectral Clustering (Ratio Cut)")
plot_graph_partition(cockroach_king, which(cockroach_cut_spectral$norm==1),
                     edge_scale=10, xlim=c(0,2*k+1), ylim=c(0,3), main="Simple Spectral Clustering (NCut)")
# For comparison consider how naive spectral clustering fails. More advanced spectral clustering methods
# (e.g. using more than one eigenvector) lead to more complicated counterexamples, see Guattery, Miller (1993).

cockroach_cut_exact <- xcut_exact(cockroach_king)
plot_graph_partition(cockroach_king, cockroach_cut_exact$rcut,
                     edge_scale=10, xlim=c(0,2*k+1), ylim=c(0,3), main="True Ratio Cut partition")
plot_graph_partition(cockroach_king, cockroach_cut_exact$ncut,
                     edge_scale=10, xlim=c(0,2*k+1), ylim=c(0,3), main="True NCut partition")
plot_graph_partition(cockroach_king, cockroach_cut_exact$ccut,
                     xlim=c(0,2*k+1), ylim=c(0,3), main="True Cheeger Cut")
# Comparing this to the exact XCuts.


# QQ-plot of RCut empirical vs. asymptotic distribution -------------------

M <- 3
t <- 1
BigM <- 2000
eps <- 0.4
lmd <- sqrt(1+1.5*eps+eps^2)
prob_vec <- rep(1/M^2, M^2)
#prob_vec <- rep(c(lmd,eps,1), 3)
#prob_vec <- prob_vec/sum(prob_vec)
# The setting.

plist <- lapply(1:(2^(M^2-1)), function(k)(1:M^2)[intToBits(k)[1:M^2]==1])
# List of all partitions of any graph with M^2 vertices.

#nvec2 <- c(50, 100, 150, 200, 500)
nvec2 <- c(50, 200, 500, 1000, 10000)
# Consider a different vector of various n.

ncut_samples <- matrix(rep(NA,BigM*5), nrow=BigM)
progress <- txtProgressBar(min = 0, max = length(nvec2)*BigM, style = 3)
for(i in 1:length(nvec2)){
  for(m in 1:BigM){
    sample_discretized <- matrix(rmultinom(1, nvec2[i], prob_vec), nrow = M)
    g_sample <- generate_graph(t, sample_discretized)
    #nvst_sample <- xcut_via_stmincut(g_sample)
    nvst_sample <- xcut_exact(g_sample)
    ncut_samples[m, i] <- nvst_sample$ccut_min
    setTxtProgressBar(progress, m + BigM*(i-1))
  }
}
close(progress)
# Generating BigM iid. samples of the NCut statistic of a uniform distribution for M = 3 and t = 1.

xcut_exact_indexlist <- function(g, all.partitions=FALSE){
  #p_min_mincut <- p_min_ratiocut <- p_min_ncut <- p_min_cheegercut <- c()
  min_mincut <- Inf
  min_ncut <- Inf
  min_cheegercut <- Inf
  min_ratiocut <- Inf
  nvert <- length(V(g$graph))
  if(nvert > 25){
    print("Error: M has to be less than 6. You wanted a result in finite time, right?")
    return(NULL)
  }
  #for(k in 1:5){
  for(k in 1:(2^(nvert-1)-1)){
    p <- as.numeric(V(g$graph)$label)[intToBits(k)[1:nvert]==1]
    xcut_p <- xcut_partition_S(p, g)
    #print(xcut_p$rcut_min)
    #print(k)
    if(isTRUE(all.equal(xcut_p$mcut_min, min_mincut))) p_min_mincut <- c(p_min_mincut, k)
    if(isTRUE(all.equal(xcut_p$rcut_min, min_ratiocut))) p_min_ratiocut <- c(p_min_ratiocut, k)
    if(isTRUE(all.equal(xcut_p$ncut_min, min_ncut))) p_min_ncut <- c(p_min_ncut, k)
    if(isTRUE(all.equal(xcut_p$ccut_min, min_cheegercut))) p_min_cheegercut <- c(p_min_cheegercut, k)
    
    if(!isTRUE(all.equal(xcut_p$mcut_min, min_mincut)) & xcut_p$mcut_min < min_mincut){
      min_mincut <- xcut_p$mcut_min
      p_min_mincut <- k
    }
    if(!isTRUE(all.equal(xcut_p$rcut_min, min_ratiocut)) & xcut_p$rcut_min < min_ratiocut){
      min_ratiocut <- xcut_p$rcut_min
      p_min_ratiocut <- k
    }
    if(!isTRUE(all.equal(xcut_p$ncut_min, min_ncut)) & xcut_p$ncut_min < min_ncut){
      #print("YATTA")
      #print(xcut_p$ncut_min/min_ncut - 1)
      min_ncut <- xcut_p$ncut_min
      p_min_ncut <- k
    }
    if(!isTRUE(all.equal(xcut_p$ccut_min, min_cheegercut)) & xcut_p$ccut_min < min_cheegercut){
      min_cheegercut <- xcut_p$ccut_min
      p_min_cheegercut <- k
    }
    #print(k)
    #print(xcut_p$ncut_min)
    #print(min_ncut)
    #print(p_min_ncut)
  }
  if(!all.partitions) return(list(mcut=p_min_mincut, rcut=p_min_ratiocut, ncut=p_min_ncut, ccut=p_min_cheegercut, 
              mcut_min=min_mincut, rcut_min=min_ratiocut, ncut_min=min_ncut, ccut_min=min_cheegercut))
  add_cmp <- function(cutvec) c(cutvec, rev(2^nvert - 1 - cutvec))
  return(list(mcut=add_cmp(p_min_mincut), rcut=add_cmp(p_min_ratiocut), ncut=add_cmp(p_min_ncut), ccut=add_cmp(p_min_cheegercut), 
              mcut_min=min_mincut, rcut_min=min_ratiocut, ncut_min=min_ncut, ccut_min=min_cheegercut))
  # The returned partitions (in $rcut, $ncut or $ccut) are given by a single number k that can be decoded by
  # transforming it into binary. For a more detailed explanation of the code see xcut_exact.
}
# Computes indeces in plist of all partitions that attain the minimal XCut. Note that one can encode a partition
# (e.g. {1,2,5} of {1,...,5}) as an integer by using its binary representation (i.e. {1,2,5} -> {1,1,0,0,1} -> 19).
# Input: Extended graph object g. Output: Vector of integers that each represent a partition.

g_true_uniform <- generate_graph(t, matrix(prob_vec, nrow=M))
g_true_uniform_cut <- xcut_exact(g_true_uniform)
# Generate and cut the true underlying graph.

Sigma_p_ncut <- outer(plist, plist, Vectorize(function(p,q) xcut_cov(prob_vec, p, g_true_uniform, q, g_true_uniform)$ccut_cov))
# Computing the covariance matrix of the limiting distribution of the NCut statistic
# sqrt(n) * (\widehat{NC}(G_{n,t,m}) - NC(P_{t,m})).

g_true_uniform_cut_indexlist <- xcut_exact_indexlist(g_true_uniform)
# Determining all partitions that attain the XCuts.

qq_line_params <- function(sample, theoretical){
  y <- quantile(sample[!is.na(sample)], c(0.25, 0.75))
  x <- quantile(theoretical[!is.na(theoretical)], c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1] - slope * x[1]
  return(as.numeric(c(slope, int)))
}
# Manual calculation of the parameters of a quantile-quantile line going through the 0.25 and the 0.75 quantile of
# 'sample' against (a list of quantiles or a sample of) the theoretical distribution 'theoretical'.

ncut_true_sample2 <- replicate(2000, min(mvtnorm::rmvnorm(1, sigma =
                            Sigma_p_ncut[g_true_uniform_cut_indexlist$ccut, g_true_uniform_cut_indexlist$ccut])))
# Generating a sample from the asymptotic distribution of the XCut statistic
ncut_true_sample <- rnorm(2000, 0, sqrt(Sigma_p_ncut[g_true_uniform_cut_indexlist$rcut, g_true_uniform_cut_indexlist$rcut]))

# instead: ncut_samples/rep(nvec2*(nvec2-1),each=2000)
ncut_samples_new <- (ncut_samples)[sapply(1:2000, function(x)!anyNA(ncut_samples[x,])),]
ncs_new <- apply(rep(sqrt(nvec2),each=2000) * (ncut_samples_new[,1:5] - g_true_uniform_cut$ccut_min),2,sort,decreasing=F)
# Compute the statistic sqrt(n) * (\widehat{NC}(G_{n,t,m}) - NC(P_{t,m})).

ncs_pivot <- tidyr::pivot_longer(as.data.frame(ncs_new), cols=1:5, names_to="cuts", values_to="values")
ncnormals <- qnorm(seq(0,1,length.out=2002), sd=sqrt(xcut_cov(prob_vec, c(1,2,3), g_true_uniform, c(1,2,3), g_true_uniform)$ccut_cov))[-c(1,2002)]
qq_abl <- qq_line_params(ncnormals, sort(ncut_true_sample))

qq_abl2 <- qq_line_params(ncut_false_sample, ncut_true_sample)
qq_abl <- qq_line_params(ncut_true_sample2, ncut_true_sample)
# Quantile-quantile line connecting the respective 0.25 and 0.75 quantiles of the sample ncnormals and (a sample from)
# the true limiting distribution ncut_true_sample.

stdnormals2 <- qnorm(seq(0,1,length.out=2002))[-c(1,2002)]
# Quantiles of standard normals that would constitute the limiting distribution if the XCut minimizing partition
# were unique. This is used to demonstrate that knowing whether uniqueness holds is critical.

ncs_plus_theoretical <- cbind(ncs_pivot, Theoretical=rep(sort(ncut_true_sample), each=5))
# Combining sample and theoretical distribution.

#gg1 <- ggplot(ncs_plus_theoretical) + geom_point(aes(x=Theoretical, y=values, colour=cuts), key_glyph="point") +
#  geom_abline(mapping=aes(slope=1,intercept=0, colour="black")) +
#  geom_abline(mapping=aes(slope=qq_abl[1],intercept=qq_abl[2], colour="red"), linetype="dashed") +
#  labs(x="Theoretical", y="Sample") + theme_minimal() +
#  scale_color_manual(name="Normalized Cut", values = c("red", "black", "green4", "lawngreen", "yellow", "darkgoldenrod2", "chocolate1"), labels=c("Asymptotic distribution", "Assuming uniqueness", "n = 50", "n = 200", "n = 500", "n = 1000", "n = 10000"),
#                     guide = guide_legend(override.aes = list(linetype = c("solid", "dashed", rep("blank", 5)), shape = c(NA, NA, rep(16, 5))))) +
#  theme(legend.position = c(0.18,0.75), legend.box.background = element_rect(colour = "black")) #+
# Plotting the aforementioned statistic and its asymptotic limiting distribution for various n.
gg1 <- ggplot(ncs_plus_theoretical) + geom_point(aes(x=Theoretical, y=values, colour=cuts), key_glyph="point") +
  geom_abline(mapping=aes(slope=1,intercept=0, colour="black")) +
  geom_abline(mapping=aes(slope=qq_abl[1],intercept=qq_abl[2], colour="purple"), linetype="dashed") +
  geom_abline(mapping=aes(slope=qq_abl2[1],intercept=qq_abl2[2], colour="red"), linetype="dashed") +
  labs(x="Theoretical", y="Sample") + theme_minimal() +
  scale_color_manual(name="Cheeger Cut", values = c("red", "purple", "black", "green4", "lawngreen", "yellow", "darkgoldenrod2", "chocolate1"), labels=c("Asymptotic distribution", "Assuming unequal vol.", "Assuming uniqueness", "n = 50", "n = 200", "n = 500", "n = 1000", "n = 10000"),
                     guide = guide_legend(override.aes = list(linetype = c("solid", "dashed", "dashed", rep("blank", 5)), shape = c(NA, NA, NA, rep(16, 5))))) +
  theme(legend.position = c(0.18,0.72), legend.box.background = element_rect(colour = "black"))
gg1
ggsave("qq_uniform_ccut_overn_1_correct.svg", plot=gg1, width=6, height=4.5)

check_for_equal_vol <- function(p1, S, gS){
  if(is.null(S)) return(FALSE)
  p <- as.vector(t(matrix(p1, sqrt(length(p1)))))
  S_C <- (1:length(p))[-S]
  q_S <- sapply(1:length(p), function(x)sum(p[S[gS$weight_mat[S,x]>0 | gS$weight_mat[x,S]>0]]))
  q_S_C <- sapply(1:length(p), function(x)sum(p[S_C[gS$weight_mat[S_C,x]>0 | gS$weight_mat[x,S_C]>0]]))
  vol_S <- sum(p*q_S)
  vol_S_C <- sum(p*q_S_C)
  return(isTRUE(all.equal(vol_S, vol_S_C)))
  # Returns TRUE if both volumes are equal (BEWARE: We have to compensate for floating-point imprecision here) and FALSE else
}

generate_ccut_limiting_dist_sample_equal_vol <- function(p1, S, gS){
  if(is.null(S)) return(NA)
  p <- as.vector(t(matrix(p1, sqrt(length(p1)))))
  Sigma_mult <- diag(p) - p%*%t(p)
  Z <- as.vector(mvtnorm::rmvnorm(1, sigma=Sigma_mult))
  S_C <- (1:length(p))[-S]
  q_S <- sapply(1:length(p), function(x)sum(p[S[gS$weight_mat[S,x]>0 | gS$weight_mat[x,S]>0]]))
  q_S_C <- sapply(1:length(p), function(x)sum(p[S_C[gS$weight_mat[S_C,x]>0 | gS$weight_mat[x,S_C]>0]]))
  aZ <- sapply(1:length(p), function(x)sum(Z[gS$weight_mat[,x]>0 | gS$weight_mat[x,]>0]))
  if(sum(Z[S]*(q_S[S]+q_S_C[S])) + sum(p[S]*aZ[S]) <= sum(Z[S_C]*(q_S[S_C]+q_S_C[S_C])) + sum(p[S_C]*aZ[S_C])){
    vol_S_min <- sum(p*q_S)
    S_min <- S
    S_C_min <- S_C
    q_S_min <- q_S
    q_S_C_min <- q_S_C
  }
  else{
    vol_S_min <- sum(p*q_S_C)
    S_min <- S_C
    S_C_min <- S
    q_S_min <- q_S_C
    q_S_C_min <- q_S
  }
  g_cc_S <- xcut_partition_S(S, gS)$ccut
  ccut_coeff <- (sum(Z[S_min]*(q_S_C_min[S_min] - g_cc_S*(q_S_C_min[S_min]+2*q_S_min[S_min])))
                 + sum(Z[S_C_min]*(q_S_min[S_C_min] - g_cc_S*q_S_min[S_C_min])))/vol_S_min
  return(ccut_coeff)
}

generate_ccut_limiting_dist_sample_equal_vol_S_list <- function(p1, Slist, gS){
  #if(is.null(S)) return(NA)
  ret_vec <- c()
  p <- as.vector(t(matrix(p1, sqrt(length(p1)))))
  Sigma_mult <- diag(p) - p%*%t(p)
  Z <- as.vector(mvtnorm::rmvnorm(1, sigma=Sigma_mult))
  for(S in Slist){
    S_C <- (1:length(p))[-S]
    q_S <- sapply(1:length(p), function(x)sum(p[S[gS$weight_mat[S,x]>0 | gS$weight_mat[x,S]>0]]))
    q_S_C <- sapply(1:length(p), function(x)sum(p[S_C[gS$weight_mat[S_C,x]>0 | gS$weight_mat[x,S_C]>0]]))
    aZ <- sapply(1:length(p), function(x)sum(Z[gS$weight_mat[,x]>0 | gS$weight_mat[x,]>0]))
    if(sum(Z[S]*(q_S[S]+q_S_C[S])) + sum(p[S]*aZ[S]) <= sum(Z[S_C]*(q_S[S_C]+q_S_C[S_C])) + sum(p[S_C]*aZ[S_C])){
      vol_S_min <- sum(p*q_S)
      S_min <- S
      S_C_min <- S_C
      q_S_min <- q_S
      q_S_C_min <- q_S_C
    }
    else{
      vol_S_min <- sum(p*q_S_C)
      S_min <- S_C
      S_C_min <- S
      q_S_min <- q_S_C
      q_S_C_min <- q_S
    }
    g_cc_S <- xcut_partition_S(S, gS)$ccut
    ccut_coeff <- (sum(Z[S_min]*(q_S_C_min[S_min] - g_cc_S*(q_S_C_min[S_min]+2*q_S_min[S_min])))
                   + sum(Z[S_C_min]*(q_S_min[S_C_min] - g_cc_S*q_S_min[S_C_min])))/vol_S_min
    ret_vec <- c(ret_vec, ccut_coeff)
  }
  return(ret_vec)
}

generate_ccut_limiting_dist_sample_list <- function(p1, Slist, gS){
  #if(is.null(S)) return(NA)
  ret_vec <- c()
  p <- as.vector(t(matrix(p1, sqrt(length(p1)))))
  Sigma_mult <- diag(p) - p%*%t(p)
  Z <- as.vector(mvtnorm::rmvnorm(1, sigma=Sigma_mult))
  for(S in Slist){
    S_C <- (1:length(p))[-S]
    q_S <- sapply(1:length(p), function(x)sum(p[S[gS$weight_mat[S,x]>0 | gS$weight_mat[x,S]>0]]))
    q_S_C <- sapply(1:length(p), function(x)sum(p[S_C[gS$weight_mat[S_C,x]>0 | gS$weight_mat[x,S_C]>0]]))
    #print(q_S)
    #print(q_S_C)
    vol_S <- sum(p*q_S)
    vol_S_C <- sum(p*q_S_C)
    aZ <- sapply(1:length(p), function(x)sum(Z[gS$weight_mat[,x]>0 | gS$weight_mat[x,]>0]))
    if((!isTRUE(all.equal(vol_S,vol_S_C)) & vol_S < vol_S_C) |
       (isTRUE(all.equal(vol_S,vol_S_C)) & sum(Z[S]*(q_S[S]+q_S_C[S])) + sum(p[S]*aZ[S]) <=
                                          sum(Z[S_C]*(q_S[S_C]+q_S_C[S_C])) + sum(p[S_C]*aZ[S_C]))){
      vol_S_min <- vol_S
      S_min <- S
      S_C_min <- S_C
      q_S_min <- q_S
      q_S_C_min <- q_S_C
    }
    else{
      vol_S_min <- vol_S_C
      S_min <- S_C
      S_C_min <- S
      q_S_min <- q_S_C
      q_S_C_min <- q_S
    }
    g_cc_S <- xcut_partition_S(S, gS)$ccut_min
    #print(paste("g_cc_S:",g_cc_S))
    #print(paste("S_min:",S_min))
    ccut_coeff <- (sum(Z[S_min]*(q_S_C_min[S_min] - g_cc_S*(q_S_C_min[S_min]+2*q_S_min[S_min])))
                   + sum(Z[S_C_min]*(q_S_min[S_C_min] - g_cc_S*q_S_min[S_C_min])))/vol_S_min
    ret_vec <- c(ret_vec, ccut_coeff)
    #print(q_S_min)
    #print(q_S_C_min)
    #print(paste("ret_vec:",ret_vec))
    #print(Z[S_min])
    #print(Z[S_C_min])
    #print(Z[S_min]*(q_S_C_min[S_min] - g_cc_S*(q_S_C_min[S_min]+2*q_S_min[S_min])))
    #print(Z[S_C_min]*(q_S_min[S_C_min] - g_cc_S*q_S_min[S_C_min]))
  }
  return(ret_vec)
}


ccut_uniform_eqvol_sample <- replicate(2000, generate_ccut_limiting_dist_sample_equal_vol_S_list(prob_vec, list(c(1,2,4,5),c(2,3,5,6),c(4,5,7,8),c(1,2,3,4,7)), g5_true))
ccut_true_uniform <- apply(ccut_uniform_eqvol_sample, 2, min)
ccut_uniform_eqvol_sample2 <- sapply(list(c(1,2,4,5),c(2,3,5,6),c(4,5,7,8),c(1,2,3,4,7)), function(S)replicate(2000, generate_ccut_limiting_dist_sample_equal_vol(prob_vec, S, g5_true)), simplify="array")
ccut_true_uniform2 <- apply(ccut_uniform_eqvol_sample2, 1, min)

ncut_false_sample <- replicate(2000, generate_ccut_limiting_dist_sample_equal_vol(prob_vec, c(1,2,4,5), g_true_uniform))
ccut_uniform_eqvol_sample2 <- sapply(list(c(1,2,4,5),c(2,3,5,6),c(4,5,7,8),c(1,2,3,4,7)), function(S)replicate(2000, generate_ccut_limiting_dist_sample_equal_vol(prob_vec, S, g_true_uniform)), simplify="array")
ncut_true_sample <- apply(ccut_uniform_eqvol_sample2, 1, min)
Sigma_p_ncut[lower.tri(Sigma_p_ncut)] <- t(Sigma_p_ncut)[lower.tri(Sigma_p_ncut)]
#list(c(1,2,4,5),c(2,3,5,6),c(4,5,7,8),c(5,6,8,9))

ncs_plus_theoretical2 <- ncs_plus_theoretical[ncs_plus_theoretical$cuts=="V1"|ncs_plus_theoretical$cuts=="V3"|ncs_plus_theoretical$cuts=="V5",]
gg1s <- ggplot(ncs_plus_theoretical2) + geom_point(aes(x=Theoretical, y=values, colour=cuts), key_glyph="point") +
  geom_abline(mapping=aes(slope=1,intercept=0, colour="black")) +
  geom_abline(mapping=aes(slope=qq_abl[1],intercept=qq_abl[2], colour="purple"), linetype="dashed") +
  geom_abline(mapping=aes(slope=qq_abl2[1],intercept=qq_abl2[2], colour="red"), linetype="dashed") +
  labs(x="Theoretical", y="Sample") + theme_minimal() +
  scale_color_manual(name="Cheeger Cut", values = c("red", "purple", "black", "green4", "yellow", "chocolate1"), labels=c("Asymptotic distribution", "Assuming unequal vol.", "Assuming uniqueness", "n = 50", "n = 500", "n = 10000"),
                     guide = guide_legend(override.aes = list(linetype = c("solid", "dashed", "dashed", rep("blank", 3)), shape = c(NA, NA, NA, rep(16, 3))))) +
  theme(legend.position = c(0.18,0.72), legend.box.background = element_rect(colour = "black"))
gg1s



# Comparing different choices of t using Kolmogorov-Smirnov ---------------

M <- 3
eps <- 0.4
prob_vec <- rep(c(sqrt(1+eps^2+1.5*eps), eps, 1), each=3)
prob_vec <- prob_vec / sum(prob_vec)
prob_vec <- rep(1/M^2, M^2)
BigM <- 2000
# The setting.

plist <- lapply(1:(2^(M^2-1)-1), function(k)(1:M^2)[intToBits(k)[1:M^2]==1])
#plist <- lapply(1:(2^(M^2)-2), function(k)(1:M^2)[intToBits(k)[1:M^2]==1]) # to account for complements also
t_vec <- sort(unique(unlist(sapply(1:(M-1), function(x)sqrt(x^2+(0:x)^2)))))
g_prob_list <- lapply(t_vec, function(t1)generate_graph(t1, matrix(prob_vec, nrow=M)))
g_prob_cut_list <- lapply(g_prob_list, xcut_exact)
g_prob_cut_index_list <- lapply(g_prob_list, function(g)xcut_exact_indexlist(g)$ccut)
# For each t in t_vec, the vector of all t that yield new connections between vertices, build the corresponding graph
# and compute its NCut value (g_prob_cut_list) as well as the partitions that attain that value (g_prob_cut_index_list).

progress <- txtProgressBar(min = 0, max = length(t_vec)*length(plist), style = 3)
Sigma_p_ncut_component <- matrix(rep(NA, (2^(M^2-1)-1)^2), 2^(M^2-1)-1)
#Sigma_p_ncut_component <- matrix(rep(NA, (2^(M^2)-2)^2), 2^(M^2)-2) # to account for complements also
Sigma_p_ncut_list <- list()
for(t1 in 1:length(t_vec)){
  for(p in 1:length(plist)){
    for(q in p:length(plist)){
      spn_term <- xcut_cov(prob_vec, plist[[p]], g_prob_list[[t1]], plist[[q]], g_prob_list[[t1]])
      Sigma_p_ncut_component[p,q] <- spn_term$ccut_cov
    }
    setTxtProgressBar(progress, p + length(plist)*(t1-1))
  }
  Sigma_p_ncut_component[lower.tri(Sigma_p_ncut_component)] <-
                                                t(Sigma_p_ncut_component)[lower.tri(Sigma_p_ncut_component)]
  Sigma_p_ncut_list[[t1]] <- Sigma_p_ncut_component
}
close(progress)
# For each t, compute the asymptotic covariance matrix Sigma_p_ncut_component, and store them in Sigma_p_ncut_list

#ncut_true_sample_list_over_t_2 <- c()
#for(i in 1:50){
ncut_true_sample_list_over_t <- lapply(2:length(t_vec), function(t1)replicate(100000,
                                            min(mvtnorm::rmvnorm(1, sigma=Sigma_p_ncut_list[[t1]])[g_prob_cut_index_list[[t1]]])))
# For each t, generate a sample of the true limiting distribution of the NCut statistic.
#ncut_true_sample_list_over_t_2 <- c(ncut_true_sample_list_over_t_2, ncut_true_sample_list_over_t_a)
#print(i)
#}
ncut_true_sample_list_over_t <- lapply(1:length(t_vec), function(t1)replicate(50000, min(mvtnorm::rmvnorm(1, sigma=matrix((Sigma_p_ncut_list[[t1]])[g_prob_cut_index_list[[t1]], g_prob_cut_index_list[[t1]]], nrow=length(g_prob_cut_index_list[[t1]]))))))

ccut_uniform_eqvol_sample <- replicate(50000, generate_ccut_limiting_dist_sample_equal_vol_S_list(prob_vec, list(c(1,2,4,5),c(2,3,5,6),c(4,5,7,8),c(1,2,3,4,7)), g_prob_list[[1]]))
ccut_uniform_eqvol_sample <- replicate(50000, generate_ccut_limiting_dist_sample_equal_vol_S_list(prob_vec, list(c(1,4,7)), g_prob_list[[1]]))

ccut_true_uniform <- apply(ccut_uniform_eqvol_sample, 2, min)

ccut_dummy <- lapply(1:length(t_vec), function(t1)replicate(1000, min(mvtnorm::rmvnorm(1, sigma=Sigma_p_ncut_list[[t1]])[c(3,6)])))
ncut_true_sample_list_over_t <- lapply(1:length(t_vec), function(t1)replicate(50000, generate_ccut_limiting_dist_sample_list(prob_vec, lapply(g_prob_cut_index_list[[t1]], function(k)(1:M^2)[intToBits(k)[1:M^2]==1]), g_prob_list[[t1]])))

#save(ncut_true_sample_list_over_t, file="ccut_uniform_true_100000.Rdata")
load("ccut_uniform_true_50000_2.Rdata")
ncut_true_sample_list_over_t <- ncut_true_sample_list_over_t_2
rm(ncut_true_sample_list_over_t_2)
gc()
#nvec3 <- c(500, 1000, 5000, 10000, 50000)
nvec3 <- c(500, 1000, 2000, 3000, seq(4000,20000,2000))
nvec3 <- c(1000, seq(2000,50000,2000))
# Defining yet another vector of sample sizes n.
nvec3 <- c(5000, 10000, 50000, seq(1e5, 1e6, 1e5))
nvec3 <- seq(1e9, 1e10, 1e9)
nvec3 <- c(1e7, 1e8, 4e8, 7e8, seq(1e9,1e10,1e9))
nvec3 <- 10^(4:13)
nvec3 <- c(1e7, 1e8, 5e8, 1e9, 3e9, 5e9, seq(1e10,1e11,1e10))
nvec3 <- 10000*c(1,2,5,seq(10,100,10))

progress <- txtProgressBar(min = 0, max = length(t_vec)*length(nvec3)*BigM, style = 3)
ncut_samples_over_t <- matrix(rep(NA,BigM*length(t_vec)*length(nvec3)), nrow=BigM)
for(t1 in 1:length(t_vec)){
  for(i in 1:length(nvec3)){
    for(m in 1:BigM){
      sample_discretized <- matrix(rmultinom(1, nvec3[i], prob_vec), nrow = M)
      g_sample <- generate_graph(t_vec[t1], sample_discretized)
      ncut_samples_over_t[m, i+length(nvec3)*(t1-1)] <- xcut_exact(g_sample)$ccut_min
      # Computes the exact NCut of g_sample and store it as one sample value in ncut_samples_over_t.
      
      #ncut_samples_over_t[m, i+5*(t1-1)] <- xcut_via_stmincut(g_sample)$ncut_min
      # Alternatively, use the Xvst algorithm. The results are very similar, though.
      
      setTxtProgressBar(progress, m + BigM*(i-1) + BigM*length(nvec3)*(t1-1))
    }
  }
}
close(progress)
# Generates BigM samples of the NCut value for various n and t.

save(ncut_samples_over_t, file="ccut_uniform_samples_over_t_100k.Rdata")
#load("ccut_uniform_samples_over_t.Rdata")

# HIER EIGENTLICH: (nvec3[i]^2 * ncut_samples_over_t[,...]
nc_ks_cbimodal <- sapply(1:length(nvec3), function(i)sapply(1:length(t_vec), function(t1)ks.test(sqrt(nvec3[i]) *
                    (ncut_samples_over_t[, i + length(nvec3) * (t1-1)] - xcut_exact(generate_graph(t1, matrix(all_samples_discretized[,m+(i-1)*BigM+(t1-1)*BigM*length(nvec3)]/nvec3[i], nrow=M)))$ccut_min),#g_prob_cut_list[[t1]]$ccut_min),
                    ncut_true_sample_list_over_t[[t1]])$statistic), simplify="array")
# Computing the Kolmogorov-Smirnov statistic of the empirical distribution of the NCut statistic against its theoretical
# distribution for n from nvec3 and over all t (in t_vec). Please ignore the warnings.

nc_ks_cbimodal_df <- tidyr::pivot_longer(as.data.frame(rbind(rep(0,5),nc_ks_cbimodal[,1:5],nc_ks_cbimodal[5,]), row.names=NULL),
                                         cols=1:5, names_to="cuts", values_to="values")
nc_ks_cbimodal_df_new <- cbind(nc_ks_cbimodal_df, x=rep(c(0,t_vec,4), each=5))
nc_ks_cbimodal_df_new$xend <- rep(c(t_vec,4,NA), each=5)
# Transforming the data and determining the coordinates of the illustratory dots for the plot below.

gg2 <- ggplot(nc_ks_cbimodal_df_new, aes(x=x, y=values, xend=xend, yend=values, colour=cuts)) +
  geom_vline(aes(xintercept=x), linetype=2, color="grey") + ylim(c(0,0.8)) +
  geom_hline(yintercept = 0.042947, colour="blue", linetype="dashed", show.legend=FALSE) +
  geom_point() + geom_point(aes(x=xend, y=values), shape=1) + geom_segment(key_glyph="path") +
  coord_cartesian(xlim=c(0,3.7)) + labs(x="t", y="Kolmogorov-Smirnov statistic") + theme_minimal() +
  scale_color_manual(name="Cheeger Cut", values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                     labels=c("n=500", "n=1000", "n=5000", "n=10000", "n=50000"),
                     guide = guide_legend(override.aes = list(linetype = rep("solid", 5), shape = rep(NA, 5)))) +
  theme(legend.position = c(0.16,0.7), legend.box.background = element_rect(colour = "black"))
# Plotting the above Kolmogorov-Smirnov statistic against t for different n.

gg2
ggsave("ks_uniform_ccut.svg", plot=gg2, width=6, height=4.5)

### NEW ONE:

manual_KSd <- function(n, m=50000, signif=0.05) sqrt(-log(signif/2)/2 * (n+m)/(n*m))
row.names(nc_ks_cbimodal) <- NULL
nc_ks_cbimodal_df_2 <- tidyr::pivot_longer(as.data.frame(as.matrix(t(nc_ks_cbimodal))), cols=1:5, names_to="cuts", values_to="values")
nc_ks_cbimodal_df_2_new <- cbind(nc_ks_cbimodal_df_2, x=rep(nvec3, each=5), Dks=rep(manual_KSd(nvec3), each=5))
#nc_ks_cbimodal_df_2_new$xend <- rep(c(t_vec,4,NA), each=5)
# Transforming the data and determining the coordinates of the illustratory dots for the plot below.
library(latex2exp)
#library(sfsmisc)

gg3 <- ggplot(nc_ks_cbimodal_df_2_new, aes(x=x, y=values, colour=cuts)) +
  #geom_vline(aes(xintercept=x), linetype=2, color="grey") + ylim(c(0,0.8)) +
  geom_line(aes(x=x, y=values, color=cuts), size=1.2) + #ylim(c(0,0.2)) +
  geom_ribbon(aes(x=x, y=values, ymin=pmax(values-Dks,0), ymax=pmin(values+Dks,1), fill=cuts), alpha=0.2, linetype=0) +
  geom_hline(yintercept = 0.03096948, colour="blue", linetype="dashed", show.legend=FALSE) +
  #geom_point() + geom_point(aes(x=xend, y=values), shape=1) + geom_segment(key_glyph="path") +
  coord_cartesian(xlim=c(0,1e10), ylim=c(0,.85)) + labs(x="n", y="Kolmogorov-Smirnov statistic") + theme_minimal() +
  scale_color_manual(name=TeX("Cheeger Cut bootstrap for $M\\,=\\,n$"), values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                     labels=unname(TeX(c("$1\\leq t \\,< \\,\\sqrt{2}$", "$\\sqrt{2}\\leq t \\,<\\;2\\,$", "$2\\leq t \\,< \\,\\sqrt{5}$",
                                         "$\\sqrt{5}\\leq t \\,< \\,\\sqrt{8}$", "$\\sqrt{8} \\leq t \\;\\;\\;\\;$"))),
                     guide = guide_legend(nrow=2, byrow=TRUE, override.aes = list(linetype = rep("solid", 5), shape = rep(NA, 5)))) +
  scale_fill_manual(name=TeX("Cheeger Cut bootstrap for $M\\,=\\,n$"), values = c("green4", "lawngreen", "yellow", "darkgoldenrod2", "orangered1"),
                    labels=unname(TeX(c("$1\\leq t \\,< \\,\\sqrt{2}$", "$\\sqrt{2}\\leq t \\,<\\;2\\,$", "$2\\leq t \\,< \\,\\sqrt{5}$",
                                        "$\\sqrt{5}\\leq t \\,< \\,\\sqrt{8}$", "$\\sqrt{8} \\leq t \\;\\;\\;\\;$")))) +
  theme(legend.position = c(0.65,0.85), legend.box.background = element_rect(fill="white",colour = "black")) + #guides(linetype = guide_legend(nrow = 2))
  scale_x_continuous(breaks=1e9*seq(0,10,2.5),
                     labels=unname(TeX(c("$0$","$10^9$","$2.5\\cdot 10^9$","$5\\cdot 10^9$","$10^{10}$")))) #trans_format("log10", math_format(.x^2)))
gg3
#ggsave("ks_uniform_ccut_overn_dampened2.svg", plot=gg3, width=6, height=4.5)


####---- Second testing range? ----####

progress <- txtProgressBar(min = 0, max = length(t_vec)*length(nvec3)*BigM, style = 3)
ncut_samples_over_t_2 <- matrix(rep(NA,BigM*length(t_vec)*length(nvec3)), nrow=BigM)
for(t1 in 1:length(t_vec)){
  for(i in 1:length(nvec3)){
    for(m in 1:BigM){
      sample_discretized <- matrix(rmultinom(1, nvec3[i], prob_vec), nrow = M)
      g_sample <- generate_graph(t_vec[t1], sample_discretized)
      ncut_samples_over_t_2[m, i+length(nvec3)*(t1-1)] <- xcut_exact(g_sample)$ccut_min
      # Computes the exact NCut of g_sample and store it as one sample value in ncut_samples_over_t.
      
      #ncut_samples_over_t[m, i+5*(t1-1)] <- xcut_via_stmincut(g_sample)$ncut_min
      # Alternatively, use the Xvst algorithm. The results are very similar, though.
      
      setTxtProgressBar(progress, m + BigM*(i-1) + BigM*length(nvec3)*(t1-1))
    }
  }
}
close(progress)
# Generates BigM samples of the NCut value for various n and t.

load("all_samples_discretized_m3_BigM2000_n26_t5.Rdata")

progress <- txtProgressBar(min = 0, max = length(t_vec)*length(nvec3)*BigM, style = 3)
#ncut_bootstrap_samples_over_t_3 <- matrix(rep(NA,BigM*length(t_vec)*length(nvec3)), nrow=BigM)
all_samples_discretized <- matrix(rep(NA,M^2*BigM*length(nvec3)*length(t_vec)), nrow=M^2)
for(t1 in 1:length(t_vec)){
  for(i in 1:length(nvec3)){
    for(m in 1:BigM){
      all_samples_discretized[,m+(i-1)*BigM+(t1-1)*BigM*length(nvec3)] <- rmultinom2(1, nvec3[i], prob_vec)
#      asd_ccuts[m+(i-1)*BigM+(t1-1)*BigM*length(nvec3)] <- xcut_exact(generate_graph(t1, matrix(all_samples_discretized[,m+(i-1)*BigM+(t1-1)*BigM*length(nvec3)]/nvec3[i], nrow=M)))$ccut_min
#      sample_discretized <- matrix(rmultinom(1, round(nvec3[i]^(4/5)), all_samples_discretized[,m+(i-1)*BigM+(t1-1)*BigM*length(nvec3)]/nvec3[i]), nrow = M)
#      second_sample_discretized <- matrix(rmultinom(1, nvec3[i], all_samples_discretized[,m+(i-1)*BigM+(t1-1)*BigM*length(nvec3)]/nvec3[i]), nrow = M)
#      g_sample <- generate_graph(t_vec[t1], sample_discretized + second_sample_discretized / sqrt(round(nvec3[i]^(4/5))))
#      ncut_bootstrap_samples_over_t_3[m, i+length(nvec3)*(t1-1)] <- xcut_exact(g_sample)$ccut_min
      # Computes the exact NCut of g_sample and store it as one sample value in ncut_samples_over_t.
      
      #ncut_samples_over_t[m, i+5*(t1-1)] <- xcut_via_stmincut(g_sample)$ncut_min
      # Alternatively, use the Xvst algorithm. The results are very similar, though.
      
      setTxtProgressBar(progress, m + BigM*(i-1) + BigM*length(nvec3)*(t1-1))
    }
  }
}
close(progress)
save(all_samples_discretized, file="all_samples_discretized_n1e4-13_BigM100.Rdata")

ncut_samples_over_t <- ncut_bootstrap_samples_over_t_3
save(ncut_bootstrap_samples_over_t_3, file="ccut_bootstrap2_samples_meq100logn.Rdata")

sapply(g_prob_cut_index_list[[3]], function(S)check_for_equal_vol(prob_vec, plist[[S]], g_prob_list[[3]]))

nc_ks_cbimodal <- sapply(1:length(nvec3), function(i)sapply(1:length(t_vec), function(t1)ks.test(sqrt(round(nvec3[i]^(4/5))) *
                             (ncut_samples_over_t[, i + length(nvec3) * (t1-1)] - rep(asd_ccuts, each=BigK)),#g_prob_cut_list[[t1]]$ccut_min),
                              ncut_true_sample_list_over_t[[t1]])$statistic), simplify="array")

BigM <- 200
BigK <- 200
test_ccuts_bootstrap3 <- matrix(rep(NA,BigM*BigK*length(t_vec)*length(nvec3)), nrow=BigM*BigK)
test_ccuts_bootstrap4 <- matrix(rep(NA,BigM*BigK*length(t_vec)*length(nvec3)), nrow=BigM*BigK)
progress <- txtProgressBar(min = 0, max = length(t_vec)*length(nvec3)*BigM, style = 3)
#ncut_bootstrap_samples_over_t_3 <- matrix(rep(NA,BigM*length(t_vec)*length(nvec3)), nrow=BigM)
for(t1 in 1:length(t_vec)){
  for(i in 1:length(nvec3)){
    for(m in 1:BigM){
      smallM <- round(nvec3[i]^(4/5))#round(50*sqrt(nvec3[i]))#round(nvec3[i]^(4/5))
      current_sample <- all_samples_discretized[,m+(i-1)*BigM+(t1-1)*BigM*length(nvec3)]/nvec3[i]
      for(k in 1:BigK){
        sample_discretized <- matrix(rmultinom2(1, smallM, current_sample), nrow = M)/smallM
        #test_ccuts_bootstrap[k+(m-1)*BigK,] <- xcut_exact(generate_graph(t_vec[3], sample_discretized))$ccut_min
        #print(xcut_exact(generate_graph(t_vec[3], current_sample + sample_discretized*sqrt(smallM)))$ccut)
        test_ccuts_bootstrap3[k+(m-1)*BigK, i+(t1-1)*length(nvec3)] <- xcut_exact(generate_graph(t_vec[t1], matrix(current_sample,nrow=M) + sample_discretized*sqrt(smallM)))$ccut_min
      }
      #test_ccuts_bootstrap2[(1:BigK)+(m-1)*BigK,] <- sqrt(smallM) * (test_ccuts_bootstrap[(1:BigK)+(m-1)*BigK,] - rep(xcut_exact(generate_graph(t_vec[3], matrix(current_sample, nrow=M)))$ccut_min, BigK))
      test_ccuts_bootstrap4[(1:BigK)+(m-1)*BigK, i+(t1-1)*length(nvec3)] <- sqrt(smallM) * (test_ccuts_bootstrap3[(1:BigK)+(m-1)*BigK, i+(t1-1)*length(nvec3)] - rep(xcut_exact(generate_graph(t_vec[t1], matrix(current_sample, nrow=M)))$ccut_min, BigK))
      setTxtProgressBar(progress, m+(i-1)*BigM+(t1-1)*BigM*length(nvec3))
    }
  }
}
close(progress)

save(test_ccuts_bootstrap4, file="ccut_bootstrap4_stat_BigMK100_meqn45_n10h4-13.Rdata")
#load("ccut_bootstrap4_stat_BigMK100_meqsqrtn_n10h4-6.Rdata")
nc_ks_cbimodal <- sapply(1:length(nvec3), function(i)sapply(1:length(t_vec), function(t1)sum(sapply(1:BigM,
                          function(m)ks.test(test_ccuts_bootstrap4[(1:BigK)+(m-1)*BigK, i+(t1-1)*length(nvec3)], ncut_true_sample_list_over_t[[t1]])$p.value) > 0.05)))

hist(test_ccuts_bootstrap4, freq=FALSE)
hist(ncut_true_sample_list_over_t[[3]], freq=FALSE)

nc_ks_cbimodal <- sapply(1:length(nvec3), function(i)sapply(1:length(t_vec), function(t1)
  ks.test(test_ccuts_bootstrap4[,i+(t1-1)*length(nvec3)], ncut_true_sample_list_over_t[[t1]])$statistic), simplify="array")

rmultinom2 <- function(n=1, size, prob) rowSums(rmultinom(n=1e4, size/1e4, prob=prob))

####---- Third testing range: FCPS

library(RANN)
library(threejs)

generate_sample_graph_general <- function(raw_sample, useGaussianWeights=FALSE, t=NA, k=5, sigm=.8){
  
  d <- dim(raw_sample)[2]
  n <- dim(raw_sample)[1]
  # Determines sample size n and number of dimensions d.
  
  adj_mat <- matrix(rep(0,n^2), n, n)
  edge_mat <- c()
  # Initializing. adj_mat will be the weight (or adjacency) matrix of the generated graph, edge list will be a vector
  # containing all edges and capacity_list will contain the same weights for using min_cuts
  
  #dist_mat <- outer(1:n, 1:n, function(i,j)euclidean_distance(raw_sample[i,], raw_sample[j,]))
  
  knn_sample <- RANN::nn2(raw_sample, k=k+1)
  # Computes the k nearest neighbors of each point and their respective distances.
  # The nearest point is always the point itself which will be ignored (therefore k+1, not k).
  
  if(is.na(t)){
    t <- median(knn_sample$nn.dists[,2])
    print(paste("Choosing t =",t))
  }
  
  for(i in 1:n){
  
    #if(dij <= t || dij == min(dist_mat[i,]) || dij == min(dist_mat[,j])){
    for(ki in 2:(k+1)){
      # Ignore the first entry of knn_sample since this is always the point i itself
      
      ki_idx <- knn_sample$nn.idx[i,ki]
      
      edge_mat <- rbind(edge_mat, c(min(i, ki_idx), max(i, ki_idx)))
      # Adds an edge between i and ki_idx. min/max to avoid counting edges double.
      
      adj_mat[i, ki_idx] <- ifelse(useGaussianWeights, exp(-knn_sample$nn.dists[i,ki] / (2 * sigm^2)), 1)
      # Should the graph be built using gaussian weights? If so, sigm is the parameter.
      # If not, the weights are given by w_{ij} = 1.
      
    }
    
    for(j in (i+1):n){

      if(i<n && euclidean_distance(raw_sample[i,], raw_sample[j,]) <= t){
        # Are nodes i and j connected?
        
        edge_mat <- rbind(edge_mat, c(min(i, j), max(i, j)))
        # Adds an edge between i and j. min/max are just for safety (since always i<j).
        
        adj_mat[i,j] <- ifelse(useGaussianWeights, exp(-euclidean_distance(raw_sample[i,], raw_sample[j,]) / (2 * sigm^2)), 1)
        # Should the graph be built using gaussian weights? If so, sigm is the parameter.
        # If not, the weights are given by w_{ij} = 1.
        
      }
    }
  }
  
  adj_mat <- matrix(mapply(max, adj_mat, t(adj_mat)), n)
  #adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
  # Since it is known that adj_mat is symmetric, this line copies all weights from the upper triangular part of the matrix
  # to the bottom one, ensuring symmetry.
  
  edge_mat <- unique(edge_mat)
  g <- graph(edges=c(t(edge_mat)), n=n, directed = FALSE)
  # Creates a graph with n vertices and only the unique edges from edge_mat.
  
  E(g)$capacity <- mapply(function(i,j)adj_mat[i,j], edge_mat[,1], edge_mat[,2])
  # Assigns capacities (=weights) of the edges from edge_list. This is needed for igraph::min_cut.
  
  E(g)$weight <- mapply(function(i,j)adj_mat[i,j], edge_mat[,1], edge_mat[,2])#rep(100/n, dim(edge_mat)[1])
  E(g)$width <- rep(10, dim(edge_mat)[1])
  # This is used solely for plotting purposes. To balance edge thickness, their capacities are normalized and scaled.
  
  E(g)$distance <- sapply(1:dim(edge_mat)[1], function(i)euclidean_distance(raw_sample[edge_mat[i,1],], raw_sample[edge_mat[i,2],]))
  # Assigns euclidean distances between vertices constituting edges
  
  V(g)$x <- raw_sample[,1]
  if(d > 1) V(g)$y <- raw_sample[,2]
  if(d > 2) V(g)$z <- raw_sample[,3]
  # Assigns the first three components (if applicable) as coordinates for plotting
  
  V(g)$label <- sapply(1:n, toString)
  V(g)$sample <- rep(1, n)
  # Assigns labels to all vertices.
  
  degree_vec <- rowSums(adj_mat)
  V(g)$degree <- degree_vec
  # Computes the degree (i.e. the sum of the weights of all neighbouring vertices) for each vertex.
  
  V(g)$local_maximum <- sapply(1:n, function(i)max((adj_mat[i,] > 0) * degree_vec) <= degree_vec[i])
  # Determines whether each vertex i is a local maximum, i.e. if deg(i) = max{deg(j) | j~i}.
  
  return(list(graph=g, weight_mat=adj_mat, removed_vertices=NULL))
  # Returns a list consisting of the graph g, associated weight matrix weight_mat
  # and a list of vertices removed due to not containing any observations
}
# Generates an extended graph object containing a graph constructed using t-radius neighbourhood,
# associated weight matrix with weights w_ij = 1 for observations X_i and X_j if i~j.
# Input: Distance t and a matrix with two columns containing the two-dim. sample.
# Output: list(graph, weight_mat, removed_vertices=NULL)


#library(FCPS)
#timevec <- c()
#maxvec <- c()
#for(nsmpls in seq(100, 1000, 100)){
engytime_sample <- FCPS::EngyTime$Data[sample(dim(FCPS::EngyTime$Data)[1], 2000, FALSE), ]
engytime <- generate_sample_graph_general(engytime_sample, useGaussianWeights=TRUE, t=.25, k=5, sigm=0.8)
#engytime <- generate_sample_graph_general(1, FCPS::EngyTime$Data, useGaussianWeights=TRUE, k=5, sigm=0.5)
st <- Sys.time()
engytime_cut <- xcut_via_stmincut_on_local_maxima(engytime)
et <- Sys.time()
et-st
#timevec <- c(timevec, et-st)
#maxvec <- c(maxvec, sum(V(engytime$graph)$local_maximum))
#print(nsmpls)
#}
atom_dataset <- FCPS::GolfBall
atom_sample <- atom_dataset$Data
atom_graph <- generate_sample_graph_general(atom_sample, k=15, t=0.05, sigm=.2, useGaussianWeights = T) #15,0.05,0.4 #15,0.05,F
atom_cut <- xvstlocmax_merge(atom_graph)
plot_graph_partition(atom_graph, atom_cut$ncut, vertex_scale=1, vertex_offset=1, edge_scale=.5, edge_offset=.5)

atom_ncut_vec <- rep(1, dim(atom_sample)[1])
atom_ncut_vec[atom_cut$ncut] <- 2
atom_cut2 <- manual_spec_clustering(atom_graph, method="kmeans")$norm
#atom_cut3 <- cluster_leiden(atom_graph$graph, resolution_parameter=1e-7, n_iterations=100)$membership
#table(atom_cut3)
resseq <- (1e-10)*1.1^(1:250)
atom_cut3 <- which.min(sapply(resseq, function(res)reorder_dbscan(cluster_leiden(atom_graph$graph, resolution_parameter=res, n_iterations=100)$membership - 1) %>%
                                (function(a)min(sum(atom_dataset$Cls!=a), sum(atom_dataset$Cls!=3-a))))) %>%
                                (function(i)cluster_leiden(atom_graph$graph, resolution_parameter=resseq[i], n_iterations=100)$membership)
#atom_cut4 <- dbscan::dbscan(atom_sample, eps=3)$cluster
#table(atom_cut4)
epsseq <- seq(0,5,0.01)
atom_cut4 <- which.min(sapply(epsseq, function(eps3)reorder_dbscan(dbscan::dbscan(atom_sample, eps=eps3)$cluster) %>%
                                (function(a)min(sum(atom_dataset$Cls!=a), sum(atom_dataset$Cls!=3-a))))) %>%
                                (function(i)reorder_dbscan(dbscan::dbscan(atom_sample, eps=epsseq[i])$cluster))

print(paste(c("Xvst:","SpecClust:","Leiden:","DBSCAN:"),
            1 - c(min(sum(atom_dataset$Cls!=atom_ncut_vec), sum(atom_dataset$Cls!=3-atom_ncut_vec)),
                  min(sum(atom_dataset$Cls!=atom_cut2), sum(atom_dataset$Cls!=3-atom_cut2)),
                  min(sum(atom_dataset$Cls!=atom_cut3), sum(atom_dataset$Cls!=3-atom_cut3)),
                  min(sum(atom_dataset$Cls!=atom_cut4), sum(atom_dataset$Cls!=3-atom_cut4)))/dim(atom_sample)[1]))

reorder_dbscan <- function(a){
  b1 <- as.integer(which.max(table(a)))
  b2 <- as.integer(which.max(table(a)[-b1]))
  a[a==b1-1] <- -2
  a[a==b2-1] <- -1
  return(as.integer(factor(a)))
}


V(atom_graph$graph)$color <- rep("green", length(V(atom_graph$graph)))
V(atom_graph$graph)$color[atom_cut$ncut] <- "blue"
#graphjs(atom_graph$graph,vertex.size = 1, showLabels=TRUE)
if(dim(atom_sample)[2]==3) colnames(atom_sample) <- c("X1", "X2", "X3")
threejs::graphjs(atom_graph$graph, layout=atom_sample, showLabels = TRUE, vertex.size=.5, edge.width=3)#%>%points3d(vertices(.), color="red",pch=V(g)$label,size=0.1)
scatterplot3d::scatterplot3d(atom_sample, color=V(atom_graph$graph)$color, box=FALSE,
                xlim=c(-1,1), ylim=c(-1,2), zlim=c(-1,1), xlab="x", ylab="y", zlab="z",
                pch=20, scale.y=0.75, asp=.5, angle=260)

n <- 1000
atom_sample <- rbind(mvtnorm::rmvnorm(n, c(2,0.2), sigma=0.01*diag(2)), cbind(runif(n,0,8), runif(n,-0.05,0)))
atom_leiden <- cluster_leiden(atom_graph$graph, resolution_parameter=1e-4, n_iterations=100)
atom_leiden <- cluster_louvain(atom_graph$graph, resolution=1e-2)
table(atom_leiden$membership)
plot_graph_partition(atom_graph, atom_leiden$membership, vertex_scale=.1, vertex_offset=2, edge_scale=2, edge_offset=.5)
atom_spec_clust <- manual_spec_clustering(atom_graph, method="kmeans")
plot_graph_partition(atom_graph, atom_spec_clust$norm, vertex_scale=.1, vertex_offset=2, edge_scale=2, edge_offset=.5)

atom_dbscan <- dbscan::dbscan(atom_sample, eps=median(RANN::nn2(atom_sample, k=2)$nn.dists[,2]))$cluster
atom_dbscan <- dbscan::dbscan(atom_sample, eps=.04)$cluster
table(atom_dbscan)
plot_graph_partition(atom_graph, atom_dbscan+1, vertex_scale=.1, vertex_offset=2, edge_scale=2, edge_offset=.5)

library(scatterplot3d)
sactter.grid <- function (x, y = NULL, z = NULL, color = par("col"), pch = NULL, 
                          main = NULL, sub = NULL, xlim = NULL, ylim = NULL, zlim = NULL, 
                          xlab = NULL, ylab = NULL, zlab = NULL, scale.y = 1, angle = 40, 
                          axis = TRUE, tick.marks = TRUE, label.tick.marks = TRUE, 
                          x.ticklabs = NULL, y.ticklabs = NULL, z.ticklabs = NULL, 
                          y.margin.add = 0, grid = TRUE, box = TRUE, lab = par("lab"), 
                          lab.z = mean(lab[1:2]), type = "p", highlight.3d = FALSE, 
                          mar = c(5, 3, 4, 3) + 0.1, bg = par("bg"), col.axis = par("col.axis"), 
                          col.grid = "grey", col.lab = par("col.lab"), cex.symbols = par("cex"), 
                          cex.axis = 0.8 * par("cex.axis"), cex.lab = par("cex.lab"), 
                          font.axis = par("font.axis"), font.lab = par("font.lab"), 
                          lty.axis = par("lty"), lty.grid = par("lty"), lty.hide = NULL, 
                          lty.hplot = par("lty"), log = "", ...) 
{
  mem.par <- par(mar = mar)
  x.scal <- y.scal <- z.scal <- 1
  xlabel <- if (!missing(x)) 
    deparse(substitute(x))
  ylabel <- if (!missing(y)) 
    deparse(substitute(y))
  zlabel <- if (!missing(z)) 
    deparse(substitute(z))
  if (highlight.3d && !missing(color)) 
    warning("color is ignored when highlight.3d = TRUE")
  if (!is.null(d <- dim(x)) && (length(d) == 2) && (d[2] >= 
                                                    4)) 
    color <- x[, 4]
  else if (is.list(x) && !is.null(x$color)) 
    color <- x$color
  xyz <- xyz.coords(x = x, y = y, z = z, xlab = xlabel, ylab = ylabel, 
                    zlab = zlabel, log = log)
  if (is.null(xlab)) {
    xlab <- xyz$xlab
    if (is.null(xlab)) 
      xlab <- ""
  }
  if (is.null(ylab)) {
    ylab <- xyz$ylab
    if (is.null(ylab)) 
      ylab <- ""
  }
  if (is.null(zlab)) {
    zlab <- xyz$zlab
    if (is.null(zlab)) 
      zlab <- ""
  }
  if (length(color) == 1) 
    color <- rep(color, length(xyz$x))
  else if (length(color) != length(xyz$x)) 
    stop("length(color) ", "must be equal length(x) or 1")
  angle <- (angle%%360)/90
  yz.f <- scale.y * abs(if (angle < 1) angle else if (angle > 
                                                      3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if (angle < 2) 
    1 - angle
    else angle - 3)
  if (angle > 2) {
    temp <- xyz$x
    xyz$x <- xyz$y
    xyz$y <- temp
    temp <- xlab
    xlab <- ylab
    ylab <- temp
    temp <- xlim
    xlim <- ylim
    ylim <- temp
  }
  angle.1 <- (1 < angle && angle < 2) || angle > 3
  angle.2 <- 1 <= angle && angle <= 3
  dat <- cbind(as.data.frame(xyz[c("x", "y", "z")]), col = color)
  if (!is.null(xlim)) {
    xlim <- range(xlim)
    dat <- dat[xlim[1] <= dat$x & dat$x <= xlim[2], , drop = FALSE]
  }
  if (!is.null(ylim)) {
    ylim <- range(ylim)
    dat <- dat[ylim[1] <= dat$y & dat$y <= ylim[2], , drop = FALSE]
  }
  if (!is.null(zlim)) {
    zlim <- range(zlim)
    dat <- dat[zlim[1] <= dat$z & dat$z <= zlim[2], , drop = FALSE]
  }
  n <- nrow(dat)
  if (n < 1) 
    stop("no data left within (x|y|z)lim")
  y.range <- range(dat$y[is.finite(dat$y)])
  if (type == "p" || type == "h") {
    y.ord <- rev(order(dat$y))
    dat <- dat[y.ord, ]
    if (length(pch) > 1) 
      if (length(pch) != length(y.ord)) 
        stop("length(pch) ", "must be equal length(x) or 1")
    else pch <- pch[y.ord]
    if (length(bg) > 1) 
      if (length(bg) != length(y.ord)) 
        stop("length(bg) ", "must be equal length(x) or 1")
    else bg <- bg[y.ord]
    if (length(cex.symbols) > 1) 
      if (length(cex.symbols) != length(y.ord)) 
        stop("length(cex.symbols) ", "must be equal length(x) or 1")
    else cex.symbols <- cex.symbols[y.ord]
    daty <- dat$y
    daty[!is.finite(daty)] <- mean(daty[is.finite(daty)])
    if (highlight.3d && !(all(diff(daty) == 0))) 
      dat$col <- rgb(red = seq(0, 1, length = n) * (y.range[2] - 
                                                      daty)/diff(y.range), green = 0, blue = 0)
  }
  p.lab <- par("lab")
  y.range <- range(dat$y[is.finite(dat$y)], ylim)
  y.prty <- pretty(y.range, n = lab[2], min.n = max(1, min(0.5 * 
                                                             lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  dat$y <- (dat$y - y.add)/y.scal
  y.max <- (max(y.prty) - y.add)/y.scal
  if (!is.null(ylim)) 
    y.max <- max(y.max, ceiling((ylim[2] - y.add)/y.scal))
  x.range <- range(dat$x[is.finite(dat$x)], xlim)
  x.prty <- pretty(x.range, n = lab[1], min.n = max(1, min(0.5 * 
                                                             lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  dat$x <- dat$x/x.scal
  x.range <- range(x.prty)/x.scal
  x.max <- ceiling(x.range[2])
  x.min <- floor(x.range[1])
  if (!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2]/x.scal))
    x.min <- min(x.min, floor(xlim[1]/x.scal))
  }
  x.range <- range(x.min, x.max)
  z.range <- range(dat$z[is.finite(dat$z)], zlim)
  z.prty <- pretty(z.range, n = lab.z, min.n = max(1, min(0.5 * 
                                                            lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  dat$z <- dat$z/z.scal
  z.range <- range(z.prty)/z.scal
  z.max <- ceiling(z.range[2])
  z.min <- floor(z.range[1])
  if (!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2]/z.scal))
    z.min <- min(z.min, floor(zlim[1]/z.scal))
  }
  z.range <- range(z.min, z.max)
  plot.new()
  if (angle.2) {
    x1 <- x.min + yx.f * y.max
    x2 <- x.max
  }
  else {
    x1 <- x.min
    x2 <- x.max + yx.f * y.max
  }
  plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
  temp <- strwidth(format(rev(y.prty))[1], cex = cex.axis/par("cex"))
  if (angle.2) 
    x1 <- x1 - temp - y.margin.add
  else x2 <- x2 + temp + y.margin.add
  plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
  if (angle > 2) 
    par(usr = par("usr")[c(2, 1, 3:4)])
  usr <- par("usr")
  title(main, sub, ...)
  if ("xy" %in% grid || grid) {
    i <- x.min:x.max
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + 
               z.min, col = col.grid, lty = lty.grid)
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min, x.max + 
               (i * yx.f), i * yz.f + z.min, col = col.grid, lty = lty.grid)
  }
  if ("xz" %in% grid) {
    i <- x.min:x.max
    segments(i + (yx.f * y.max), yz.f * y.max + z.min, 
             i + (yx.f * y.max), yz.f * y.max + z.max, 
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.max + temp,temp1 + i , col = col.grid, lty = lty.grid)
    
  }
  
  if ("yz" %in% grid) {
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min,  
             x.min + (i * yx.f) ,i * yz.f + z.max,  
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.min, i , col = col.grid, lty = lty.grid)
    
    
    
  }
  
  if (axis) {
    xx <- if (angle.2) 
      c(x.min, x.max)
    else c(x.max, x.min)
    if (tick.marks) {
      xtl <- (z.max - z.min) * (tcl <- -par("tcl"))/50
      ztl <- (x.max - x.min) * tcl/50
      mysegs <- function(x0, y0, x1, y1) segments(x0, 
                                                  y0, x1, y1, col = col.axis, lty = lty.axis)
      i.y <- 0:y.max
      mysegs(yx.f * i.y - ztl + xx[1], yz.f * i.y + z.min, 
             yx.f * i.y + ztl + xx[1], yz.f * i.y + z.min)
      i.x <- x.min:x.max
      mysegs(i.x, -xtl + z.min, i.x, xtl + z.min)
      i.z <- z.min:z.max
      mysegs(-ztl + xx[2], i.z, ztl + xx[2], i.z)
      if (label.tick.marks) {
        las <- par("las")
        mytext <- function(labels, side, at, ...) mtext(text = labels, 
                                                        side = side, at = at, line = -0.5, col = col.lab, 
                                                        cex = cex.axis, font = font.lab, ...)
        if (is.null(x.ticklabs)) 
          x.ticklabs <- format(i.x * x.scal)
        mytext(x.ticklabs, side = 1, at = i.x)
        if (is.null(z.ticklabs)) 
          z.ticklabs <- format(i.z * z.scal)
        mytext(z.ticklabs, side = if (angle.1) 
          4
          else 2, at = i.z, adj = if (0 < las && las < 
                                      3) 
            1
          else NA)
        temp <- if (angle > 2) 
          rev(i.y)
        else i.y
        if (is.null(y.ticklabs)) 
          y.ticklabs <- format(y.prty)
        else if (angle > 2) 
          y.ticklabs <- rev(y.ticklabs)
        text(i.y * yx.f + xx[1], i.y * yz.f + z.min, 
             y.ticklabs, pos = if (angle.1) 
               2
             else 4, offset = 1, col = col.lab, cex = cex.axis/par("cex"), 
             font = font.lab)
      }
    }
    mytext2 <- function(lab, side, line, at) mtext(lab, 
                                                   side = side, line = line, at = at, col = col.lab, 
                                                   cex = cex.lab, font = font.axis, las = 0)
    lines(c(x.min, x.max), c(z.min, z.min), col = col.axis, 
          lty = lty.axis)
    mytext2(xlab, 1, line = 1.5, at = mean(x.range))
    lines(xx[1] + c(0, y.max * yx.f), c(z.min, y.max * yz.f + 
                                          z.min), col = col.axis, lty = lty.axis)
    mytext2(ylab, if (angle.1) 
      2
      else 4, line = 0.5, at = z.min + y.max * yz.f)
    lines(xx[c(2, 2)], c(z.min, z.max), col = col.axis, 
          lty = lty.axis)
    mytext2(zlab, if (angle.1) 
      4
      else 2, line = 1.5, at = mean(z.range))
    if (box) {
      if (is.null(lty.hide)) 
        lty.hide <- lty.axis
      temp <- yx.f * y.max
      temp1 <- yz.f * y.max
      lines(c(x.min + temp, x.max + temp), c(z.min + temp1, 
                                             z.min + temp1), col = col.axis, lty = lty.hide)
      lines(c(x.min + temp, x.max + temp), c(temp1 + z.max, 
                                             temp1 + z.max), col = col.axis, lty = lty.axis)
      temp <- c(0, y.max * yx.f)
      temp1 <- c(0, y.max * yz.f)
      lines(temp + xx[2], temp1 + z.min, col = col.axis, 
            lty = lty.hide)
      lines(temp + x.min, temp1 + z.max, col = col.axis, 
            lty = lty.axis)
      temp <- yx.f * y.max
      temp1 <- yz.f * y.max
      lines(c(temp + x.min, temp + x.min), c(z.min + temp1, 
                                             z.max + temp1), col = col.axis, lty = if (!angle.2) 
                                               lty.hide
            else lty.axis)
      lines(c(x.max + temp, x.max + temp), c(z.min + temp1, 
                                             z.max + temp1), col = col.axis, lty = if (angle.2) 
                                               lty.hide
            else lty.axis)
    }
  }
  x <- dat$x + (dat$y * yx.f)
  z <- dat$z + (dat$y * yz.f)
  col <- as.character(dat$col)
  if (type == "h") {
    z2 <- dat$y * yz.f + z.min
    segments(x, z, x, z2, col = col, cex = cex.symbols, 
             lty = lty.hplot, ...)
    points(x, z, type = "p", col = col, pch = pch, bg = bg, 
           cex = cex.symbols, ...)
  }
  else points(x, z, type = type, col = col, pch = pch, bg = bg, 
              cex = cex.symbols, ...)
  if (axis && box) {
    lines(c(x.min, x.max), c(z.max, z.max), col = col.axis, 
          lty = lty.axis)
    lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + 
            z.max, col = col.axis, lty = lty.axis)
    lines(xx[c(1, 1)], c(z.min, z.max), col = col.axis, 
          lty = lty.axis)
  }
  ob <- ls()
  rm(list = ob[!ob %in% c("angle", "mar", "usr", "x.scal", 
                          "y.scal", "z.scal", "yx.f", "yz.f", "y.add", "z.min", 
                          "z.max", "x.min", "x.max", "y.max", "x.prty", "y.prty", 
                          "z.prty")])
  rm(ob)
  invisible(list(xyz.convert = function(x, y = NULL, z = NULL) {
    xyz <- xyz.coords(x, y, z)
    if (angle > 2) {
      temp <- xyz$x
      xyz$x <- xyz$y
      xyz$y <- temp
    }
    y <- (xyz$y - y.add)/y.scal
    return(list(x = xyz$x/x.scal + yx.f * y, y = xyz$z/z.scal + 
                  yz.f * y))
  }, points3d = function(x, y = NULL, z = NULL, type = "p", 
                         ...) {
    xyz <- xyz.coords(x, y, z)
    if (angle > 2) {
      temp <- xyz$x
      xyz$x <- xyz$y
      xyz$y <- temp
    }
    y2 <- (xyz$y - y.add)/y.scal
    x <- xyz$x/x.scal + yx.f * y2
    y <- xyz$z/z.scal + yz.f * y2
    mem.par <- par(mar = mar, usr = usr)
    on.exit(par(mem.par))
    if (type == "h") {
      y2 <- z.min + yz.f * y2
      segments(x, y, x, y2, ...)
      points(x, y, type = "p", ...)
    } else points(x, y, type = type, ...)
  }, plane3d = function(Intercept, x.coef = NULL, y.coef = NULL, 
                        lty = "dashed", lty.box = NULL, ...) {
    if (!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
    if (is.null(lty.box)) lty.box <- lty
    if (is.null(x.coef) && length(Intercept) == 3) {
      x.coef <- Intercept[if (angle > 2) 3 else 2]
      y.coef <- Intercept[if (angle > 2) 2 else 3]
      Intercept <- Intercept[1]
    }
    mem.par <- par(mar = mar, usr = usr)
    on.exit(par(mem.par))
    x <- x.min:x.max
    ltya <- c(lty.box, rep(lty, length(x) - 2), lty.box)
    x.coef <- x.coef * x.scal
    z1 <- (Intercept + x * x.coef + y.add * y.coef)/z.scal
    z2 <- (Intercept + x * x.coef + (y.max * y.scal + y.add) * 
             y.coef)/z.scal
    segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, 
             lty = ltya, ...)
    y <- 0:y.max
    ltya <- c(lty.box, rep(lty, length(y) - 2), lty.box)
    y.coef <- (y * y.scal + y.add) * y.coef
    z1 <- (Intercept + x.min * x.coef + y.coef)/z.scal
    z2 <- (Intercept + x.max * x.coef + y.coef)/z.scal
    segments(x.min + y * yx.f, z1 + y * yz.f, x.max + y * 
               yx.f, z2 + y * yz.f, lty = ltya, ...)
  }, box3d = function(...) {
    mem.par <- par(mar = mar, usr = usr)
    on.exit(par(mem.par))
    lines(c(x.min, x.max), c(z.max, z.max), ...)
    lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + 
            z.max, ...)
    lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + 
            z.max, ...)
    lines(c(x.max, x.max), c(z.min, z.max), ...)
    lines(c(x.min, x.min), c(z.min, z.max), ...)
    lines(c(x.min, x.max), c(z.min, z.min), ...)
  }))
}
sactter.grid(atom_sample, color=V(atom_cut_iterated$g$graph)$color, box=FALSE, y.margin.add=.05,
             xlim=c(-4,4), ylim=c(-4,4), zlim=c(-4,4), xlab="", ylab="", zlab="",
             pch=20, scale.y=0.75, asp=.5, angle=45, grid=c('xy','xz','yz'), cex.symbols=1.5)

xvst_locmax_iterated <- function(g, min_size=5, cut="ncut"){
  
  while(length(partition) >= 2*min_size && T){
    # Only accept partitions that are of min_size - we always divide into two.
    
    cuts <- xcut_via_stmincut_on_local_maxima(g)
    
    if(cut=="ccut"){
      cut_value <- cuts$ccut_value
      cut <- cuts$ccut
    }
    else if(cut=="rcut"){
      cut_value <- cuts$rcut_value
      cut <- cuts$rcut
    }
    else if(cut=="mcut"){
      cut_value <- cuts$mcut_value
      cut <- cuts$mcut
    }
    else{
      if(cut!="ncut") print("Invalid cut value. Defaulting to ncut.")
      cut_value <- cuts$ncut_value
      cut <- cuts$ncut
    }
    # Selecting the correct cut and corresponding cut value.
    
    all_divisions[[index_best_cut]] <- cut
    all_divisions <- append(all_divisions, list(cut))
    all_cut_values <- c(all_cut_values, cut_value)
    # 
    
    index_best_cut <- which.max(all_cut_values)
    
  }
  
  
}

atom_cut_iterated <- xvst_iterated(4, atom_graph)
threejs::graphjs(atom_cut_iterated$g$graph, layout=atom_sample, showLabels = TRUE, vertex.size=.5, edge.width=3)#%>%points3d(vertices(.), color="red",pch=V(g)$label,size=0.1)

n <- 100
BigN <- 100
distrange <- seq(2,6,0.1)
#distrange <- seq(0.1,1,0.05)
#distrange <- seq(0.1,3,0.1)
#distrange <- c(1)
resseq <- (1e-8)*1.1^(1:150)
epsseq <- seq(0.02,5,0.02)
vmmat <- vmmat2 <- vmmat3 <- vmmat4 <- matrix(nrow=length(distrange), ncol=3*BigN)
#cmmat <- cmmat2 <- cmmat3 <- cmmat4 <- matrix(nrow=length(distrange), ncol=3*BigN)
progress <- txtProgressBar(min = 0, max = length(distrange)*BigN, style = 3)
#st <- Sys.time()
for(dst in 1:length(distrange)){
  for(i in 1:BigN){
    n1 <- rbinom(1, n, 1/2)
#    gauss_cut1 <- list(ncut=NULL,ncut_min=Inf)
#  while(is.null(gauss_cut1) || is.null(gauss_cut1$ncut)){
    #gauss_mix_sample <- mvtnorm::rmvnorm(n, c(0,0), sigma=distrange[dst]*diag(2))
    gauss_mix_sample <- rbind(mvtnorm::rmvnorm(n1, c(0,0), sigma=diag(2)), mvtnorm::rmvnorm(n-n1, c(0,distrange[dst]), sigma=diag(2)))
    #gauss_mix_sample <- rbind(mvtnorm::rmvnorm(n1, c(0,0), sigma=matrix(c(2,-1,-1,2),nrow=2)), mvtnorm::rmvnorm(n-n1, distrange[dst]*c(1,1), sigma=matrix(c(2,1,1,2),nrow=2)))
    #gauss_mix_sample <- rbind(mvtnorm::rmvnorm(n1, c(2,0.2), sigma=0.01*diag(2)), cbind(runif(n-n1,0,8), runif(n-n1,-0.05,0)))
    gauss_mix_clusters <- c(rep(1,n1),rep(2,n-n1))
    #gauss_mix_clusters <- rep(1,n)
    gauss_graph <- generate_sample_graph(gauss_mix_sample, t=.2, sigm=.2, nearest_neighbors=5)
    #gauss_graph <- generate_sbm_graph(c(n1,n-n1), distrange[dst], 0.1)
    
    gauss_cut1 <- xist_gusfield(gauss_graph, locmax_by_degree=T, progress_bar=F)
#  }
    
    gauss_cut <- gauss_cut1$ncut
    gauss_cut_trafo <- rep(1,n) %>% (function(a){a[gauss_cut] <- 2;a})
    #gauss_cut_trafo[gauss_cut] <- 2
    gauss_cut_spec <- manual_multiway_spec_clustering(gauss_graph, k=2)[[1]] #manual_spec_clustering(gauss_graph, method="kmeans")$norm

    vmmat[dst, i] <- FDRSeg::v_measure(gauss_mix_clusters, gauss_cut_trafo)
    vmmat2[dst, i] <- FDRSeg::v_measure(gauss_mix_clusters, gauss_cut_spec)
    vmmat3[dst, i] <- max(sapply(resseq, function(res)FDRSeg::v_measure(gauss_mix_clusters,
                            reorder_dbscan(cluster_leiden(gauss_graph, resolution_parameter=res)$membership - 1))))
    vmmat4[dst, i] <- max(sapply(epsseq, function(eps3)FDRSeg::v_measure(gauss_mix_clusters,
                            reorder_dbscan(dbscan::dbscan(gauss_mix_sample, eps=eps3)$cluster))))
    
    vmmat[dst, i+BigN] <- max(mean(gauss_mix_clusters==gauss_cut_trafo), mean(gauss_mix_clusters==3-gauss_cut_trafo))
    vmmat2[dst, i+BigN] <- max(mean(gauss_mix_clusters==gauss_cut_spec), mean(gauss_mix_clusters==3-gauss_cut_spec))
    vmmat3[dst, i+BigN] <- max(sapply(resseq, function(res)reorder_dbscan(cluster_leiden(gauss_graph, resolution_parameter=res)$membership - 1) %>%
                                    (function(a)max(mean(gauss_mix_clusters==a), mean(gauss_mix_clusters==3-a)))))
    vmmat4[dst, i+BigN] <- max(sapply(epsseq, function(eps3)reorder_dbscan(dbscan::dbscan(gauss_mix_sample, eps=eps3)$cluster) %>%
                                    (function(a)max(mean(gauss_mix_clusters==a), mean(gauss_mix_clusters==3-a)))))
    
    #print(i)
    
    vmmat[dst, i+2*BigN] <- xcut_partition_S(gauss_cut, gauss_graph)$ncut_min
    vmmat2[dst, i+2*BigN] <- xcut_partition_S(which(gauss_cut_spec==1), gauss_graph)$ncut_min
    vmmat3[dst, i+2*BigN] <- min(sapply(resseq, function(res)xcut_multipartition_S(gauss_graph, cluster_leiden(gauss_graph, resolution_parameter=res)$membership)$ncut_min))
    vmmat4[dst, i+2*BigN] <- min(sapply(epsseq, function(eps3)xcut_multipartition_S(gauss_graph, dbscan::dbscan(gauss_mix_sample, eps=eps3)$cluster)$ncut_min))
    
    #print(Sys.time()-st)
    setTxtProgressBar(progress, i+(dst-1)*BigN)
  }
}
close(progress)
#load("4algs_gaussians_k5_t0p2_sigm0p2_n1000_BigN100.RData")
save(vmmat, vmmat2, vmmat3, vmmat4, file="4algs_gaussians_k5_t0p2_sigm0p2_n1000_BigN50_1_part1.RData")

plot_graph_partition(gauss_graph, gauss_cut, vertex_scale=1, vertex_offset=1, edge_scale=.5, edge_offset=.5)

plot(distrange, sapply(1:length(distrange), function(i)mean(vmmat2[i,])), type="l")

algs_gauss <- tidyr::pivot_longer(data.frame(cbind(sapply(1:length(distrange), function(i)median(vmmat[i,2*BigN+1:BigN],na.rm=T)),
                                                  sapply(1:length(distrange), function(i)median(vmmat2[i,2*BigN+1:BigN],na.rm=T)),
                                                  sapply(1:length(distrange), function(i)median(vmmat3[i,2*BigN+1:BigN],na.rm=T)),
                                                  sapply(1:length(distrange), function(i)median(vmmat4[i,2*BigN+1:BigN],na.rm=T))), row.names=NULL),
                               cols=1:4, names_to="algs", values_to="values")
algs_gauss_new <- cbind(algs_gauss, x=rep(distrange, each=4))
library(scales)
library(latex2exp)
gg3 <- ggplot(algs_gauss_new, aes(x=x, y=values, colour=algs)) +
  geom_line(aes(x=x, y=values, color=algs), size=2.5) +
  #geom_vline(aes(xintercept=x), linetype=2, color="grey") +
#  scale_y_continuous(trans = scales::trans_new(name="almost_log", function(x)log10(x+1e-5), function(x)10^x-1e-5),
#                     breaks=c(1e-6,1e-5,1e-4,1e-3), labels=c(expression(10^{-6}),expression(10^{-5}),expression(10^{-4}),expression(10^{-3}))) + #unname(TeX(c("$10^{-6}$","$10^{-5}$","$10^{-4}$","$10^{-3}$")))) + #ylim(c(0,0.0015)) +
  #geom_hline(yintercept = 0.042947, colour="blue", linetype="dashed", show.legend=FALSE) +
  #geom_point() + geom_point(aes(x=xend, y=values), shape=1) + geom_segment(key_glyph="path") + coord_cartesian(xlim=c(3,6)) +
  #labs(x=unname(TeX("Distance \\delta")), y="Classification rate") + 
  labs(x=expression("Distance"~delta), y="NCut value") + 
  theme(text = element_text(size=rel(4)), axis.text=element_text(size=14), panel.background = element_rect(fill="gray98", colour="white")) +
        #panel.grid.major = element_line(color="grey"), panel.grid.minor = element_line(color="lightgrey")) + #theme_minimal() +
  scale_color_manual(name="Algorithm", values = c("orangered1", "orange1", "yellow2", "mediumpurple2"),
                     labels=c("Xist", "Spectral clustering", "Leiden (oracle)", "DBSCAN (oracle)"),
                     guide = guide_legend(override.aes = list(linetype = rep("solid", 4), shape = rep(NA, 4)))) + #theme(legend.position="none")
  #theme(legend.position = "bottom", legend.box.background = element_rect(colour="black"))
  #scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),
  #                   labels = trans_format("log2", math_format(2^.x))) +
  #ylim(c(0,7e-6))+
  theme(legend.position = c(0.75,0.2), legend.box.background = element_rect(colour = "black"), legend.text=element_text(size=14))
gg3

st <- Sys.time()
a <- xvstlocmax_merge(atom_graph)$ncut
et <- Sys.time()
et-st

generate_sbm_graph <- function(Mvec, p, q){
  
  partitionvec <- unlist(sapply(1:length(Mvec),function(i)rep(i,Mvec[i])))
  M <- length(partitionvec)
  
  adj_mat <- matrix(rep(0,M^2),nrow=M)
  edge_list <- c()
  #degree_vec <- rep(0,M)
  for(i in 1:(M-1)){
    for(j in (i+1):M){
      if(partitionvec[i]==partitionvec[j]) decider <- rbinom(1,1,p)
      else decider <- rbinom(1,1,q)
      if(decider){
        edge_list <- c(edge_list, i,j)
        adj_mat[i,j] <- adj_mat[j,i] <- 1
        #degree_vec[c(i,j)] <- degree_vec[c(i,j)] + c(1,1)
      }
    }
  }
  
  g <- graph(edges=edge_list, n=M, directed = FALSE)
  V(g)$label = 1:M
  
  E(g)$capacity <- 1
  #E(g)$weight <- 100 * capacity_list / sum(capacity_list)
  #E(g)$width <- 10*E(g)$weight/max(E(g)$weight)+1
  # This is used solely for plotting purposes. To balance edge thickness, their capacities are normalized and scaled.
  #V(g)$x <- all_squares[1,]
  #V(g)$y <- all_squares[2,]
  #V(g)$label <- str_squares
  #V(g)$sample <- as.numeric(t(sample_discretized))
  # Assigning coordinates and labels to all vertices but those removed due to not containing any observation.

  adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
  # Since it is known that adj_mat is symmetric, this line copies all weights from the upper triangular part of the matrix
  # to the bottom one, ensuring symmetry.
  
  degreevec <- sapply(1:M, function(i)sum(adj_mat[i,]))
  V(g)$degree <- degreevec
  V(g)$local_maximum <- sapply(1:M, function(i)max((adj_mat[i,] > 0) * degreevec) <= degreevec[i])
  #for(i in 1:M) V(g)$local_maximum[i] <- ifelse(all(adj_mat[i,]==0 | degreevec[i] >= degreevec), degreevec[i], 0)
  #print(V(g)$local_maximum)
  # Determines whether each vertex i is a local maximum, i.e. if deg(i) = max{deg(j) | j~i}.
  
  return(list(graph=g, weight_mat=adj_mat, removed_vertices=c()))
}

g1 <- generate_sbm_graph(c(10,10), 0.5, 0.1)
#g1$weight_mat
a <- xvstlocmax_merge(g1)
length(a$ncut)

algs_gauss <- tidyr::pivot_longer(data.frame(cbind(sapply(1:length(distrange), function(i)median(vmmat[i,2*BigN+1:BigN],na.rm=T)),
                                                   sapply(1:length(distrange), function(i)median(vmmat2[i,2*BigN+1:BigN],na.rm=T)),
                                                   sapply(1:length(distrange), function(i)median(vmmat3[i,2*BigN+1:BigN],na.rm=T)),
                                                   sapply(1:length(distrange), function(i)median(vmmat4[i,2*BigN+1:BigN]))), row.names=NULL),
                                  cols=1:2, names_to="algs", values_to="values")
algs_gauss_new <- cbind(algs_gauss, x=rep(distrange, each=2))
library(scales)
gg3 <- ggplot(algs_gauss_new, aes(x=x, y=values, colour=algs)) +
  geom_line(aes(x=x, y=values, color=algs), size=1.2) +
  #geom_vline(aes(xintercept=x), linetype=2, color="grey") +
  #scale_y_continuous(trans = scales::trans_new(name="almost_log", function(x)log10(x+1e-5), function(x)10^x-1e-5), breaks=c(1e-5,1e-4,1e-3)) + #ylim(c(0,0.0015)) +
  #geom_hline(yintercept = 0.042947, colour="blue", linetype="dashed", show.legend=FALSE) +
  #geom_point() + geom_point(aes(x=xend, y=values), shape=1) + geom_segment(key_glyph="path") + coord_cartesian(xlim=c(3,6)) +
  labs(x="Distance between Gaussians", y="Median NCut") + theme_minimal() +
  scale_color_manual(name="Algorithm", values = c("green4", "lawngreen"),
                     labels=c("Xvst", "Spec. clust."),
                     guide = guide_legend(override.aes = list(linetype = rep("solid", 2), shape = rep(NA, 2)))) +
  #theme(legend.position = "bottom", legend.box.background = element_rect(colour="black"))
  #scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),
  #                   labels = trans_format("log2", math_format(2^.x))) +
  #ylim(c(0,7e-6))+
  theme(legend.position = c(0.8,0.6), legend.box.background = element_rect(colour = "black"))
gg3

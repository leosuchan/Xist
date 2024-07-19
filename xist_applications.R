############################################################################################################################
# R code supplementary to the paper "A scalable clustering algorithm to approximate graph cuts" by Suchan, Li, Munk (2023) #
############################################################################################################################
# Part II: Applications and examples


# Multiway clustering comparison between Xist and spectral clustering -------------


M <- 128
# Defining grid size M; this can be either integer or 2d integer vector of grid lengths the image (below) will be discretized into.

k <- 10
# Defining the number of partitions k the raster image will be divided into.

setwd("/path/to/NIH3T3/data")
# Sets the working directory to the folder containing the cell dataset.

confocal <- raster::raster("./s_C001Tubulin13.tif")
# Uses the "raster" package to read the .tif file into a raster object. Please ignore the "unknown extent" warning.

rasterVis::levelplot(flip(t(confocal), 2), margin = FALSE, colorkey = FALSE, par.settings = list(axis.line=list(col='transparent')),
                     scales = list(draw=FALSE), col.regions = colorRampPalette(RColorBrewer::brewer.pal(9, 'Greys')), at = seq(0,5000,length.out=101))
# Plotting the cell sample with proper scaling and b/w coloring. If colors are desired use e.g. 'YlOrRd' instead of 'Greys'.

confocal_sample <- as.matrix(raster::aggregate(confocal, fact=dim(confocal)[1:2]/M, fun=sum, expand=FALSE)) %>% (function(a)a/sum(a))
# Discretizes confocal to a matrix of dimensions M and normalizes it. Again, please ignore the warning.

g_confocal <- generate_grid_graph(confocal_sample, t=1.5)
# Generates a graph out of confocal_sample.

confocal_multixist <- multixist(g, k=k)
# Applies Multi-Xist to the image graph to obtain a list of k partitions representing the output of Xist after each of the 9 partitions.

colvec <- c("orangered1", "lawngreen", "yellow", "purple","green4", "darkgoldenrod2", "blue", "cyan", "brown3", "magenta2")
# Specifies the color vector.

plot_graph_partition(g_confocal, colvec[confocal_multixist[[k]]])
# Plots the final partition Multi-Xist output using colvec.

confocal_spec <- manual_multiway_spec_clustering(g, k=k)
# Computes the multiway spectral clustering partition consisting of k=10 parts.

plot_partition_overlay(confocal, matrix(confocal_multixist[[k]], M), colvec)
plot_partition_overlay(confocal, matrix(confocal_spec[[k-1]], M), colvec)
# Compares the partitions of Xist with those of spectral clustering, depicted as an overlay over the original image cf1_13.
# Note that multixist() returns a list, with the first component being the trivial partition 1,...,1.
# manual_multiway_spec_clustering doesn't do this, so the last index of the list here is actually k-1.


# Application example on a large network dataset ----------------


setwd("/path/to/SNAP/data")
# Sets the working directory to the folder containing the SNAP datasets.

sample_hep <- read.csv("musae_squirrel_edges.csv") + 1
# Reads in the dataset, adds 1 to all indices to let them start at 1 instead of 0 which igraph::graph() doesn't like.

g_hep <- simplify(largest_component(graph(c(t(sample_hep)), directed=FALSE)))
# Builds a graph out of the largest connected component only (to remove obvious cuts).

E(g_hep)$weight <- 1
V(g_hep)$degree <- sapply(1:length(V(g_hep)), function(i)sum(E(g_hep)[.inc(i)]$weight))
# Makes all edges unweighted and computes the vertex degree. Note that this is rather inefficient and may take a while for large graphs.

st <- Sys.time()
hep_xist <- xist(g_hep, locmax_by_degree=TRUE)
et <- Sys.time()
print(paste("Xist returned an NCut value of", hep_xist$ncut_min, "in", round(difftime(et, st, units="secs")[[1]],2), "seconds."))
# Computes and measures the time of Xist run on the graph built from the dataset.

resseq <- (1e-7)*1.1^(1:100)
resleiden <- which.min(sapply(resseq, function(res)xcut_multipartition_S(g_hep,cluster_leiden(g_hep, resolution_parameter=res)$membership)$ncut_min))
# Tries to find the optimal (i.e. minimal NCut) parameter for Leiden. Uses xcut_multipartition_S to compute XCut values of k-fold partitions.

st <- Sys.time()
hep_leiden <- xcut_multipartition_S(g_hep, cluster_leiden(g_hep, resolution_parameter=resseq[resleiden])$membership)
et <- Sys.time()
print(paste("The Leiden oracle returned an NCut value of", hep_leiden$ncut_min, "in", round(difftime(et, st, units="secs")[[1]],2), "seconds."))
# Computes and measures the time of the Leiden oracle. Note that the time taken to find the optimal resolution parameter is not measured.


# Runtime comparison of Xvst and Xist with other algorithms ---------------------


setwd("/path/to/NIH3T3/data")
# Sets the working directory to the folder containing the cell dataset.

Mvec <- c(8,9,12,14,18,21,24,28)
# Resolution levels for the confocal images, so that the discretized images will be squares containing Mvec^2 vertices.

timearray <- data.frame(n = rep(Mvec^2,each=21*5), value = rep(NA,length(Mvec)*21*5), alg = rep(paste0("alg",1:5),length(Mvec)*21))
# Initializes the dataframe containing the empirical runtimes of all algorithms, for each image and at every resolution level.

progress <- txtProgressBar(min = 0, max = length(Mvec)*21, style = 3)
# Initializes progess bar.

for(i in 1:length(Mvec)){
  
  for(j in c(1:21)){
  # Iterates over all resolution levels (i through Mvec) and confocal images (j through the .tif images).
    
    cfr1 <- suppressWarnings(raster::raster(paste0("./s_C001Tubulin", ifelse(j<17,j,j+1), ".tif")))
    cfr <- suppressWarnings(crop(cfr1, extent(cfr1, 9, 512, 1, 504)))
    # Raster and crop a confocal image to a 504 x 504 image. "Unknown extent" warnings are suppressed here; I'm not sure why they appear.
    
    cfr_sample <- as.matrix(raster::aggregate(cfr, fact=dim(cfr)[1:2]/Mvec[i], fun=sum, expand=FALSE))
    cfr_sample <- cfr_sample / sum(cfr_sample)
    # Discretizes confocal to a matrix of dimensions M and normalizes it.
    
    gc <- generate_grid_graph(cfr_sample, t = 1.5)
    gc <- delete.vertices(gc, which(V(gc)$sample == 0))
    # Removes vertices without edges (i.e. all edge weights are 0) to prevent degenerated eigenvalues.
    
    wm_gc <- matrix(rep(0,length(V(gc))^2), length(V(gc)))
    wm_gc[ends(gc, E(gc))] <- E(gc)$weight
    wm_gc <- wm_gc + t(wm_gc)
    # Computes the weight matrix wm_gc of gc for use in sClust::VonLuxburgSC.
    
    timevec <- Sys.time()
    test <- xcut_via_stmincut(gc, progress_bar = FALSE)
    timevec[2] <- Sys.time()
    test <- xist(gc, progress_bar = FALSE)
    timevec[3] <- Sys.time()
    test <- dbscan::dbscan(wm_gc, eps = 1)
    timevec[4] <- Sys.time()
    test <- sClust::VonLuxburgSC(wm_gc, K = 2) # Alternative: test <- manual_multiway_spec_clustering(gc, k=2)
    timevec[5] <- Sys.time()
    test <- cluster_leiden(gc, resolution_parameter = 1e-6)
    timevec[6] <- Sys.time()
    # Calls each algorithm separately and saves the time between each call.
    
    timearray$value[1:5 + (j-1)*5 + (i-1)*5*21] <- sapply(1:5, function(i)difftime(timevec[i+1], timevec[i], units="secs")[[1]])
    # Saves the time differences of the respective algorithms in timearray.
    
    setTxtProgressBar(progress, j + (i-1)*21)
    # Updates the progress bar.

  }
}
# Computes the empirical runtime of five clustering algorithms on 21 discretized confocal cell images at different resolution levels M.

#save(timearray, "timearray_5algs.RData")
#load("timearray_5algs.RData")
# Saves/Loads timearray for later use.

lmcoeff <- sapply(1:5, function(i)lm(log2(value) ~ log2(n), timearray, timearray$alg==paste0("alg",i))$coefficients)
lmcoeff_df <- data.frame(slope = round(lmcoeff[2,], 2), intercept = lmcoeff[1,])
# Computes the coefficients necessary to plot the lines indicating empirical runtime for each algorithm.

ggplot(timearray, aes(x=n, y=value, group=interaction(n,alg), fill=alg)) + geom_boxplot(position="identity") +
  scale_x_continuous(trans='log2', breaks=Mvec^2, labels=trans_format("sqrt", math_format(.x^2))) + theme_minimal() +
  scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x="Number of vertices n", y="Running time (sec.)") + coord_cartesian(ylim=c(2^(-10),2^11)) +
  geom_abline(aes(slope=slope, intercept=intercept, color=factor(1:5)), lmcoeff_df, linetype="dashed", size=.8) +
  scale_fill_manual(name="Algorithm", values=c("blue", "orangered1", "orange1", "yellow2", "mediumpurple2"),
                    labels=c("Xvst", "Xist", "Spec. Clust.", "Leiden", "DBSCAN")) +
  scale_color_manual(name="Runtime", values=c("blue", "orangered1", "orange1", "yellow2", "mediumpurple2"),
                     labels=as.expression(sapply(1:5, function(i)bquote(n^.(lmcoeff_df[bquote(.(i)),1]))))) +
  theme(legend.text.align = 0, legend.position = c(0.2,0.81), legend.box="horizontal",
        legend.box.background = element_rect(color = "black",fill="white"), legend.text=element_text(size=14),
        text = element_text(size=rel(4)), axis.text=element_text(size=14), panel.background = element_rect(fill="gray98", color="white")) +
  guides(fill=guide_legend(order=1), color=guide_legend(order=2))
# Plots the contents of timearray as boxplots at each resolution level, with added runtime lines for comparison.


# Algorithm comparison on a weighted sum of Gaussians -----------------


n <- 100
BigN <- 100
distrange <- seq(2,6,0.1)
# Sets sample size n, number of samples BigN, and vector of distance between Gaussians distrange.

resseq <- (1e-8)*1.1^(1:150)
epsseq <- seq(0.02,5,0.02)
# Sets sequence of resolution parameters for Leiden (resseq) and DBSCAN (epsseq) to optimize over.

mat_rate <- mat_ncut <- matrix(nrow=length(distrange), ncol=BigN) |> (function(a)list(xist=a, specclust=a, leiden=a, dbscan=a))()
# Initializing the matrix lists for the classification rate (mat_rate) and NCut values (mat_ncut) for all four algorithms.

progress <- txtProgressBar(min = 0, max = length(distrange)*BigN, style = 3)
# Initializing the progress bar.

for(dst in 1:length(distrange)){
  
  for(i in 1:BigN){
  # Iterates over all distances dst and, for each distance, takes BigN samples.
    
    n1 <- rbinom(1, n, 1/2)
    gauss_mix_sample <- rbind(matrix(rnorm(2*n1), ncol=2), matrix(c(rnorm(n-n1, mean=distrange[dst]), rnorm(n-n1)), ncol=2))
    # Determines a sample out of n1 standard 2d Gaussians and n-n1 shifted standard Gaussians, where n1 is Bernoulli random variable.
    
    gauss_mix_clusters <- c(rep(1,n1),rep(2,n-n1))
    # Represents the underlying true cluster assignment.
    
    gauss_graph <- generate_sample_graph(gauss_mix_sample, t=.2, sigm=.2, nearest_neighbors=5)
    # Computes a 0.2-neighbourhood graph with additional 5-nearest-neighbourhood and weights with exponential decay (with scaling parameter 0.2)
    
    gauss_cut <- xist(gauss_graph, locmax_by_degree=TRUE, progress_bar=FALSE)
    gauss_spec <- manual_multiway_spec_clustering(gauss_graph, k=2)[[1]]
    # Applies Xist and spectral clustering to gauss_graph.
    
    gauss_cut_trafo <- rep(1,n)
    gauss_cut_trafo[gauss_cut$ncut] <- 2
    # Transforms the output of Xist from vertex indices to cluster assignments for all vertices..
    
    mat_rate$xist[dst, i] <- classification_rate(gauss_cut_trafo, gauss_mix_clusters)
    mat_rate$specclust[dst, i] <- classification_rate(gauss_spec, gauss_mix_clusters)
    mat_rate$leiden[dst, i] <- max(sapply(resseq, function(res)classification_rate(cluster_leiden(gauss_graph, resolution_parameter=res)$membership, gauss_mix_clusters)))
    mat_rate$dbscan[dst, i] <- max(sapply(epsseq, function(eps)classification_rate(dbscan::dbscan(gauss_mix_sample, eps=eps)$cluster, gauss_mix_clusters)))
    # Computes and stores the classification rate of the four algorithms. Leiden and DBSCAN are oracles (i.e. selecting the best rate over resseq/epsseq).
    
    mat_ncut$xist[dst, i] <- gauss_cut$ncut_min
    mat_ncut$specclust[dst, i] <- xcut_partition_S(gauss_graph, which(gauss_spec==1))$ncut_min
    mat_ncut$leiden[dst, i] <- min(sapply(resseq, function(res)xcut_multipartition_S(gauss_graph, cluster_leiden(gauss_graph, resolution_parameter=res)$membership)$ncut_min))
    mat_ncut$dbscan[dst, i] <- min(sapply(epsseq, function(eps)xcut_multipartition_S(gauss_graph, dbscan::dbscan(gauss_mix_sample, eps=eps)$cluster)$ncut_min))
    # Computes and stores the XCut value of the four algorithms. Leiden and DBSCAN are oracles (i.e. selecting the best XCut over resseq/epsseq).
    
    setTxtProgressBar(progress, i+(dst-1)*BigN)
    # Updates the progress bar.
    
  }
  
}

close(progress)
# Closes the progress bar after the matrices are filled.

#save(mat_rate, mat_ncut, file="4algs_gaussians_k5_t0p2_sigm0p2_n100_BigN100.RData")
#load("4algs_gaussians_k5_t0p2_sigm0p2_n100_BigN100.RData")
# Saves/Loads mat_rate and mat_ncut for later use.

plot_graph_partition(gauss_graph, gauss_cut, vertex_scale=1, vertex_offset=1, edge_scale=.5, edge_offset=.5)
# Plots the graph and Xist partition from the final cut, just to show an example.

gauss_rate_trafo <- data.frame(x=rep(distrange,4), values=unlist(lapply(mat_rate,rowMeans)), alg=rep(paste0("alg",1:4),each=length(distrange)))
gauss_ncut_trafo <- data.frame(x=rep(distrange,4), values=unlist(lapply(mat_ncut,rowMeans)), alg=rep(paste0("alg",1:4),each=length(distrange)))
# Transforms the classification rate and NCut matrices for plotting.

ggplot(gauss_rate_trafo, aes(x=x, y=values, color=alg, group=alg)) +
  geom_line(linewidth=1.5) + labs(x=expression("Distance"~delta), y="Classification rate") + theme_minimal() +
  theme(text = element_text(size=rel(4)), axis.text=element_text(size=14), panel.background = element_rect(fill="gray98", color="white")) +
  scale_color_manual(name="Algorithm", values = c("orangered1", "orange1", "yellow2", "mediumpurple2"),
                     labels=c("Xist", "Spectral clustering", "Leiden (oracle)", "DBSCAN (oracle)"),
                     guide = guide_legend(override.aes = list(linetype = rep("solid", 4), shape = rep(NA, 4)))) +
  theme(legend.position = c(0.8,0.2), legend.box.background = element_rect(color = "black"), legend.text=element_text(size=14))
# Plots the classification rate of the four algorithms over the distance delta between the Gaussians.

ggplot(gauss_ncut_trafo, aes(x=x, y=values, color=alg, group=alg)) +
  geom_line(linewidth=1.5) + labs(x=expression("Distance"~delta), y="NCut value") + theme_minimal() +
  theme(text = element_text(size=rel(4)), axis.text=element_text(size=14), panel.background = element_rect(fill="gray98", color="white")) +
  scale_color_manual(name="Algorithm", values = c("orangered1", "orange1", "yellow2", "mediumpurple2"),
                     labels=c("Xist", "Spectral clustering", "Leiden (oracle)", "DBSCAN (oracle)"),
                     guide = guide_legend(override.aes = list(linetype = rep("solid", 4), shape = rep(NA, 4)))) +
  scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  theme(legend.position = c(0.25,0.2), legend.box.background = element_rect(color = "black"), legend.text=element_text(size=14))
# Plots the (log-scaled) NCut value of the four algorithms over the distance delta between the Gaussians.

###########################################################################################################################
# R code supplementary to the paper "Distributional limits of graph cuts on discretized grids" by Suchan, Li, Munk (2024) #
###########################################################################################################################
# Part II: Applications and examples


source("graph_cut_limits.R")
# Loads the necessary functions. Requires the "graph_cut_limits.R" file from https://github.com/leosuchan/Xist.


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

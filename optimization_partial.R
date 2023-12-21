library(plyr)
library(dplyr)
library(tidyr)
library(bootnet)
library(ggplot2)
library(ComplexHeatmap)


#' compute partial correlation networks on single tissue
#'
#' @param tissue
#' @param gamma tuning parameter
#' @param traits traits to add to the network
#' @param pval_thr indicates pval of minimum significance
#' @param subset_trait boolean that indicates if to filter by significance to traits
#' @param thresholded boolean that indicates if to threshold the network
#'
#' @return save the network and the adjacency matrix
part.cor.net.tissue <- function(tissues, gamma=0.5, traits=c("cogng_demog_slope", "motor10_lv", "tdp_st4", "amylsqrt", "tangsqrt"),
                                pval_thr=0.05, subset_trait=F, thresholded=F, penalize=F, boot=F, save=F){
  Results <- compute_network(tissues, traits=traits, pval_thr=pval_thr, subset_trait = subset_trait, gamma=gamma, threshold = thresholded, penalize=penalize, boot=boot)
  Network <- Results$Network
  adj_matrix <- as.matrix(Network$graph)
  tissues <- sort(tissues)
  tissue_path <- ""
  for (tissue in tissues){
    tissue_path <- paste0(tissue_path, substr(tissue,1,1))
  } 
  subset_path <- ifelse(subset_trait, "subset_net/", "full_net/")
  thr_path <- ifelse(thresholded, "thresholded/", "")
  penalize_path <- ifelse(penalize, "penalized/", "")
  boot_path <- ifelse(boot, "boot/", "")
  traits <- sort(traits)
  trait_path <- ""
  for (trait in traits){
    if (trait_path == ""){
      trait_path <- "_"
    }
    trait_path <- paste0(trait_path, substr(trait, 1, 1))
  }
  heatmap_obj <- Heatmap(adj_matrix,
                         name = "Partial Correlation",
                         show_row_names = TRUE,
                         cluster_rows = F,
                         cluster_columns = F)
  if (save){
    dir.create("plots/partial_cor/", showWarnings = FALSE)
    dir.create(paste0("plots/partial_cor/",subset_path), showWarnings = FALSE)
    dir.create(paste0("plots/partial_cor/", subset_path, thr_path), showWarnings = FALSE)
    dir.create(paste0("plots/partial_cor/", subset_path, thr_path, penalize_path), showWarnings = FALSE)
    dir.create(paste0("plots/partial_cor/", subset_path, thr_path, penalize_path,boot_path), showWarnings = FALSE)
    dir.create(paste0("plots/partial_cor/", subset_path, thr_path, penalize_path, boot_path, pval_thr, "/"), showWarnings = FALSE)
    
    pdf(file = paste0("plots/partial_cor/",subset_path, thr_path, penalize_path, boot_path, pval_thr, "/", tissue_path, trait_path, "_", pval_thr, ".pdf"))
    plot(heatmap_obj)
    plot(Network)
    dev.off()
  }
  else {
    plot(heatmap_obj)
    plot(Network)
  }
}

tissues <- c("Muscle", "SpinalCord")
for (tissue in tissues){
  part.cor.net.tissue(tissue, subset_trait=T, thresholded=T)
}


#' Compute a curves plot that show the influence of each composite of the EBIC score 
#' for different values of lambda
#'
#' @param tissue 
#' @param traits to add to the network 
#' @param pval_thr indicates pval of minimum significance
#' @param subset_trait boolean that indicates if to filter by significance to traits
#' @param gamma tuning parameter
#'
#' @return save the curves plot
ebic_partition <- function(tissue, traits= c("cogng_demog_slope", "motor10_lv", "tdp_st4", "amylsqrt", "tangsqrt"), pval_thr=0.1, subset_trait=T, gamma=1){
  Results <- compute_network(tissue, traits=traits, pval_thr=0.1, subset_trait = subset_trait, gamma=gamma)
  lambda_values <- c()
  Network <- Results$Network
  bar_data <- matrix(0, nrow = 10, ncol = 3)
  
  # Loop through each row and compute the values
  for (i in 1:10) {
    # Replace these calculations with your actual computations
    #sapply(seq_along(lambda), function(i) {}
    j <- 1 + (i-1)*10
    lambda_values <- c(lambda_values,round(Network$results$results$rholist[j],4))
    S <-  bootnet:::bootnet_correlate(data = Network$data, corMethod = Network$arguments$corMethod,
                                      corArgs = list(), verbose = T,
                                      nonPositiveDefinite = c("stop", "continue"))
    n <- Network$nPerson
    gamma <- Network$arguments$tuning
    invSigma <- Network$results$results$wi[,,j]
    L <- qgraph:::logGaus(S, invSigma, n)
    E <- sum(invSigma[lower.tri(invSigma, diag = F)]!= 0)
    p <- Network$nNode
    # Assign the computed values to the corresponding row in the result matrix
    bar_data[i, ] <- c(-2*L, E*log(n) + 4*E*gamma*log(p), -2*L + E*log(n) + 4*E*gamma*log(p))
  }
  legend_labels <- c("BIC", "gamma penalization", "EBIC")
  curve_df <- data.frame(x = lambda_values, bar_data)
  curve_df <- curve_df %>%
    gather(key = "Curve", value = "y", -x)
  custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c")
  curves <- ggplot(data = curve_df, aes(x = x, y = y, color = Curve)) +
    geom_line() +
    labs(x = "lambda", color = "Courbe") +
    scale_color_manual(values = custom_colors, labels = legend_labels)  # Personnalisez les couleurs si nécessaire
  ggtitle("EBIC decomposition analysis")
  dir.create("plots/part_cor/final/ebic", showWarnings = FALSE)
  ggsave(paste0("plots/part_cor/final/ebic/", tissue ,"_ebic_partition.png"), plot = curves, width = 6, height = 6, dpi = 300)
}


#' compute curves plot that show the evolution of the EBIC score given lambda 
#' for different values of gamm
#'
#' @param tissue 
#' @param traits to add to the network
#' @param pval_thr indicates pval of minimum significance
#' @param subset_trait boolean that indicates if to filter by significance to traits
#'
#' @return save the curves plot
gamma_analysis <- function(tissue, traits, pval_thr=0.1, subset_trait=T){
  MEs <- load.data(tissue, traits, subset_trait)
  gamma_part_data <- matrix(0, nrow=10, ncol=5)
  lambda_values <- c()
  gamma_values <- c()
  for (i in 1:5){
    gamma <- 0.25 * (i-1)
    gamma_values <- c(gamma_values, gamma)
    Network <- compute_network(MEs, gamma)
    
    S <-  bootnet:::bootnet_correlate(data = MEs, corMethod = Network$arguments$corMethod,
                                      corArgs = list(), verbose = T,
                                      nonPositiveDefinite = c("stop", "continue"))
    n <- Network$nPerson
    for (j in 1:10){
      if (i==1){
        k <- 1 + (j-1)*10
        lambda_values <- c(lambda_values, round(Network$results$results$rholist[k],4))
      }
      invSigma <- Network$results$results$wi[,,j]
      L <- qgraph:::logGaus(S, invSigma, n)
      E <- sum(invSigma[lower.tri(invSigma, diag = F)]!= 0)
      p <- Network$nNode
      gamma_part_data[j,i] <- -2*L + E*log(n) + 4*E*gamma*log(p)
    }
  }
  
  curve_df <- data.frame(x = lambda_values, gamma_part_data)
  curve_df <- curve_df %>%
    gather(key = "Curve", value = "y", -x)
  custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
  curves <- ggplot(data = curve_df, aes(x = x, y = y, color = Curve)) +
    geom_line() +
    labs(x = "lambda", color = "Courbe") +
    scale_color_manual(values = custom_colors, labels = gamma_values)  # Personnalisez les couleurs si nécessaire
  ggtitle("EBIC comparison over gamma values")
  dir.create("plots/part_cor/final/ebic", showWarnings = FALSE)
  ggsave(paste0("plots/part_cor/final/ebic/", tissue ,"_gamma_analysis.png"), plot = curves, width = 6, height = 6, dpi = 300)
}



#' compute scatter plot that compare the correlation values to the partial 
#' correlation values in a network
#'
#' @param tissue 
#' @param comp boolean that indicates if compare to edges involved with trait 
#' (else edges modules-modules) 
#' @param trait if comp, indicates the trait involved in the edges to analyse
#' @param pval_thr indicates pval of minimum significance
#' @param subset_trait 
#' @param gamma 
#' @param threshold 
#' @param penalize 
#'
#' @return
#' @export
#'
#' @examples
part_cor.cor.comp <- function(tissue, comp=F, trait=NULL, pval_thr=0.05, subset_trait=F, gamma=0.5, threshold=F, penalize=F){
  Results <- compute_network(tissue, traits=trait, pval_thr=pval_thr, subset_trait=subset_trait, gamma=gamma, threshold=threshold, penalize=penalize)
  group1 <- colnames(Results$data)
  if (comp){
    group2 <- trait
  } else {
    group2 <- group1
  }
  correlation_matrix <- qgraph::cor_auto(Results$Network$data)[group1, group2]
  partial_correlation_matrix <- Results$Network$graph[group1, group2]
  
  correlation_data <- data.frame(
    Node1 = rep(group1, each = length(group2)),
    Node2 = rep(group2, times = length(group1)),
    Correlation = as.vector(correlation_matrix),
    Partial_Correlation = as.vector(partial_correlation_matrix)
  )
  
  if (!comp){
    correlation_data <- correlation_data %>% filter("Node1" < "Node2")
  }
  scatterplot <- ggplot(data = correlation_data, aes(x = Correlation, y = Partial_Correlation)) +
    geom_point() +
    labs(x = "Correlation", y = "Partial Correlation") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    ggtitle(paste(tissue, trait,cor(correlation_data$Correlation, correlation_data$Partial_Correlation))) +
    lims(x = c(min(min(correlation_data$Correlation), -max(correlation_data$Correlation)),
               max(-min(correlation_data$Correlation), max(correlation_data$Correlation))),
         y = c(min(min(correlation_data$Partial_Correlation), -max(correlation_data$Partial_Correlation)),
               max(-min(correlation_data$Partial_Correlation), max(correlation_data$Partial_Correlation))))
  filename <- paste0("plots/part_cor/final/comp/",tissue, "_", trait, ".png")
  ggsave(filename)
}

tissues <- c("Brain","SpinalCord", "Muscle")
for (tissue in tissues){
  part_cor.cor.comp(tissue, trait="cogng_demog_slope", comp=T)
}

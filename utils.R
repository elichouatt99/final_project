library(plyr)
library(dplyr)
library(tidyr)
library(bootnet)
library(ggplot2)
library(ComplexHeatmap)


#' @param tissue from which tissue we need to get the MEs
#' @param level level of the SpeakEasy algorithm (as 'lx' where x is the level)
#'
#' @return Return the MES of a given pair (tissue, level)
load.MEs.l4 <- function(tissue, level='l4'){
  infos <- read.csv("datasets/samples_infos.csv")
  path <- paste0("datasets/partial_cor/", tissue, "_modules.rds")
  dataset <- readRDS(path)
  MEs <- dataset$MEs[,grepl(paste0("^", level), colnames(dataset$MEs))]
  rownames(MEs) <- plyr::mapvalues(rownames(MEs), infos$SampleID_BatchID, infos$projid)
  return(MEs)
}


#' remove the modules that have less than a certain number of genes 
#'
#' @param MEs 
#' @param tissue 
#' @param min_size minimal number of genes for a module to be kept
#'
#' @return the MEs that have more than min_size genes
filter.length.MEs <- function(MEs, tissue, min_size=30){
  filtered.MEs <- MEs
  modules <- readRDS(paste0("datasets/partial_cor/",tissue,'_modules.rds'))
  for (module in colnames(MEs)){
    if (length(modules$genes[modules$genes$module == module,]$gene) < min_size){
      filtered.MEs[[module]] <- NULL
    }
  }
  return(filtered.MEs)
}


#' Remove the modules with pval not significant for of the traits.
#' pval.signed is -log10(FDR) * sign.beta
#'
#' @param MEs 
#' @param SE dataset with signicance info for the modules
#' @param traits traits to associate to the modules
#' @param pval_thr maximal pval for significance
#'
#' @return the modules that are significantly correlated  with at least one trait
filter.trait.MEs <- function(MEs, SE, traits, pval_thr){
  filtered.MEs <- MEs
  for (module in colnames(MEs)){
    pval <- 10^(-SE$pval.signed[module, traits] * sign(SE$cor[module, traits]))
    if(all(pval > pval_thr)){
      filtered.MEs[[module]] <- NULL
      print(module)
    }
  }
  return(filtered.MEs)
}


#' Load MEs of given tissue, effectuate preprocessing and add phenotypes data.
#'
#' @param tissue 
#' @param traits traits add as phenotypes
#' @param subset_trait boolean if we need to subset given association 
#' significance with the traits
#'
#' @return filtered MEs with phenotypes data
load.data <- function(tissue, traits=NULL, subset_trait=F){
  MEs <- load.MEs.l4(tissue)
  SE <- readRDS(paste0("datasets/traits/SE_", tissue,"_l4.rds"))
  
  MEs <- filter.length.MEs(MEs, tissue)
  if (subset_trait){
    MEs <- filter.trait.MEs(MEs, SE, traits)
  }
  if (length(MEs) == 0){
    return()
  }
  modules.names <- colnames(MEs)
  phenotypes <- readRDS("datasets/traits/phenotypes.rds")
  MEs.merged <- merge(MEs,phenotypes[traits],by="row.names")
  rownames(MEs.merged) <- MEs.merged$Row.names
  MEs.merged$Row.names <- NULL
  return(MEs.merged)
}


#' Create block matrix from shape matrix and blocks coordinates
#' @param shape shape of the final matrix
#' @param blocks coordinates of the blocks
#'
#' @return block matrix
create_block_matrix <- function(shape, blocks) {
  result <- matrix(0, nrow = shape, ncol = shape)
  for (block in blocks) {
    coordinates <- block$coordinates
    value <- block$value
    result[coordinates[1]:coordinates[3], coordinates[2]:coordinates[4]] <- value
  }
  return(result)
}


#' Return clusters from correlation matrix between elements, given the final 
#' number of clusters or the cut height of the cluster tree
#'
#' @param matrix_cor correlation matrix between elements 
#' @param num_clusters final number of clusters
#' @param cut_height cut height of the cluster tree (dendogram)
#'
#' @return the clusters
get_clusters <- function(matrix_cor, num_clusters=5, cut_height=NULL){
  distance_matrix <- dist(matrix_cor)
  hierarchical_clustering <- hclust(distance_matrix, method = "complete")
  
  if (!is.null(cut_height)){
    clusters <- cutree(hierarchical_clustering, h = cut_height)
  }
  else {
    clusters <- cutree(hierarchical_clustering, k = num_clusters)
  }
  return(clusters)
}


#' compute the network for given tissues and traits
#'
#' @param tissues 
#' @param cluster modules of a given cluster
#' @param pval_thr maximal pval for significance
#' @param subset_trait boolean that indicance if it is needed to filter by significance with traits
#' @param gamma tuning parameter of partial correlation networks
#' @param threshold boolean that indicates if thresholded network required
#' @param penalize boolean that indicates if partial penalization.
#' 
#' @return @param Network 
#' @aliases @param MEs the Mean expression of the modules of the final networks
#' @aliases @param data Mean expression of the modules of the tissues (only filtered by length)
compute_network <- function(tissues,
                            traits=c("cogng_demog_slope", "motor10_lv", "tdp_st4", "amylsqrt", "tangsqrt"),
                            cluster=NULL,
                            pval_thr=0.1,
                            subset_trait=F,
                            gamma=0.5,
                            threshold=F, 
                            penalize=F){
  tissues <- sort(tissues)
  MEs <- NULL
  
  blocks <- list()
  current_index <- 1
  for (tissue in tissues){
    temp.data <- load.MEs.l4(tissue)
    temp.data <- filter.length.MEs(temp.data, tissue)
    temp.SE <- readRDS(paste0("datasets/traits/SE_", tissue,"_l4.rds"))
    if (subset_trait){
      temp.MEs <- filter.trait.MEs(temp.data, temp.SE, traits, pval_thr)
    }
    else {
      temp.MEs <- temp.data
    }
    if (!is.null(cluster)){
      temp.MEs <- temp.data[, cluster]
    }
    coordinates <- list(list(coordinates = c(current_index, current_index, current_index + length(temp.MEs) - 1, current_index + length(temp.MEs)-1), value = 1))
    blocks <- append(blocks, coordinates)
    current_index <- current_index + length(temp.MEs)
    
    if(is.null(MEs)){
      MEs <- temp.MEs
      data <- temp.data
    } 
    else {
      MEs <- merge(MEs, temp.MEs, by="row.names")
      rownames(MEs) <- MEs$Row.names
      MEs$Row.names <- NULL
      
      data <- merge(data, temp.data, by="row.names")
      rownames(data) <- data$Row.names
      data$Row.names <- NULL
    }
  }
  
  modules.names <- colnames(MEs)
  
  phenotypes <- readRDS("datasets/traits/phenotypes.rds")
  MEs.traits <- merge(MEs, phenotypes[traits], by="row.names")
  rownames(MEs.traits) <- MEs$Row.names
  MEs.traits$Row.names <- NULL
  
  coordinates <- list(list(coordinates = c(length(MEs)+1, length(MEs)+1, length(MEs.traits), length(MEs.traits)), value = 1))
  blocks <- append(blocks, coordinates)
  
  size <- length(MEs.traits)  
  if (penalize){
    penalizeMatrix <- create_block_matrix(size, blocks)
  }
  else {
    penalizeMatrix <- matrix(1, nrow = size, ncol = size)
  }
  
  
  Network <- estimateNetwork(
    MEs.traits,
    default = "EBICglasso",
    corMethod = "cor_auto",
    tuning = gamma,
    threshold = threshold,
    penalizeMatrix =  penalizeMatrix)
  
  return(list(Network=Network, MEs=MEs, data=data))
}

plot_figure <- function(tissues, 
                        filename,
                        pval_thr=0.1,
                        subset_trait=F,
                        gamma=0.5,
                        threshold=F, 
                        penalize=F,
                        clusters=NULL){
  Results <- compute_network(tissues, pval_thr, subset_trait, gamma, threshold, penalize, clusters)
  png(filename, 2048,2048)
  plot(Results$Network)
  dev.off()
}

